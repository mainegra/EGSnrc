
/*
###############################################################################
#
#  EGSnrc egs++ egs_shield application
#  Copyright (C) 2026 National Research Council Canada
#
#  This file is part of EGSnrc.
#
#  EGSnrc is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Affero General Public License as published by the
#  Free Software Foundation, either version 3 of the License, or (at your
#  option) any later version.
#
#  EGSnrc is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#  FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for
#  more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with EGSnrc. If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################
#
#  Author:        Ernesto Mainegra-Hing, 2026
#
###############################################################################
#
#  Buildup-factor code using staged forced-detection transport with
#  weight equalization (combing) for population control.
#
#  Two distinct sets of shells are decoupled:
#
#    Scoring shells  (N shells, user-specified):
#      Thin air shells where kerma is estimated via the FD estimator.
#      Photons pass through these freely in actual transport — only the
#      FD weight is accumulated.  No killing at scoring shells.
#
#    Combing shells  (M shells, M <= N, subset of scoring shell positions):
#      Population-control checkpoints.  When the FD ray-trace reaches a
#      combing shell for a scattered photon: the shell is scored, a virtual
#      photon is saved at the shell EXIT (iron side), and the actual photon
#      is killed upon entry.  Containers are snapshot-swapped and combed
#      between iterations.
#
#  Algorithm (Population-controlled batch FD):
#
#    For each of n_bunches independent bunches of n_per_bunch source photons:
#
#      Stage 0:
#        Transport source photons (latch = 0).
#        selectPhotonMFP fires at each free path start:
#          Primary   (latch == 0): scoreFD_all to all N scoring shells.
#          Scattered (latch != 0): scoreFD_all to shells up to first combing
#                                  shell; save virtual photon at shell EXIT.
#        ausgab kills actual photon upon combing-shell entry.
#
#      Iterative backscatter loop (per bunch):
#        Snapshot all containers, comb them, clear containers.
#        Replay every snapshot particle (same FD + kill logic as Stage 0).
#        Repeat until all containers are empty.
#
#      Record K_total[k] and K_pri[k] for this bunch.
#
#    Output: BUF[k]      = mean(K_tot/K_pri) across bunches.
#            sigma(BUF)  = between-bunch stddev / sqrt(n_bunches).
#
#  Reliable uncertainty: each bunch is one statistically independent event.
#
#  Input block:
#
#    :start scoring options:
#        geometry name     = iron_sphere
#        scoring regions   = 1 3 5 7 9 ...   # 0-based, thin air scoring shells
#        scoring volumes   = 0.5 0.5 ...     # cm^3, one per scoring shell
#        combing regions   = 41 49 ...       # thin iron regions between scoring positions (NOT scoring shells)
#        emuen file        = $EGS_HOME/egs_kerma/emuen_rho_air_1keV-20MeV.data
#        n bunches         = 100
#        photons per bunch = 1000000
#    :stop scoring options:
#
###############################################################################
*/

#include "egs_advanced_application.h"
#include "egs_interface2.h"
#include "egs_functions.h"
#include "egs_input.h"
#include "egs_base_source.h"
#include "egs_rndm.h"
#include "egs_interpolator.h"
#include "egs_timer.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <string>

using std::string;
using std::vector;

#define calculatePhotonMFP \
    F77_OBJ_(calculate_photon_mfp,CALCULATE_PHOTON_MFP)
extern __extc__ void calculatePhotonMFP(EGS_Float *, EGS_Float *);

static const EGS_Float TSTEP_MAX = 1e35;


/*============================================================================
  TmpPhsp — growable in-memory particle buffer with weight equalization.
============================================================================*/
static const int MAXBUF = 40000000;

class TmpPhsp {
public:
    TmpPhsp() : ntot(4096), np(0), p(new EGS_Particle[4096]) {}
    ~TmpPhsp() { delete [] p; }

    void save(const EGS_Particle &par) {
        if (np >= ntot) grow();
        p[np++] = par;
    }

    // Weight equalization: w_bar = sum(wi)/n_target (default: current np).
    // Each particle receives floor(wi/w_bar) deterministic copies plus one
    // Russian-roulette copy for the fractional remainder.
    // n_target < 0 means use current np (no population change expected).
    // n_target > 0 caps output to that expected count (population control).
    void comb(EGS_RandomGenerator *rndm, int n_target = -1) {
        if (np < 2) return;
        double wtot = 0.0;
        for (int j = 0; j < np; j++) wtot += p[j].wt;
        if (wtot <= 0.0) { np = 0; return; }
        if (n_target <= 0) n_target = np;
        EGS_Float wbar = (EGS_Float)(wtot / n_target);

        int          tmp_ntot = ntot;
        EGS_Particle *tmp     = new EGS_Particle[tmp_ntot];
        int           nout    = 0;

        for (int j = 0; j < np; j++) {
            double ratio   = p[j].wt / wbar;
            int    ncopies = (int)ratio;
            if (rndm->getUniform() < (ratio - ncopies)) ncopies++;
            for (int k = 0; k < ncopies; k++) {
                if (nout >= tmp_ntot) {
                    int nn = 2 * tmp_ntot;
                    if (nn > MAXBUF)
                        egsFatal("TmpPhsp::comb(): temporary buffer exceeded "
                                 "MAXBUF\n");
                    EGS_Particle *tmp2 = new EGS_Particle[nn];
                    for (int m = 0; m < nout; m++) tmp2[m] = tmp[m];
                    delete [] tmp; tmp = tmp2; tmp_ntot = nn;
                }
                tmp[nout]    = p[j];
                tmp[nout].wt = wbar;
                nout++;
            }
        }
        if (nout > ntot) {
            delete [] p;
            p = new EGS_Particle[nout];
            ntot = nout;
        }
        for (int j = 0; j < nout; j++) p[j] = tmp[j];
        np = nout;
        delete [] tmp;
    }

    void          clean()           { np = 0; }
    int           size()  const     { return np; }
    EGS_Particle &operator[](int j) { return p[j]; }
    double        weight() const {
        double w = 0.0;
        for (int j = 0; j < np; j++) w += p[j].wt;
        return w;
    }

private:
    int          ntot, np;
    EGS_Particle *p;

    void grow() {
        int nn = 2 * ntot;
        if (nn > MAXBUF) egsFatal("TmpPhsp::grow(): exceeded MAXBUF\n");
        EGS_Particle *q = new EGS_Particle[nn];
        for (int j = 0; j < np; j++) q[j] = p[j];
        delete [] p; p = q; ntot = nn;
    }
};


/*============================================================================
  EGS_ShieldApplication
============================================================================*/
class APP_EXPORT EGS_ShieldApplication : public EGS_AdvancedApplication {
public:

    EGS_ShieldApplication(int argc, char **argv)
        : EGS_AdvancedApplication(argc, argv),
          nreg(0),
          n_scoring(0), scoring_reg(nullptr), scoring_idx(nullptr),
          is_scoring(nullptr), V_shell(nullptr),
          n_combing(0), combing_reg(nullptr), combing_idx(nullptr),
          is_combing(nullptr),
          containers(nullptr), snapshots(nullptr),
          E_Muen_Rho(nullptr),
          cascade_diag(false), bunch_stats(false),
          cascade_cutoff(1e-10), comb_target_in(0),
          n_bunches(100), n_per_bunch(1000000LL),
          current_bunch(0), n_completed(0),
          sum_ratio(nullptr), sum_ratio2(nullptr),
          sum_Kt(nullptr), sum_Kp(nullptr), sum_Kp0(nullptr),
          max_ratio(nullptr),
          n_valid_b(nullptr),
          total_cpu_time_(0.0),
          n_jobs_comb(0),
          fc_analog_w(0), fc_analog_n(0),
          forced_collision(false), fc_found(false),
          fc_Lambda(0), fc_ireg(-1), fc_m(-1),
          pri_diag(false), deepest_pri_b(-1)
    {}

    ~EGS_ShieldApplication();

    /*------------------------------------------------------------------------
      ausgab

      (a) Latch bookkeeping: after each scatter/absorption event, increment
          latch for all produced particles so that latch==0 identifies a
          primary photon throughout its history.
      (b) Kill photons that enter a combing shell.  The FD contribution for
          this free path was already scored inside selectPhotonMFP.
    ------------------------------------------------------------------------*/
    /* latch layout.  Bits 0..29 count interactions, so a photon that has not
       yet scattered has them all clear.  Bit 30 is NO_FD: "this photon's
       primary FD contribution is already booked".  It is set on the source
       photon straight after its single full-range scoreFD_all(true), and is
       inherited by every forced-collision descendant -- which is what keeps
       K_pri a single deterministic ray from the source however many segments
       the primary is subsequently staged through.

       Every "is this still a primary" test must mask bit 30 off.  Using
       latch != 0 directly would classify a staged primary as scattered and
       kill it at the next combing surface.  ++latch at an interaction leaves
       bit 30 untouched (the counter cannot reach 2^30), so the flag survives
       scattering -- harmlessly, since a scattered photon is never asked. */
    static const int NO_FD_FLAG = (1 << 30);
    static bool isPrimary(int latch) {
        return (latch & ~NO_FD_FLAG) == 0;
    }

    int ausgab(int iarg) override {
        int np = the_stack->np - 1;

        if (iarg == AfterCompton  || iarg == AfterPhoto  ||
            iarg == AfterPair     || iarg == AfterRayleigh) {
            int npold = the_stack->npold - 1;
            for (int ip = npold; ip <= np; ip++)
                ++the_stack->latch[ip];
            return 0;
        }

        if (iarg != BeforeTransport) return 0;
        if (the_stack->iq[np] != 0)  return 0;   // only photons

        int ir    = the_stack->ir[np] - 2;
        int latch = the_stack->latch[np];
        // Kill scattered photons at combing shells; primaries pass through freely.
        // If we killed primaries here, they could not scatter in the outer iron
        // and their secondaries (heading inward) would never contribute to the
        // combing-shell scoring — causing a systematic downward bias in BUF.
        //
        // With 'forced collision = yes' a primary can no longer get here at all:
        // its collided branch is forced to interact before the surface and its
        // uncollided branch is banked AT the surface, so this exemption becomes
        // unreachable rather than merely unused.  It is kept because the option
        // defaults off.
        if (ir >= 0 && ir < nreg && is_combing[ir] && !isPrimary(latch)) {
            the_stack->wt[np]    = 0;
            the_epcont->idisc = -1;
        }
        // Diagnostic: a primary reaching a combing surface is the one particle
        // that is neither banked nor killed there.  Count them, and record how
        // deep any primary got in this bunch.  Under forced collision this
        // must read zero, which is how you confirm the option took effect.
        else if (pri_diag && ir >= 0 && ir < nreg && is_combing[ir]
                 && isPrimary(latch)) {
            int m = combing_idx[ir];
            pri_cross_b[m]++;
            pri_cross_wsum[m] += the_stack->wt[np];
            if (m > deepest_pri_b) deepest_pri_b = m;
        }
        return 0;
    }

    /*------------------------------------------------------------------------
      selectPhotonMFP — called by Fortran at the start of each photon free
      path via the $SELECT-PHOTON-MFP macro override in egs_shield.macros.
    ------------------------------------------------------------------------*/
    void selectPhotonMFP(EGS_Float &dpmfp) {
        int  np      = the_stack->np - 1;
        int  latch   = the_stack->latch[np];
        bool primary = isPrimary(latch);

        // FD scoring.  A primary scores once, on its first free path, over the
        // full ray to the geometry boundary; NO_FD then suppresses every later
        // call for that lineage.  Scattered photons score each free path and
        // stop at the first combing surface, as before.
        bool score_as_primary = primary && !(latch & NO_FD_FLAG);
        if (primary) {
            // Always TRACE, even when the scoring is suppressed.  The trace is
            // what fills the fc_* record, and a replayed primary carries NO_FD
            // precisely so that it does not re-score -- if that also skipped
            // the trace, fc_found and fc_x/fc_Lambda/fc_m would still hold the
            // PREVIOUS particle's values, and this one would bank a copy at a
            // foreign exit point and force its collision against a foreign
            // Lambda.  Measured as a coherent positive bias in BUF growing with
            // depth (+0.03% at 0.5 mfp to +0.09% at 2.5 mfp) before the trace
            // and the scoring were separated.
            //
            // A trace-only call returns at the first combing surface, so it
            // costs a fraction of a full scoring pass.
            scoreFD_all(true, score_as_primary);
            if (score_as_primary) {
                the_stack->latch[np] = latch | NO_FD_FLAG;
            }
        }
        else {
            scoreFD_all(false);
        }

        // ---- Forced collision between combing surfaces ---------------------
        //
        // Analog transport lets a primary reach a combing surface with its full
        // weight and probability exp(-Lambda) -- at 10 mfp spacing, ~4 arrivals
        // per 10^6-photon bunch carrying weight 1 into a container whose
        // population sits at wbar ~ 10^-9.  comb() then splits that one particle
        // toward n_target copies and roulettes everything else away, so a bunch's
        // deep shells become clones of a single ancestor.  Measured: 4.0% of the
        // weight arriving at the first surface, in ~1 particle out of ~10^5.
        //
        // Forced collision replaces the sampling with the split it is estimating:
        //
        //     uncollided   w * T          banked AT the surface, T = exp(-Lambda)
        //     collided     w * (1 - T)    interacts inside the segment, at
        //                                 lambda = -ln(1 - xi(1 - T))
        //
        // The weights sum to w, so it is unbiased.  What changes is that the
        // uncollided branch now arrives with certainty at weight w*T instead of
        // with probability T at weight w -- the same expected weight, at the
        // local scale, so comb() roulettes it down instead of splitting it up.
        //
        // scoreFD_all already ray-traced this direction, so Lambda to the first
        // combing surface and the exit point come back from that trace for free.
        if (forced_collision && primary && fc_found) {
            EGS_Float T = std::exp(-fc_Lambda);
            // T == 1 would leave the collided branch weightless and the sampling
            // below degenerate; T == 0 leaves nothing to bank.  Both mean the
            // segment is not worth forcing, so fall through to analog.
            if (T > 1e-300 && T < 1.0 - 1e-9) {
                EGS_Float w = the_stack->wt[np];

                EGS_Particle vp;
                vp.q     = 0;
                vp.latch = the_stack->latch[np];   // primary, NO_FD set
                vp.E     = the_stack->E[np];
                vp.wt    = w * T;
                vp.x     = fc_x;
                vp.u     = EGS_Vector(the_stack->u[np], the_stack->v[np],
                                      the_stack->w[np]);
                vp.ir    = fc_ireg;
                containers[fc_m]->save(vp);

                fc_in_w[fc_m]   += w;
                fc_bank_w[fc_m] += w * T;
                fc_coll_w[fc_m] += w * (1.0 - T);
                fc_n_split[fc_m]++;

                the_stack->wt[np] = w * (1.0 - T);
                dpmfp = -std::log(1.0 - rndm->getUniform() * (1.0 - T));
                return;
            }
        }

        if (forced_collision && primary) {
            fc_analog_w += the_stack->wt[np];
            fc_analog_n++;
        }
        dpmfp = -std::log(1.0 - rndm->getUniform());
    }

    int  initScoring()         override;
    int  initRunControl()      override;
    int  runSimulation()       override;
    int  finishSimulation()    override;
    void outputResults()       override;
    int  outputData()          override;
    int  readData()            override;
    int  addState(istream &d)  override;
    int  combineResults()      override;
    void resetCounter()        override;

protected:
    int startNewShower() override {
        return EGS_AdvancedApplication::startNewShower();
    }

private:

    int nreg;   // total regions in geometry

    // ---- Scoring shells (N shells, FD only, no killing) ----
    int        n_scoring;
    int       *scoring_reg;   // scoring_reg[k]  = 0-based region index
    int       *scoring_idx;   // scoring_idx[ir] = k, or -1
    bool      *is_scoring;    // [nreg]
    EGS_Float *V_shell;       // shell volume [cm^3], one per scoring shell

    // ---- Combing shells (M <= N shells, subset of scoring shell positions) ----
    int   n_combing;
    int  *combing_reg;        // combing_reg[m]  = 0-based region index
    int  *combing_idx;        // combing_idx[ir] = m, or -1
    bool *is_combing;         // [nreg]

    // ---- Staged transport buffers ----
    TmpPhsp **containers;     // [n_combing]: filled during current replay pass
    TmpPhsp **snapshots;      // [n_combing]: particles being replayed this iteration

    // ---- Per-bunch kerma accumulators (reset at the end of each bunch) ----
    // Plain sums, normalized by n_per_bunch in endBunch().  Using
    // EGS_ScoringArray here was wrong: currentResult() divides by the last
    // setHistory() argument (cumulative source-photon index), not by
    // n_per_bunch, so the denominator grew with every bunch and included
    // replay-photon increments — giving K values that are ~5× too small
    // while BUF was unaffected (same wrong denominator in numerator/denominator).
    vector<double> sum_tot;
    vector<double> sum_pri;
    vector<double> sum_pri0;   // as sum_pri but with the exp(-Λ) factor omitted

    // ---- E*mu_en/rho interpolator (log E abscissa) ----
    EGS_Interpolator *E_Muen_Rho;

    // ---- Bunch bookkeeping ----
    bool      cascade_diag;    // print the replay-cascade trace (input option)
    bool      bunch_stats;     // print per-shell bunch-contribution stats (input)
    double    cascade_cutoff;  // stop cascade below this fraction of w_initial
    int       comb_target_in;  // combing population target; <=0 -> derive/preserve
    int       n_bunches;
    long long n_per_bunch;
    int       current_bunch;
    int       n_completed;     // bunches scored so far (all jobs combined)
    // Sufficient statistics for combining parallel jobs (all [n_scoring]):
    double   *sum_ratio;       // Σ (K_tot_b / K_pri_b) over valid bunches
    double   *sum_ratio2;      // Σ (K_tot_b / K_pri_b)² over valid bunches
    double   *sum_Kt;          // Σ K_tot_b
    double   *sum_Kp;          // Σ K_pri_b
    double   *max_ratio;       // largest single-bunch K_tot/K_pri, per shell.
                               // Not derivable from the other sums; shows how
                               // much one bunch dominates a shell.
    double   *sum_Kp0;         // Σ K_pri_b computed with exp(-Λ) omitted;
                               // K_pri/K_pri0 is the primary transmission, so
                               // η_eff = -ln(ΣK_pri / ΣK_pri0).  See outputResults().
    int      *n_valid_b;       // # bunches with K_pri_b > 0
    double    total_cpu_time_; // accumulated CPU time over all bunches [s]

    // ---- Between-job (batch) variance, built during combineResults() --------
    // The bunch-level sigma above and this one estimate the same quantity by
    // different batchings -- 6400 bunches vs ~800 job means.  They agree only
    // if the per-bunch distribution is sampled well enough for its variance to
    // be determined; a heavy tail makes both underestimates, and makes them
    // disagree.  This is a batch-size scaling test, not an independent
    // estimator: disagreement proves non-convergence, agreement does not prove
    // convergence.  Combine-time only, so it never enters the .egsdat and the
    // file format is unchanged.
    vector<double> job_nr2;    // Σ_j n_j * rbar_j^2  = Σ_j (ΔΣratio)^2 / Δn_j
    vector<int>    job_used;   // # jobs contributing >= 1 valid bunch, per shell
    int            n_jobs_comb;

    // ---- Top-k bunch ratios, for leave-k-out sensitivity -------------------
    // sigma and N_eff both describe the spread; neither says how much of the
    // ANSWER rests on a few bunches.  Keeping the k largest per-bunch ratios
    // lets outputResults() report how far the mean moves when the top 1, 2 ...
    // k bunches are dropped, which is the number that matters: at 40 mfp the
    // single largest of 6400 bunches carries 4.2% of the shell.
    //
    // Combines by merging sorted lists, so the top-k over all jobs is exact.
    // Every per-bunch ratio, not just the top few.  38 shells x 8 bunches is
    // 2.4 kB per job file and ~1.9 MB combined at 6400 bunches -- nothing --
    // and it is strictly more informative than any summary: sigma, N_eff,
    // top-k, batch scaling, tail index and bootstrap intervals are all
    // recoverable from it offline, and new questions can be asked of an old
    // run without re-running it.  Dumped to <output>.bunchdat on the combine.
    //
    // This matters because the deep shells turn out to have a tail index below
    // 2, i.e. infinite variance, where sigma is not a convergent quantity and
    // the honest report is an interval estimated from the sample itself.
    vector< vector<double> > bunch_r;   // [n_scoring][valid bunches]

    static const int N_TOPK = 5;
    vector<double> topk;        // [n_scoring * N_TOPK], descending per shell
    void insertTopK(int k, double v) {
        double *a = &topk[(size_t)k * N_TOPK];
        if (v <= a[N_TOPK - 1]) return;
        int i = N_TOPK - 1;
        while (i > 0 && a[i - 1] < v) { a[i] = a[i - 1]; --i; }
        a[i] = v;
    }

    // ---- Forced collision for primaries between combing surfaces -----------
    // scoreFD_all fills these on a primary trace: the first combing surface the
    // ray meets, the number of mean free paths to it, and where the ray leaves it.  selectPhotonMFP
    // then banks the uncollided branch there and forces the collided branch to
    // interact short of it.  Single-threaded, one trace at a time, so a member
    // is as safe as an out-parameter and keeps scoreFD_all's signature.
    // Weight-balance audit of the staged primary.  The split itself is exact by
    // construction (w = w*T + w*(1-T)), so that is not what this measures.  The
    // informative comparison is ACROSS stages: the weight entering the split at
    // surface m should equal the weight banked at surface m-1, because combing
    // and replay are weight-preserving in expectation.  A ratio above 1 that
    // compounds per stage localises a leak to the bank -> comb -> replay path
    // rather than to the split, which is what the 20-bunch 9-surface comparison
    // pointed at (exact before the first surface, then ~1.1-1.25 per stage).
    vector<double>    fc_in_w, fc_bank_w, fc_coll_w;
    vector<long long> fc_n_split;
    double            fc_analog_w;   // primary weight with no surface ahead
    long long         fc_analog_n;

    bool       forced_collision;   // input: 'forced collision'
    bool       fc_found;
    EGS_Float  fc_Lambda;          // number of mean free paths, source position -> surface exit
    EGS_Vector fc_x;               // exit point of the combing region
    int        fc_ireg;            // region just beyond it
    int        fc_m;               // combing-surface index

    // ---- Primary-crossing diagnostic (per job; not combined) ---------------
    // Primaries are the one species that is neither banked at a combing surface
    // nor killed there: they transport analog straight through.  Their survival
    // probability to the first surface is exp(-spacing) ~ 4.5e-5 at 10 mfp, and
    // a survivor carries weight 1 against a combed population whose total weight
    // has decayed by the same factor.  That is a rare high-weight event, and the
    // hypothesis under test is that it is what produces the heavy per-bunch tail
    // (N_eff/n ~ 0.015 at 55-65 mfp, one bunch in 6400 carrying ~10% of a shell).
    bool           pri_diag;         // input: 'primary crossing diagnostic'
    vector<double> pri_cross_sum;    // Σ over bunches of crossings, per surface
    vector<double> pri_cross_wsum;   // Σ over bunches of crossing weight
    vector<double> bank_wsum;        // Σ of weight BANKED at each surface, i.e.
                                     // the scattered channel the primaries are
                                     // competing with.  Their ratio is the
                                     // decisive number and it does not require
                                     // observing the rare deep crossing: it says
                                     // what fraction of the weight arriving at a
                                     // surface is carried by analog primaries.
    vector<long long> pri_cross_b;   // crossings this bunch, per surface
    int            deepest_pri_b;    // deepest surface a primary reached, this bunch
    vector<int>    deep_hist;        // # bunches by deepest surface reached (+1 slot)
    vector<int>    max_deep;         // deepest_pri of the bunch holding max_ratio[k]

    // ---- Scattered-channel comb()/replay weight ledger ---------------------
    // bank_wsum (above) answers "how much weight arrives at a surface".  These
    // answer a different question: does comb()'s stochastic rounding actually
    // conserve that weight once it starts passing through the container ->
    // comb -> replay pipeline, many hundreds to thousands of times per job?
    //
    // comb() is unbiased PER PARTICLE by construction (E[ncopies_j] = w_j/wbar
    // exactly), so comb_after_w should equal comb_before_w in expectation at
    // every call, with a per-call relative fluctuation of order 1/sqrt(n_target)
    // from the stochastic rounding.  Summed over many calls, that fluctuation
    // should average down -- if the CUMULATIVE ratio does not approach 1 as
    // n_comb_calls grows, that is a real bug in comb() or its call sites, not
    // sampling noise, and it is exactly the kind of placement-dependent
    // artifact that could explain why A, B and C disagree with each other
    // while each one's own primary-channel weight ledger balances exactly.
    //
    // Instrumented at BOTH call sites: the periodic Stage-0 downsizing
    // (containers[m]->comb() during source transport) and the cascade replay
    // loop (snapshots[m]->comb() before replayContainer()).  replay_in_w is
    // the post-comb weight actually handed to shower() -- the input side of
    // the one operation (replay) that the primary-channel audit cannot see at
    // all, since replayed particles are ordinary scattered photons by then.
    vector<double>    comb_before_w;   // Σ container weight immediately before comb()
    vector<double>    comb_after_w;    // Σ container weight immediately after
    vector<long long> n_comb_calls;    // how many times comb() ran on this surface
    vector<double>    replay_in_w;     // Σ weight actually fed into shower() via replay


    /*------------------------------------------------------------------------
      scoreFD_all

      Ray-trace from the current photon position along direction u,
      accumulating Lambda, the number of mean free paths traversed.

      At each scoring shell encountered:
        score += w * exp(-Lambda) * (E*muen/rho) * t_shell / V_shell[k]
        accumulated in sum_pri (primary) or sum_tot (total).

      For scattered photons only — stop at the first combing shell hit:
        save a virtual photon at the shell EXIT (iron side) with weight
        w * exp(-Lambda), then return the combing shell index m.

      For primary photons — continue through ALL scoring shells and return -1.

      Key correctness property: virtual photons are placed at the shell EXIT
      so that replayed particles begin in iron and do not re-score the combing
      shell on the first step of their next free path.

      Returns -1 for primaries (always) or when geometry is exited without
      hitting a combing shell.
    ------------------------------------------------------------------------*/
    int scoreFD_all(bool is_primary, bool do_score = true);

    // Replay every particle in src as an independent shower.
    void replayContainer(TmpPhsp *src);

    // Read scoring arrays into per-bunch tables then reset arrays.
    void endBunch();

    static string revision;
};

string EGS_ShieldApplication::revision = "$Revision: 0.3 $";


/*----------------------------------------------------------------------------
  scoreFD_all
----------------------------------------------------------------------------*/
int EGS_ShieldApplication::scoreFD_all(bool is_primary, bool do_score) {
    int np_idx = the_stack->np - 1;
    if (the_stack->E[np_idx] < the_bounds->pcut) return -1;

    EGS_Vector x(the_stack->x[np_idx],
                  the_stack->y[np_idx],
                  the_stack->z[np_idx]);
    EGS_Vector u(the_stack->u[np_idx],
                  the_stack->v[np_idx],
                  the_stack->w[np_idx]);
    int       ireg  = the_stack->ir[np_idx] - 2;
    EGS_Float gle   = the_epcont->gle;
    EGS_Float E     = the_stack->E[np_idx];
    EGS_Float wt    = the_stack->wt[np_idx];
    int       latch = the_stack->latch[np_idx];

    fc_found = false;        // reset the forced-collision record for this trace

    EGS_Float Lambda = 0.0;
    int       imed   = -2;   // sentinel so first medium triggers refresh
    EGS_Float sigma  = 0.0;
    int       newmed = geometry->medium(ireg);

    while (true) {
        // Refresh attenuation coefficient whenever medium changes.
        if (newmed != imed) {
            imed = newmed;
            if (imed >= 0) {
                EGS_Float gmfp = i_gmfp[imed].interpolateFast(gle);
                if (the_xoptions->iraylr)
                    gmfp *= i_cohe[imed].interpolateFast(gle);
                sigma = 1.0 / gmfp;
            } else {
                sigma = 0.0;
            }
        }

        // Advance to next geometry boundary.
        EGS_Float tstep = TSTEP_MAX;
        int inew = geometry->howfar(ireg, x, u, tstep, &newmed);

        Lambda += tstep * sigma;

        if (inew < 0) return -1;   // photon left the geometry

        x    += u * tstep;
        ireg  = inew;
        if (ireg < 0 || ireg >= nreg) return -1;   // safety

        // ---- Entered a scoring shell ----
        if (is_scoring[ireg]) {
            int k = scoring_idx[ireg];

            // Determine the path length through the shell.
            EGS_Float t_shell = TSTEP_MAX;
            int       med_after;
            int ireg_after = geometry->howfar(ireg, x, u, t_shell, &med_after);

            if (t_shell < TSTEP_MAX && do_score) {
                EGS_Float emuen_rho = E_Muen_Rho->interpolateFast(gle);
                EGS_Float unatt = wt * emuen_rho * t_shell / V_shell[k];
                EGS_Float score = unatt * std::exp(-Lambda);
                sum_tot[k] += score;               // total = primary + scattered
                if (is_primary) {
                    sum_pri[k]  += score;
                    // Same score with the attenuation factor removed.  The
                    // ratio of the two sums is the primary transmission along
                    // the FD rays, from which outputResults() recovers the
                    // effective number of mean free paths.  Costs one multiply-add per
                    // primary FD score.
                    sum_pri0[k] += unatt;
                }
            }

            // Traverse the scoring shell (air attenuation negligible, no Lambda update).
            if (t_shell >= TSTEP_MAX || ireg_after < 0) return -1;
            x     += u * t_shell;
            ireg   = ireg_after;
            newmed = med_after;
            // imed refreshed at top of next iteration.
        }

        // ---- Combing check — independent of scoring ----
        // Fires whenever a photon's ray reaches any combing region, whether or
        // not that region is also a scoring shell.  The virtual photon is
        // placed at the combing region's EXIT so the replayed particle starts
        // in the correct medium on the far side.
        if (is_combing[ireg]) {
            EGS_Float t_comb = TSTEP_MAX;
            int       med_after;
            int ireg_after = geometry->howfar(ireg, x, u, t_comb, &med_after);
            // Include attenuation through the combing shell itself.
            // sigma here is for the medium just entered (same iron as surroundings).
            if (t_comb < TSTEP_MAX) Lambda += t_comb * sigma;

            if (is_primary) {
                // Record the first surface for the forced-collision split and
                // KEEP GOING: a primary's FD trace must run the full ray, since
                // it is the single deterministic pass that books K_pri at every
                // shell out to 100 mfp.  Stopping here would truncate K_pri
                // exactly the way the missing FD key truncated it on 2026-08-12.
                if (!fc_found && t_comb < TSTEP_MAX && ireg_after >= 0) {
                    fc_found  = true;
                    fc_Lambda = Lambda;
                    fc_x      = x + u * t_comb;   // combing region EXIT
                    fc_ireg   = ireg_after;
                    fc_m      = combing_idx[ireg];
                    // Trace-only call: the record is all it was for.
                    if (!do_score) return fc_m;
                }
                // Traverse the combing region and carry on down the ray.
                if (t_comb >= TSTEP_MAX || ireg_after < 0) return -1;
                x     += u * t_comb;
                ireg   = ireg_after;
                newmed = med_after;
                continue;
            }

            if (t_comb < TSTEP_MAX && ireg_after >= 0) {
                EGS_Particle vp;
                vp.q     = 0;
                vp.latch = latch;
                vp.E     = E;
                vp.wt    = wt * std::exp(-Lambda);
                vp.x     = x + u * t_comb;   // combing region EXIT
                vp.u     = u;
                vp.ir    = ireg_after;
                containers[combing_idx[ireg]]->save(vp);
                if (pri_diag) bank_wsum[combing_idx[ireg]] += vp.wt;
            }
            return combing_idx[ireg];
        }
    }
}


/*----------------------------------------------------------------------------
  initScoring
----------------------------------------------------------------------------*/
int EGS_ShieldApplication::initScoring() {
    EGS_Input *options = input->takeInputItem("scoring options");
    if (!options) {
        egsWarning("egs_shield::initScoring: no 'scoring options' block\n");
        return 1;
    }

    // ---- E*mu_en/rho table ----
    string emuen_file;
    if (options->getInput("emuen file", emuen_file)) {
        egsWarning("egs_shield::initScoring: 'emuen file' not specified\n");
        delete options; return 1;
    }
    emuen_file = egsExpandPath(emuen_file);
    {
        std::ifstream mf(emuen_file.c_str());
        if (!mf) {
            egsWarning("egs_shield::initScoring: cannot open '%s'\n",
                       emuen_file.c_str());
            delete options; return 1;
        }
        int ndat; mf >> ndat;
        if (ndat < 2 || mf.fail())
            egsFatal("egs_shield::initScoring: bad emuen file header\n");
        EGS_Float *xm = new EGS_Float[ndat], *fm = new EGS_Float[ndat];
        for (int j = 0; j < ndat; j++) mf >> xm[j] >> fm[j];
        if (mf.fail())
            egsFatal("egs_shield::initScoring: read error in emuen file\n");
        E_Muen_Rho = new EGS_Interpolator(ndat,
                         std::log(xm[0]), std::log(xm[ndat-1]), fm);
        delete [] xm; delete [] fm;
    }

    // ---- Geometry ----
    string geom_name;
    if (options->getInput("geometry name", geom_name)) {
        egsWarning("egs_shield::initScoring: 'geometry name' not specified\n");
        delete options; return 1;
    }
    geometry = EGS_BaseGeometry::getGeometry(geom_name);
    if (!geometry) {
        egsWarning("egs_shield::initScoring: geometry '%s' not found\n",
                   geom_name.c_str());
        delete options; return 1;
    }
    nreg = geometry->regions();

    // ---- Scoring regions ----
    vector<int> s_regs;
    options->getInput("scoring regions", s_regs);
    n_scoring = (int)s_regs.size();
    if (n_scoring < 1) {
        egsWarning("egs_shield::initScoring: 'scoring regions' not given\n");
        delete options; return 1;
    }

    // ---- Scoring volumes ----
    vector<EGS_Float> s_vols;
    options->getInput("scoring volumes", s_vols);
    if ((int)s_vols.size() != n_scoring) {
        egsWarning("egs_shield::initScoring: "
                   "'scoring volumes' must have %d entries\n", n_scoring);
        delete options; return 1;
    }

    scoring_reg = new int[n_scoring];
    V_shell     = new EGS_Float[n_scoring];
    scoring_idx = new int[nreg];
    is_scoring  = new bool[nreg];
    std::fill(scoring_idx, scoring_idx + nreg, -1);
    std::fill(is_scoring,  is_scoring  + nreg, false);

    for (int k = 0; k < n_scoring; k++) {
        int r = s_regs[k];
        if (r < 0 || r >= nreg)
            egsFatal("egs_shield::initScoring: scoring region %d out of "
                     "range [0,%d)\n", r, nreg);
        scoring_reg[k] = r;
        scoring_idx[r] = k;
        is_scoring[r]  = true;
        V_shell[k]     = s_vols[k];
    }

    // ---- Combing regions (thin iron shells; need NOT be scoring regions) ----
    vector<int> c_regs;
    options->getInput("combing regions", c_regs);
    n_combing = (int)c_regs.size();
    if (n_combing < 1) {
        egsWarning("egs_shield::initScoring: 'combing regions' not given\n");
        delete options; return 1;
    }

    combing_reg = new int[n_combing];
    combing_idx = new int[nreg];
    is_combing  = new bool[nreg];
    std::fill(combing_idx, combing_idx + nreg, -1);
    std::fill(is_combing,  is_combing  + nreg, false);

    for (int m = 0; m < n_combing; m++) {
        int r = c_regs[m];
        if (r < 0 || r >= nreg)
            egsFatal("egs_shield::initScoring: combing region %d out of range "
                     "(nreg=%d)\n", r, nreg);
        combing_reg[m] = r;
        combing_idx[r] = m;
        is_combing[r]  = true;
    }

    // ---- Kerma accumulators ----
    sum_tot.assign(n_scoring, 0.0);
    sum_pri.assign(n_scoring, 0.0);
    sum_pri0.assign(n_scoring, 0.0);

    // ---- Containers and snapshot buffers ----
    containers = new TmpPhsp*[n_combing];
    snapshots  = new TmpPhsp*[n_combing];
    for (int m = 0; m < n_combing; m++) {
        containers[m] = new TmpPhsp();
        snapshots[m]  = new TmpPhsp();
    }

    // ---- Bunch parameters ----
    int nb = 100;
    options->getInput("n bunches", nb);
    n_bunches = nb;

    int npb = 1000000;
    options->getInput("photons per bunch", npb);
    n_per_bunch = (long long)npb;

    // Opt-in replay-cascade trace; off by default.  Useful when adding combing
    // shells or changing their spacing, to confirm every container is actually
    // reached before the loop's weight threshold stops it.  Prints ~10 lines
    // per job, so leave it off for large parallel runs.
    vector<string> choice;
    choice.push_back("no");
    choice.push_back("yes");
    cascade_diag = options->getInput("cascade diagnostic", choice, 0) ? true : false;

    // Per-shell bunch-contribution statistics: how concentrated a result is in
    // a few bunches.  Cheap, off by default; turn on when a deep shell's sigma
    // looks untrustworthy.
    bunch_stats  = options->getInput("bunch statistics", choice, 0) ? true : false;

    // Primary-crossing diagnostic.  Per-job output only: it is deliberately NOT
    // written to the .egsdat, because adding a field there breaks combining
    // against every file produced by an older binary (as the sum_Kp0 field did
    // on 2026-08-10).  A diagnostic run needs tens of bunches, not thousands,
    // so per-job tables are sufficient.
    pri_diag     = options->getInput("primary crossing diagnostic", choice, 0)
                   ? true : false;

    // Forced collision for primaries between combing surfaces.  Off by default:
    // it changes the transport, so an existing input must keep reproducing its
    // existing answer unless the change is explicitly asked for.  With it off
    // every code path below is the one that produced the validated results.
    forced_collision = options->getInput("forced collision", choice, 0)
                       ? true : false;
    if (forced_collision && n_combing < 1) {
        egsWarning("\n*** 'forced collision = yes' but no combing regions are"
                   " defined:\n*** there is no surface to force against."
                   "  Ignoring.\n\n");
        forced_collision = false;
    }

    // Replay-cascade weight cutoff.  The cascade stops once the surviving weight
    // falls below this fraction of its initial value.
    //
    // CRITICAL: this must be smaller than exp(-eta_max) for the deepest shell,
    // or the cascade is truncated before it reaches the outer combing shells and
    // the deep results are silently wrong.  At 100 mfp exp(-100) = 4e-44, so a
    // cutoff of 1e-10 stops the cascade at roughly 44 mfp -- with containers
    // still full, so their pending contributions are discarded as well and the
    // damage reaches inward.  Symptom: BUF -> exactly 1.000 (K_tot == K_pri) at
    // the deep shells, with the shells just inside them biased low.
    //
    // The default 1e-10 is kept for backward compatibility with the 9-shell
    // inputs, where the cascade completed in ~10 iterations and happened to
    // finish just as the threshold fired.  Any geometry with more combing
    // shells, or reaching deeper, must set this explicitly.
    cascade_cutoff = 1e-10;
    options->getInput("cascade weight cutoff", cascade_cutoff);
    if (cascade_cutoff <= 0 || cascade_cutoff >= 1) {
        egsWarning("\n*** 'cascade weight cutoff' = %g is outside (0,1);"
                   " using 1e-10\n", cascade_cutoff);
        cascade_cutoff = 1e-10;
    }

    // Population target for the combing operation.  n_target <= 0 selects the
    // Divide-et-Impera rule proper: w_bar = sum(w_i)/N, which preserves the
    // arriving population.  A positive value caps the population, at the cost
    // of culling arrivals above it.
    //
    // Default 0 means "derive from photons per bunch" (comb_interval*10), the
    // historical behaviour.  Lower it when running many combing shells: peak
    // memory is roughly 2 * n_combing * (container size) * 80 bytes, and the
    // container size tracks this target.
    comb_target_in = 0;
    options->getInput("comb target", comb_target_in);

    sum_ratio  = new double[n_scoring]();
    sum_ratio2 = new double[n_scoring]();
    sum_Kt     = new double[n_scoring]();
    sum_Kp     = new double[n_scoring]();
    sum_Kp0    = new double[n_scoring]();
    max_ratio  = new double[n_scoring]();
    n_valid_b  = new int[n_scoring]();
    n_completed = 0;

    job_nr2.assign(n_scoring, 0.0);
    job_used.assign(n_scoring, 0);
    topk.assign((size_t)n_scoring * N_TOPK, 0.0);
    bunch_r.assign(n_scoring, vector<double>());
    fc_in_w.assign(n_combing, 0.0);
    fc_bank_w.assign(n_combing, 0.0);
    fc_coll_w.assign(n_combing, 0.0);
    fc_n_split.assign(n_combing, 0);
    max_deep.assign(n_scoring, -1);
    pri_cross_sum.assign(n_combing, 0.0);
    pri_cross_wsum.assign(n_combing, 0.0);
    bank_wsum.assign(n_combing, 0.0);
    pri_cross_b.assign(n_combing, 0);
    deep_hist.assign(n_combing + 1, 0);   // slot 0 = no surface reached
    comb_before_w.assign(n_combing, 0.0);
    comb_after_w.assign(n_combing, 0.0);
    n_comb_calls.assign(n_combing, 0);
    replay_in_w.assign(n_combing, 0.0);

    // Enable ausgab calls needed for latch bookkeeping.
    // By default the framework only enables BeforeTransport..AfterTransport.
    // AfterCompton/Photo/Pair/Rayleigh are above AfterTransport in the enum
    // and must be turned on explicitly so that scattered photons get latch > 0.
    setAusgabCall(AfterCompton,  true);
    setAusgabCall(AfterPhoto,    true);
    setAusgabCall(AfterPair,     true);
    setAusgabCall(AfterRayleigh, true);

    // ---- Variance-reduction summary --------------------------------------
    // egs_shield now layers three independent VR mechanisms, and a run's
    // .egslog previously gave no single place to confirm which were active.
    // That is exactly the gap that let egs_kerma silently fall back to
    // track-length scoring for every input this generator produced until
    // 2026-08-13 (12.6.17/12.6.20 in the research log) -- a missing or
    // misread input key with no confirming banner cost a full production
    // run before anyone noticed.  Print the equivalent confirmation here so
    // 'was X actually on' is answered by a `grep` of the log, not by an
    // inference from the input file.
    egsInformation(
        "\n===========================================================\n"
        " Variance reduction techniques\n"
        "===========================================================\n\n");
    egsInformation(
        " Forced detection (FD):        ON  (always -- FD is egs_shield's\n"
        "                                     core estimator; fires at every\n"
        "                                     free path via scoreFD_all)\n\n");
    egsInformation(
        " Population control (combing): %s\n", n_combing > 0 ? "ON" : "OFF");
    if (n_combing > 0) {
        egsInformation("   combing surfaces    = %d\n", n_combing);
        egsInformation("   combing regions     =");
        for (int m = 0; m < n_combing; m++) egsInformation(" %d", combing_reg[m]);
        egsInformation("\n");
        if (comb_target_in > 0)
            egsInformation("   comb target         = %d (fixed population cap)\n",
                           comb_target_in);
        else
            egsInformation("   comb target         = preserve arriving population"
                           " (D&I rule, wbar = sum(w)/N)\n");
        egsInformation("   cascade weight cutoff = %g\n", cascade_cutoff);
    }
    egsInformation("\n Forced collision (primaries): %s\n",
                   forced_collision ? "ON" : "OFF");
    if (forced_collision) {
        egsInformation("   Between each pair of combing surfaces, a primary's\n"
                       "   weight is split deterministically into an\n"
                       "   uncollided branch (banked at the far surface) and a\n"
                       "   collided branch (forced to interact inside the\n"
                       "   segment).  K_pri is unaffected -- see 12.6.20/12.6.22\n"
                       "   in the research log for the mechanism.  CONFIRM the\n"
                       "   RNG in the block above is xoshiro256++, not a silent\n"
                       "   fallback -- that exact failure mode (12.6.22) produced\n"
                       "   a spurious combing-surface-placement dependence that\n"
                       "   looked like a real bias in this estimator until traced\n"
                       "   to the RNG on 2026-08-17.\n");
    }
    egsInformation("\n===========================================================\n\n");

    egsInformation("egs_shield: %d scoring shells, %d combing shells, "
                   "%d bunches of %lld photons\n",
                   n_scoring, n_combing, n_bunches, n_per_bunch);

    delete options;
    return 0;
}


/*----------------------------------------------------------------------------
  replayContainer — transport each particle in src as a separate shower.
----------------------------------------------------------------------------*/
void EGS_ShieldApplication::replayContainer(TmpPhsp *src) {
    for (int j = 0; j < src->size(); j++) {
        EGS_Particle &par = (*src)[j];
        p.E     = par.E;
        p.wt    = par.wt;
        p.q     = par.q;
        p.latch = par.latch;
        p.x     = par.x;
        p.u     = par.u;
        p.ir    = par.ir;
        ++current_case;
        startNewShower();
        shower();
        finishShower();
    }
}


/*----------------------------------------------------------------------------
  initRunControl — force EGS_UniformRunControl (URCO) for parallel runs.

  egs_shield distributes bunches among jobs itself and needs no JCF lock
  file.  If the framework created a JCF (the default when "rco type" is
  absent from the run control input block), swap it out for URCO.

  Preferred usage: add "rco type = uniform" (plus timeout settings) to the
  run control input block so URCO is created with correct configuration.
  This override is a safety net for egsinp files that omit that key.
----------------------------------------------------------------------------*/
int EGS_ShieldApplication::initRunControl() {
    int err = EGS_AdvancedApplication::initRunControl();
    if (err) return err;
    if (dynamic_cast<EGS_JCFControl *>(run)) {
        EGS_I64 nc = run->getNcase();
        delete run;
        run = new EGS_UniformRunControl(this);
        run->setNcase(nc);
        egsInformation("egs_shield: replaced JCF with URCO "
                       "(add 'rco type = uniform' to run control block "
                       "to configure watcher timeout)\n");
    }
    return 0;
}


/*----------------------------------------------------------------------------
  runSimulation — staged + iterative backscatter-aware loop.

  APP_MAIN calls finishSimulation() which calls outputResults() after we
  return, so we do NOT call outputResults() here.
----------------------------------------------------------------------------*/
int EGS_ShieldApplication::runSimulation() {
    // Let the RCO do its startup work (deletes stale .egsdat, handles
    // resume/combine/analyze modes).  Without this, stale .egsdat files
    // from a previous parallel run are counted as done by URCO's watcher.
    int start_err = run->startSimulation();
    if (start_err != 0) {
        if (start_err < 0) {
            egsWarning("egs_shield: run control start failed (%d)\n", start_err);
        }
        return 0;   // combine/analyze mode: nothing to simulate
    }

    // In a parallel run (getNparallel() > 0, getIparallel() > 0), divide
    // n_bunches evenly across jobs so total work equals a serial run.
    // Each job uses a different RNG sequence (set by the parallel framework),
    // so bunches from different jobs are independent and can be combined.
    int my_n_bunches = n_bunches;
    if (getNparallel() > 0 && getIparallel() > 0) {
        int np = getNparallel();
        int ip = getIparallel() - getFirstParallel();  // 0-based
        my_n_bunches = n_bunches / np + (ip < n_bunches % np ? 1 : 0);
        egsInformation("egs_shield: parallel job %d/%d — running %d of %d bunches\n",
                       getIparallel(), np, my_n_bunches, n_bunches);
    }

    EGS_Timer timer;
    for (current_bunch = 0; current_bunch < my_n_bunches; current_bunch++) {
        timer.start();
        egsInformation("egs_shield: bunch %d / %d  (%lld photons)",
                       current_bunch + 1, my_n_bunches, n_per_bunch);

        for (int m = 0; m < n_combing; m++) containers[m]->clean();

        if (pri_diag) {
            std::fill(pri_cross_b.begin(), pri_cross_b.end(), 0);
            deepest_pri_b = -1;
        }

        // ---- Stage 0: source photons ----
        // Periodic combing every comb_interval photons keeps containers from
        // growing without bound.  We cap the population at comb_interval×10 so
        // that the container never grows much larger than that ceiling.
        const long long comb_interval = std::max(10000LL, n_per_bunch / 100);
        const int       comb_target   = comb_target_in > 0
                                      ? comb_target_in
                                      : (int)(comb_interval * 10);

        for (long long ih = 0; ih < n_per_bunch; ih++) {
            EGS_Vector x0, u0;
            current_case = source->getNextParticle(rndm,
                               p.q, p.latch, p.E, p.wt, x0, u0);
            p.x = x0; p.u = u0;
            p.latch = 0;   // mark as primary

            // Place particle inside geometry (handles isotropic sources).
            int ireg = geometry->isWhere(p.x);
            if (ireg < 0) {
                EGS_Float t = 1e30;
                ireg = geometry->howfar(-1, p.x, p.u, t);
                if (ireg >= 0) p.x += p.u * t;
            }
            if (ireg < 0) continue;
            p.ir = ireg;

            startNewShower();
            shower();
            finishShower();

            if ((ih + 1) % comb_interval == 0) {
                for (int m = 0; m < n_combing; m++) {
                    if (pri_diag) {
                        comb_before_w[m] += containers[m]->weight();
                        n_comb_calls[m]++;
                    }
                    containers[m]->comb(rndm, comb_target);
                    if (pri_diag) comb_after_w[m] += containers[m]->weight();
                }
            }
        }

        // ---- Iterative backscatter loop ----
        //
        // Each iteration:
        //   1. Check whether any containers hold photons; stop if all empty.
        //   2. Snapshot <- containers (swap), clear containers.
        //   3. Comb each snapshot.
        //   4. Replay all snapshot particles; this may refill containers
        //      (forward scatter to deeper combing shells, backward scatter
        //      to shallower combing shells — both handled identically).
        //
        // Convergence is geometric: each pass attenuates by exp(-Λ_min)
        // where Λ_min is the number of mean free paths to the nearest combing shell.
        // Track total weight to exit early once the remaining weight is
        // negligible relative to its initial value (< 1e-10).
        //
        double w_initial = 0;
        for (int m = 0; m < n_combing; m++)
            for (int j = 0; j < containers[m]->size(); j++)
                w_initial += (*containers[m])[j].wt;

        // Cascade diagnostic (first bunch only, so cost and log noise are
        // negligible).  Prints, per replay iteration, the surviving weight
        // relative to w_initial and the occupancy of every combing container.
        // Its purpose is to show how deep the cascade actually propagates:
        // a scoring shell beyond combing shell m can only be scored once
        // container m has been replayed, so a container that never fills is
        // a scoring hole, not merely a slow tail.
        const bool diag = cascade_diag && (current_bunch == 0);
        if (diag) {
            egsInformation("\n  [cascade] bunch 0, %d combing shells, "
                           "w_initial = %.6g\n", n_combing, w_initial);
            egsInformation("  [cascade] %4s %12s %12s   %s\n",
                           "iter", "w_cur", "w_cur/w_ini", "container sizes");
        }

        static const int MAX_ITER = 500;
        /* Has the cascade frontier ever reached each container?  Used only for
         * the truncation warning below. */
        vector<bool> ever_filled(n_combing, false);
        for (int iter = 0; iter < MAX_ITER; iter++) {
            for (int m = 0; m < n_combing; m++)
                if (containers[m]->size() > 0) ever_filled[m] = true;
            double w_cur = 0;
            for (int m = 0; m < n_combing; m++)
                for (int j = 0; j < containers[m]->size(); j++)
                    w_cur += (*containers[m])[j].wt;

            if (diag) {
                string occ;
                char b[32];
                for (int m = 0; m < n_combing; m++) {
                    sprintf(b, "%d", containers[m]->size());
                    occ += b;
                    if (m + 1 < n_combing) occ += " ";
                }
                egsInformation("  [cascade] %4d %12.4g %12.3g   %s\n",
                               iter, w_cur,
                               w_initial > 0 ? w_cur / w_initial : 0.0,
                               occ.c_str());
            }

            if (w_cur <= 0 || w_cur < cascade_cutoff * w_initial) {
                if (diag)
                    egsInformation("  [cascade] STOP at iter %d: w_cur/w_ini = %.3g"
                                   " < %g cutoff\n\n", iter,
                                   w_initial > 0 ? w_cur / w_initial : 0.0,
                                   cascade_cutoff);
                /* Truncation test: did the cascade frontier ever REACH every
                 * combing container?  A container that was never populated
                 * means no particle ever got that deep, so every scoring shell
                 * beyond it received primaries only and reads BUF = 1 exactly.
                 *
                 * Do NOT test for containers still holding particles at the
                 * stop: they always do.  Backscatter refills the shallower
                 * containers indefinitely, so the cascade settles into a
                 * quasi-steady population while its weight decays -- occupancy
                 * never falls to zero even when fully converged. */
                if (current_bunch == 0) {
                    int unreached = 0, first = -1;
                    for (int m = 0; m < n_combing; m++)
                        if (!ever_filled[m]) {
                            unreached++;
                            if (first < 0) first = m;
                        }
                    if (unreached > 0) {
                        /* Self-calibrating recommendation.  The frontier advances
                         * about one container per iteration, and the weight falls
                         * by a roughly constant factor per iteration -- measure
                         * that factor from this run rather than assuming one, as
                         * it depends strongly on the combing-shell spacing. */
                        double per_iter = (iter > 0 && w_initial > 0 && w_cur > 0)
                                        ? std::pow(w_cur / w_initial, 1.0 / iter)
                                        : 0.1;
                        double need = std::pow(per_iter, (double)(n_combing + 6));
                        egsWarning("\n*** egs_shield: replay cascade TRUNCATED. %d of %d combing"
                                   " container(s) were never\n"
                                   "*** reached -- the frontier stopped at container %d after %d"
                                   " iteration(s).\n"
                                   "*** Every scoring shell beyond that combing shell receives"
                                   " primaries only and\n"
                                   "*** will read BUF = 1 exactly (K_tot == K_pri); shells just"
                                   " inside it read low.\n"
                                   "***\n"
                                   "*** The frontier advances ~1 container per iteration while the"
                                   " weight falls %.3g per\n"
                                   "*** iteration in this run, so reaching all %d containers needs"
                                   " 'cascade weight\n"
                                   "*** cutoff' below about %.0e.  It is currently %g.\n\n",
                                   unreached, n_combing, first - 1, iter,
                                   per_iter, n_combing, need, cascade_cutoff);
                    }
                }
                break;
            }

            for (int m = 0; m < n_combing; m++) {
                std::swap(containers[m], snapshots[m]);
                containers[m]->clean();
                // Cap at current size: don't split when np < comb_target.
                // Splitting would create a fixed-cost 100K-shower replay
                // regardless of n_per_bunch, breaking time scaling.
                if (pri_diag) {
                    comb_before_w[m] += snapshots[m]->weight();
                    n_comb_calls[m]++;
                }
                snapshots[m]->comb(rndm, std::min(snapshots[m]->size(), comb_target));
                if (pri_diag) comb_after_w[m] += snapshots[m]->weight();
            }
            for (int m = 0; m < n_combing; m++)
                if (snapshots[m]->size() > 0) {
                    if (pri_diag) replay_in_w[m] += snapshots[m]->weight();
                    replayContainer(snapshots[m]);
                }
        }

        double bunch_time = timer.time();
        total_cpu_time_ += bunch_time;
        endBunch();
        egsInformation("  %.1f s\n", bunch_time);

        // We replace the standard batch loop, so EGS_RunControl::finishBatch()
        // is never reached.  Call it here, once per bunch, exactly where the
        // batch loop would: it sets the run control's cpu_time and then calls
        // outputData() for us.  Doing this for EVERY parallel job -- the
        // watcher included -- is what gets every job's results into the
        // combine.
        //
        // Per bunch rather than once at the end, so a job that dies partway
        // through still contributes its completed bunches.  The .egsdat holds
        // cumulative sums and is simply overwritten, so the newest file always
        // reflects every bunch finished so far.  This matters on a cluster:
        // OOM kills land mid-cascade, at peak memory, and previously took every
        // bunch of that job with them.
        //
        // It must happen here and not in finishSimulation(): for the watcher
        // the whole combine runs inside EGS_AdvancedApplication::
        // finishSimulation(), and the partial combine inside URCO's poll loop
        // calls resetCounter(), which zeroes our accumulators.  Writing
        // afterwards would store zeros into a file nothing reads any more.
        //
        // run_dir is still set at this point, so the file lands in this job's
        // egsrun_* directory: invisible to howManyJobsDone() (which scans the
        // application directory only), so the watcher's poll loop still waits
        // for npar-1 *other* jobs.  finishRun() then moves it up to the
        // application directory before combineResults() runs.
        //
        // finishBatch() also prints a batch summary line and applies the
        // "statistical accuracy sought" early-termination test.  We do not
        // override getCurrentResult(), so that line's result/uncertainty
        // columns print as 0 and 100.00; the return value is only meaningful
        // for loop control, which we do not have, so it is discarded.
        if (getNparallel() > 0 && getIparallel() > 0) {
            run->finishBatch();
        }
    }
    return 0;
}


/*----------------------------------------------------------------------------
  endBunch
----------------------------------------------------------------------------*/
void EGS_ShieldApplication::endBunch() {
    for (int k = 0; k < n_scoring; k++) {
        double Kt  = sum_tot[k]  / n_per_bunch;
        double Kp  = sum_pri[k]  / n_per_bunch;
        double Kp0 = sum_pri0[k] / n_per_bunch;
        sum_tot[k]  = 0.0;
        sum_pri[k]  = 0.0;
        sum_pri0[k] = 0.0;
        if (Kp > 0.0) {
            double ratio   = Kt / Kp;
            sum_ratio[k]  += ratio;
            sum_ratio2[k] += ratio * ratio;
            sum_Kt[k]     += Kt;
            sum_Kp[k]     += Kp;
            sum_Kp0[k]    += Kp0;
            insertTopK(k, ratio);
            bunch_r[k].push_back(ratio);
            if (ratio > max_ratio[k]) {
                max_ratio[k] = ratio;
                // Which bunch owns the maximum, characterised by how deep its
                // primaries got.  If the outlier bunches are systematically the
                // ones where a primary reached a deep surface, the analog
                // treatment of primaries is the tail's source.
                if (pri_diag) max_deep[k] = deepest_pri_b;
            }
            n_valid_b[k]++;
        }
    }
    if (pri_diag) {
        for (int m = 0; m < n_combing; m++)
            pri_cross_sum[m] += (double)pri_cross_b[m];
        deep_hist[deepest_pri_b + 1]++;
    }
    n_completed++;
}


/*----------------------------------------------------------------------------
  outputResults

  Effective number of mean free paths (eta/mfp column)
  -----------------------------------------------------
  The FD estimator already carries the exact quantity we want.  Every primary
  score is w·exp(-Λ)·(E μen/ρ)·t/V, where Λ is the number of mean free paths
  accumulated along the ray from the source to the shell (path length divided
  by the local mean free path length, summed across any media crossed).
  Accumulating the same score without exp(-Λ) (sum_Kp0) makes the ratio
  ΣK_pri/ΣK_pri0 the primary transmission, so

      η_eff = -ln( ΣK_pri / ΣK_pri0 )

  For a monoenergetic point source and concentric shells every primary ray to
  a given shell is identical, so η_eff is exactly that shell's number of mean
  free paths, to machine precision and with no user input.

  For a polyenergetic source the per-photon η varies with energy and η_eff
  becomes -ln⟨exp(-η)⟩, the transmission-weighted mean — which is the right
  characteristic depth here, because it is precisely the attenuation that
  forms the denominator of the buildup factor.  It is not the arithmetic mean
  number of mean free paths, and for a broad spectrum it will sit below it
  (low-energy components are attenuated away and stop contributing).  The same expression
  also handles extended or off-axis sources, where different primaries reach a
  shell along different chords.  Deriving it beats an input-supplied mfp,
  which would additionally assume a spherical geometry and a source at the
  centre.
----------------------------------------------------------------------------*/
void EGS_ShieldApplication::outputResults() {
    static const double MeVtoGy = 1.6021773e-10;

    // Warn when the reported result is not the requested one.  Fires on serial
    // runs and on the final parallel combine (i_parallel reset to 0 before
    // combineResults()), but not on an individual job's own log, where
    // n_completed is only that job's share of the bunches.
    bool is_single_job = (getNparallel() > 0 && getIparallel() > 0);
    if (!is_single_job && n_completed != n_bunches) {
        egsWarning("\n"
                   "  ****************************************************************\n"
                   "  *** INCOMPLETE: %d of %d bunches present in this result.\n"
                   "  *** %d bunch(es) are missing -- %.1f%% of the requested work\n"
                   "  *** is not included below.\n"
                   "  *** In a parallel run this normally means the watcher job\n"
                   "  *** stopped waiting before the other jobs finished.  Raise\n"
                   "  *** 'number of intervals' and/or 'interval wait time' in the\n"
                   "  *** run control block and re-combine.\n"
                   "  ****************************************************************\n",
                   n_completed, n_bunches, n_bunches - n_completed,
                   100.0 * (n_bunches - n_completed) / std::max(n_bunches, 1));
    }

    egsInformation("\n");
    egsInformation("  %-6s  %-9s  %-14s  %-14s  %-12s  %-10s  %-12s\n",
                   "Region", "eta/mfp", "K_tot/[Gy]", "K_pri/[Gy]", "BUF",
                   "sigma/%", "FOM/[s^-1]");
    egsInformation("  %s\n", string(85, '-').c_str());

    for (int k = 0; k < n_scoring; k++) {
        int nv = n_valid_b[k];

        // eta is deterministic, so report it even when the statistics are not
        // yet usable.
        double eta = (sum_Kp[k] > 0.0 && sum_Kp0[k] > 0.0)
                     ? -std::log(sum_Kp[k] / sum_Kp0[k]) : 0.0;

        if (nv < 2) {
            egsInformation("  %-6d  %-9.4f  (%d bunch(es) in this job -- "
                           "statistics come from the combined run)\n",
                           scoring_reg[k], eta, nv);
            continue;
        }

        double mean    = sum_ratio[k] / nv;
        double var     = (sum_ratio2[k] - sum_ratio[k] * sum_ratio[k] / nv) / (nv - 1);
        double sem     = std::sqrt(std::max(var, 0.0) / nv);
        double Kt_mean = sum_Kt[k] / nv;
        double Kp_mean = sum_Kp[k] / nv;

        double sigma_frac = sem / std::max(mean, 1e-30);
        double fom = (total_cpu_time_ > 0 && sigma_frac > 0)
                     ? 1.0 / (sigma_frac * sigma_frac * total_cpu_time_)
                     : 0.0;

        egsInformation("  %-6d  %-9.4f  %-14.6g  %-14.6g  %-12.6g  %-10.4f  %-12.4g\n",
                       scoring_reg[k], eta,
                       Kt_mean * MeVtoGy, Kp_mean * MeVtoGy, mean,
                       100.0 * sigma_frac, fom);
    }

    // ---- Bunch-contribution diagnostics ------------------------------------
    //
    // sigma above is the standard error of a mean over bunches.  It is only
    // trustworthy if many bunches genuinely contribute.  With a heavy-tailed
    // per-bunch distribution -- which is what deep shells produce, since a
    // shell's score can hinge on one rare high-weight history -- a handful of
    // bunches carry the result, and both the mean and its error are poorly
    // determined however many bunches were run.
    //
    // Kish's effective sample size measures this:
    //
    //     N_eff = (sum x)^2 / sum x^2 = n / (1 + CV^2)
    //
    // with CV the coefficient of variation of the per-bunch values.  N_eff = n
    // when all bunches contribute equally; N_eff << n when few do.  Both sums
    // are already accumulated, so this costs nothing.
    //
    // max_ratio gives the complementary picture the sums cannot: the single
    // largest bunch, and how far above the mean it sits.
    const double neff_warn = 0.05;   // flag shells below 5% effective usage

    int n_flagged = 0;
    for (int k = 0; k < n_scoring; k++) {
        if (n_valid_b[k] < 2 || sum_ratio2[k] <= 0.0) continue;
        double neff = sum_ratio[k] * sum_ratio[k] / sum_ratio2[k];
        if (neff < neff_warn * n_valid_b[k]) n_flagged++;
    }
    if (n_flagged > 0) {
        egsWarning("\n  *** %d shell(s) have an effective bunch count below %.0f%% of\n"
                   "  *** the bunches run: their sigma is unreliable and the mean may\n"
                   "  *** be biased by a few high-weight bunches.  Set\n"
                   "  *** 'bunch statistics = yes' for the per-shell breakdown.\n",
                   n_flagged, 100.0 * neff_warn);
    }

    if (bunch_stats) {
        egsInformation("\n  Per-shell bunch-contribution statistics\n");
        egsInformation("  %-6s  %-9s  %-8s  %-9s  %-8s  %-12s  %-8s\n",
                       "Region", "eta/mfp", "CV", "N_eff", "N_eff/n",
                       "max bunch", "max/mean");
        egsInformation("  %s\n", string(72, '-').c_str());
        for (int k = 0; k < n_scoring; k++) {
            int nv = n_valid_b[k];
            if (nv < 2 || sum_ratio2[k] <= 0.0) continue;
            double mean = sum_ratio[k] / nv;
            double eta  = (sum_Kp[k] > 0.0 && sum_Kp0[k] > 0.0)
                          ? -std::log(sum_Kp[k] / sum_Kp0[k]) : 0.0;
            double neff = sum_ratio[k] * sum_ratio[k] / sum_ratio2[k];
            double var  = (sum_ratio2[k] - sum_ratio[k] * sum_ratio[k] / nv) / (nv - 1);
            double cv   = std::sqrt(std::max(var, 0.0)) / std::max(mean, 1e-30);
            egsInformation("  %-6d  %-9.4f  %-8.3f  %-9.1f  %-8.4f  %-12.6g  %-8.1f\n",
                           scoring_reg[k], eta, cv, neff, neff / nv,
                           max_ratio[k], max_ratio[k] / std::max(mean, 1e-30));
        }
        egsInformation("\n  N_eff = (sum r)^2 / sum r^2, r = per-bunch K_tot/K_pri (Kish).\n"
                       "  N_eff ~ n means every bunch contributes; N_eff << n means a few\n"
                       "  dominate and sigma understates the true uncertainty.\n");
    }

    // ---- Leave-k-out sensitivity -------------------------------------------
    //
    // The question sigma cannot answer: how much of the ANSWER rests on a few
    // bunches?  Dropping the largest and re-forming the mean answers it
    // directly.  A shell where removing 1 bunch in 6400 moves the mean by 4%
    // is not a shell whose 4% sigma means what it looks like it means.
    //
    // This is the diagnostic the between-job comparison cannot provide.  With a
    // single dominant bunch of value X among n, the between-bunch and
    // between-job variances are algebraically forced to agree -- both reduce to
    // X^2/n at any batch size -- so their ratio is blind to exactly this.
    // Absent from .egsdat written before 2026-08-13.  Without the guard the
    // table silently reports n/(n-j)-1 -- subtracting nothing while shrinking
    // the denominator -- which looks like a real and alarming sensitivity.
    bool have_topk = false;
    for (int k = 0; k < n_scoring; k++)
        if (topk[(size_t)k * N_TOPK] > 0.0) { have_topk = true; break; }
    if (bunch_stats && !have_topk) {
        egsInformation("\n  Leave-k-out sensitivity unavailable: no top-k data in the\n"
                       "  .egsdat file(s), which means they were written before this\n"
                       "  diagnostic existed.  Re-run to get it.\n");
    }
    if (bunch_stats && have_topk) {
        egsInformation("\n  Leave-k-out sensitivity of the mean"
                       "   [top %d bunches per shell]\n", N_TOPK);
        egsInformation("  %-6s  %-9s  %-12s", "Region", "eta/mfp", "BUF");
        for (int j = 1; j <= N_TOPK; j++) egsInformation("  drop%-6d", j);
        egsInformation("\n  %s\n", string(50 + 10 * N_TOPK, '-').c_str());
        for (int k = 0; k < n_scoring; k++) {
            int nv = n_valid_b[k];
            if (nv <= N_TOPK + 1) continue;
            double mean = sum_ratio[k] / nv;
            if (mean <= 0.0) continue;
            egsInformation("  %-6d  %-9.4f  %-12.6g", scoring_reg[k],
                           (sum_Kp[k] > 0.0 && sum_Kp0[k] > 0.0)
                           ? -std::log(sum_Kp[k] / sum_Kp0[k]) : 0.0, mean);
            double s = sum_ratio[k];
            for (int j = 0; j < N_TOPK; j++) {
                s -= topk[(size_t)k * N_TOPK + j];
                double m = s / (nv - j - 1);
                egsInformation("  %+9.3f", 100.0 * (m - mean) / mean);
            }
            egsInformation("\n");
        }
        egsInformation("\n  Percent shift in BUF when the largest 1, 2 ... %d bunches are\n"
                       "  removed.  Shifts of order 1/n are what an untailed distribution\n"
                       "  gives; anything far larger means the shell is carried by those\n"
                       "  few bunches and its sigma is optimistic however many were run.\n",
                       N_TOPK);
    }

    // ---- Between-job (batch) variance --------------------------------------
    //
    // Same quantity, coarser batching: sigma_bunch pools n_valid_b bunches,
    // sigma_job pools J job means of ~n_valid_b/J bunches each.  With
    //
    //     sigma_b^2 estimated by  SS/(J-1),  SS = sum_j n_j (rbar_j - rbar)^2
    //                                           = sum_j n_j rbar_j^2 - N rbar^2
    //
    // and Var(rbar) = sigma_b^2 / N, the two agree when the per-bunch variance
    // is actually determined.  They diverge when it is not -- which is the
    // point of printing them side by side.
    //
    // Read the ratio, not either column alone, and mind two things:
    //   * sigma on sigma is 1/sqrt(2(J-1)) -- 2.5% at J=800.  A ratio inside
    //     ~0.9-1.1 is agreement; do not read structure into it.
    //   * BOTH are biased low under a heavy tail.  Agreement rules out one
    //     specific failure; it does not demonstrate convergence.  Same trap as
    //     N_eff: necessary, not sufficient.
    if (n_jobs_comb >= 2) {
        egsInformation("\n  Between-job (batch) vs between-bunch variance"
                       "   [%d jobs combined]\n", n_jobs_comb);
        egsInformation("  %-6s  %-9s  %-8s  %-12s  %-12s  %-9s\n",
                       "Region", "eta/mfp", "jobs", "sigma_bunch%", "sigma_job%",
                       "job/bunch");
        egsInformation("  %s\n", string(66, '-').c_str());
        for (int k = 0; k < n_scoring; k++) {
            int nv = n_valid_b[k], J = job_used[k];
            if (nv < 2 || J < 2 || sum_ratio2[k] <= 0.0) continue;
            double mean = sum_ratio[k] / nv;
            if (mean <= 0.0) continue;
            double eta  = (sum_Kp[k] > 0.0 && sum_Kp0[k] > 0.0)
                          ? -std::log(sum_Kp[k] / sum_Kp0[k]) : 0.0;

            double var_b = (sum_ratio2[k] - sum_ratio[k] * sum_ratio[k] / nv)
                           / (nv - 1);
            double s_bunch = std::sqrt(std::max(var_b, 0.0) / nv) / mean;

            double ss = job_nr2[k] - sum_ratio[k] * sum_ratio[k] / nv;
            double var_j = std::max(ss, 0.0) / (J - 1);
            double s_job = std::sqrt(var_j / nv) / mean;

            egsInformation("  %-6d  %-9.4f  %-8d  %-12.4f  %-12.4f  %-9.2f\n",
                           scoring_reg[k], eta, J,
                           100.0 * s_bunch, 100.0 * s_job,
                           s_bunch > 0 ? s_job / s_bunch : 0.0);
        }
        egsInformation("\n  Two batchings of the same data, not independent estimators.\n"
                       "  Agreement (ratio ~1) is consistent with a converged variance;\n"
                       "  a ratio well above 1 means the per-bunch estimate is missing\n"
                       "  tail mass and the quoted sigma is too small.  Both are biased\n"
                       "  low under a heavy tail, so agreement is necessary, not sufficient.\n"
                       "  Uncertainty on sigma_job itself is 1/sqrt(2(J-1)) = %.1f%%.\n",
                       100.0 / std::sqrt(2.0 * std::max(n_jobs_comb - 1, 1)));
    }

    const bool is_combined_output = (getNparallel() > 0 && getIparallel() == 0);

    // ---- Forced-collision weight audit -------------------------------------
    //
    // Per stage the split is exact: w = w*T + w*(1-T).  What is NOT guaranteed
    // exact is the hand-off between stages.  A primary banked at surface m is
    // dropped into a container, combed with everything else there, and replayed;
    // comb() preserves total weight only in EXPECTATION.  So
    //
    //     w_in[m]  should equal  w_bank[m-1]
    //
    // and the last column is the whole point of this table.  A value near 1 at
    // every surface clears the bank -> comb -> replay path; a value above 1 that
    // compounds is the leak, and its size per stage should reproduce the ~1.1 to
    // 1.25 per-stage excess seen in BUF at 25-50 mfp.
    //
    // w_in[0] is a special case: those primaries come straight from the source,
    // not from a container, so its ratio column is blank.
    if (forced_collision && pri_diag && !is_combined_output) {
        egsInformation("\n  Forced-collision weight audit"
                       "   [cumulative, this job]\n");
        egsInformation("  %-6s  %-6s  %-12s  %-12s  %-12s  %-12s  %-10s\n",
                       "surf", "region", "splits", "w_in", "w_bank", "w_coll",
                       "w_in/w_bank(prev)");
        egsInformation("  %s\n", string(88, '-').c_str());
        double coll_tot = 0, bank_last = 0;
        for (int m = 0; m < n_combing; m++) {
            coll_tot += fc_coll_w[m];
            egsInformation("  %-6d  %-6d  %-12lld  %-12.6g  %-12.6g  %-12.6g",
                           m, combing_reg[m], fc_n_split[m],
                           fc_in_w[m], fc_bank_w[m], fc_coll_w[m]);
            if (m == 0) egsInformation("  %-10s\n", "(source)");
            else if (fc_bank_w[m-1] > 0)
                egsInformation("  %-10.4f\n", fc_in_w[m] / fc_bank_w[m-1]);
            else egsInformation("  %-10s\n", "-");
            bank_last = fc_bank_w[m];
        }
        egsInformation("\n  Primary weight accounted for: collided %.6g + banked at last"
                       " surface %.6g\n  + fell through to analog %.6g (%lld free paths)"
                       "  =  %.6g\n",
                       coll_tot, bank_last, fc_analog_w, fc_analog_n,
                       coll_tot + bank_last + fc_analog_w);
        egsInformation("  Source weight this job: %.6g"
                       "   (equal only if no primary escaped the geometry)\n",
                       (double)n_completed * (double)n_per_bunch);
        egsInformation("\n  w_in[m] should equal w_bank[m-1]: comb() and replay preserve\n"
                       "  weight in expectation, so a ratio that sits above 1 and\n"
                       "  compounds per stage is a leak in bank -> comb -> replay, not\n"
                       "  in the split.\n");
    }

    // ---- Scattered-channel comb()/replay weight ledger ---------------------
    //
    // The primary-channel audit above traces one specific quantity (the staged
    // primary split) and has come back exact in three independent placements.
    // It says nothing about the ordinary scattered population, which passes
    // through comb() and replayContainer() far more often -- hundreds to
    // thousands of times per job -- and is the one part of the cascade the
    // primary audit cannot see, since a replayed particle is an ordinary
    // scattered photon by the time it is transported again.
    //
    // comb() is unbiased per particle by construction: E[ncopies_j] = w_j/wbar
    // exactly, so the cumulative after/before ratio should approach 1 as
    // n_comb_calls grows, with per-call noise of order 1/sqrt(comb_target).  A
    // ratio that does NOT approach 1 -- especially one that differs between
    // placements with different comb-call densities -- would be a real bug in
    // comb() or a call site, not sampling noise, and is a live candidate for
    // the A/B/C disagreement now that the primary ledger is cleared.
    //
    // Runs whenever 'primary crossing diagnostic = yes', independent of
    // 'forced collision', so the same table is available for the FC-off
    // control arm -- comb() and replay exist regardless of forced collision.
    if (pri_diag && n_completed > 0 && !is_combined_output) {
        egsInformation("\n  Scattered-channel comb()/replay weight ledger"
                       "   [cumulative, this job]\n");
        egsInformation("  %-6s  %-6s  %-10s  %-13s  %-13s  %-13s  %-13s  %-13s\n",
                       "surf", "region", "comb calls", "bank(scat)", "comb_before",
                       "comb_after", "after/before", "replay_in");
        egsInformation("  %s\n", string(96, '-').c_str());
        for (int m = 0; m < n_combing; m++) {
            double ratio = comb_before_w[m] > 0
                          ? comb_after_w[m] / comb_before_w[m] : 0.0;
            egsInformation("  %-6d  %-6d  %-10lld  %-13.6g  %-13.6g  %-13.6g  %-13.7f  %-13.6g\n",
                           m, combing_reg[m], n_comb_calls[m],
                           bank_wsum[m] / n_completed,
                           comb_before_w[m], comb_after_w[m], ratio,
                           replay_in_w[m]);
        }
        egsInformation("\n  after/before should approach 1 as comb calls accumulate; the\n"
                       "  per-call noise scale is roughly 1/sqrt(comb_target), so with\n"
                       "  comb_target ~1e5 a single call fluctuates by ~0.3%%, and the\n"
                       "  CUMULATIVE ratio over N calls should tighten toward 1 by\n"
                       "  roughly 1/sqrt(N).  A ratio that stays away from 1 as N grows\n"
                       "  is not noise.\n"
                       "  replay_in is the post-comb weight actually handed to shower():\n"
                       "  the one quantity in this cascade the primary-channel audit\n"
                       "  above cannot see at all.\n");
    }

    // ---- Primary-crossing diagnostic ---------------------------------------
    // Suppressed on the final combine: max_ratio[] is merged across jobs but
    // max_deep[] is not (it is not in the .egsdat), so on a combine the two
    // columns would describe different bunches.  Print only where this process
    // ran the bunches itself -- a serial run, or an individual parallel job.
    if (pri_diag && n_completed > 0 && !is_combined_output) {
        egsInformation("\n  Primary crossings of the combing surfaces"
                       "   [this job only, %d bunches]\n", n_completed);
        egsInformation("  %-8s  %-6s  %-14s  %-13s  %-13s  %-9s  %-9s\n",
                       "Surface", "region", "crossings/bnch", "pri wt/bnch",
                       "bank wt/bnch", "pri frac", "bnch deep");
        egsInformation("  %s\n", string(86, '-').c_str());
        for (int m = 0; m < n_combing; m++) {
            double pw = pri_cross_wsum[m] / n_completed;
            double bw = bank_wsum[m]      / n_completed;
            egsInformation("  %-8d  %-6d  %-14.6g  %-13.6g  %-13.6g  %-9.4f  %-9d\n",
                           m, combing_reg[m],
                           pri_cross_sum[m] / n_completed, pw, bw,
                           (pw + bw) > 0 ? pw / (pw + bw) : 0.0,
                           deep_hist[m + 1]);
        }
        egsInformation("  %-8s  %-6s  %-14s  %-13s  %-13s  %-9s  %-9d\n",
                       "none", "-", "-", "-", "-", "-", deep_hist[0]);

        egsInformation("\n  Deepest surface reached by any primary, in the bunch holding\n"
                       "  each shell's largest K_tot/K_pri:\n");
        egsInformation("  %-6s  %-9s  %-10s  %-10s  %-12s\n",
                       "Region", "eta/mfp", "max/mean", "max bunch", "deepest surf");
        egsInformation("  %s\n", string(56, '-').c_str());
        for (int k = 0; k < n_scoring; k++) {
            if (n_valid_b[k] < 2) continue;
            double mean = sum_ratio[k] / n_valid_b[k];
            double eta  = (sum_Kp[k] > 0.0 && sum_Kp0[k] > 0.0)
                          ? -std::log(sum_Kp[k] / sum_Kp0[k]) : 0.0;
            egsInformation("  %-6d  %-9.4f  %-10.3f  %-10.6g  %-12d\n",
                           scoring_reg[k], eta,
                           max_ratio[k] / std::max(mean, 1e-30), max_ratio[k],
                           max_deep[k]);
        }
        egsInformation("\n  A primary is neither banked at a combing surface nor killed\n"
                       "  there; it transports analog, surviving with probability\n"
                       "  exp(-spacing) and carrying weight 1 into a population whose\n"
                       "  weight has decayed by the same factor.\n"
                       "\n"
                       "  'pri frac' is the decisive column and needs no rare event to\n"
                       "  be observed: it is the share of the weight arriving at a\n"
                       "  surface that is carried by analog primaries rather than by\n"
                       "  banked scatter.  A share that grows with depth means the\n"
                       "  primary channel dominates there, and since it is sampled\n"
                       "  analog at probability exp(-spacing) it is then also the\n"
                       "  dominant variance source.\n"
                       "\n"
                       "  'bnch deep'/'deepest surf' are the direct but expensive test:\n"
                       "  they need enough bunches for a primary to actually reach a\n"
                       "  deep surface -- at 10 mfp spacing that is one bunch in ~500\n"
                       "  for surface 1 and far rarer beyond.  If 'deepest surf' runs\n"
                       "  high exactly where 'max/mean' does, the rare high-weight\n"
                       "  event is confirmed as the tail's source.  -1 means no primary\n"
                       "  reached any surface.\n"
                       "\n"
                       "  Per-job only: this is not written to the .egsdat.\n");
    }

    // ---- Speed and efficiency summary ----
    long long total_source = (long long)n_completed * n_per_bunch;
    double speed = (total_cpu_time_ > 0)
                   ? total_source / total_cpu_time_ : 0.0;
    egsInformation("\n");
    // Dump every per-bunch ratio for offline analysis.  Written on a serial run
    // or the final combine only; an individual parallel job's 8 bunches are not
    // worth a file.
    if (bunch_stats && !(getNparallel() > 0 && getIparallel() > 0)) {
        size_t nb = 0;
        for (int k = 0; k < n_scoring; k++) nb += bunch_r[k].size();
        if (nb > 0) {
            string fn = egsJoinPath(getAppDir(),
                                    getFinalOutputFile() + ".bunchdat");
            ofstream f(fn.c_str());
            if (f) {
                f << "# per-bunch K_tot/K_pri.  one block per shell.\n";
                f << "# region eta n_bunches, then the values\n";
                for (int k = 0; k < n_scoring; k++) {
                    double eta = (sum_Kp[k] > 0.0 && sum_Kp0[k] > 0.0)
                                 ? -std::log(sum_Kp[k] / sum_Kp0[k]) : 0.0;
                    f << scoring_reg[k] << " " << eta << " "
                      << bunch_r[k].size() << "\n";
                    for (size_t m = 0; m < bunch_r[k].size(); m++)
                        f << bunch_r[k][m] << (((m+1) % 10) ? " " : "\n");
                    f << "\n";
                }
                egsInformation("\n  Per-bunch ratios written to %s\n"
                               "  (%lu values; sigma, N_eff, tail index and bootstrap\n"
                               "  intervals are all recoverable from this offline)\n",
                               fn.c_str(), (unsigned long)nb);
            }
        }
    }

    egsInformation("  Bunches: %d of %d  |  Source photons: %lld  |  CPU: %.1f s"
                   "  |  Speed: %.0f photons/s\n",
                   n_completed, n_bunches, total_source, total_cpu_time_, speed);
    egsInformation("  eta = -ln(K_pri/K_pri_unattenuated), the transmission-weighted\n"
                   "  number of mean free paths along the forced-detection rays.\n");
    egsInformation("  FOM = 1/(sigma_BUF^2 * T);  higher is better.  Expect a plateau in\n"
                   "  the deep, cascade-driven shells; it falls steeply before that.\n");
    egsInformation("\n");
}


/*----------------------------------------------------------------------------
  finishSimulation

  The .egsdat needed for URCO parallel combining is written at the end of
  runSimulation() (see the comment there) — every job writes one, including
  the watcher, so every job's results reach the combine.

  The watcher (last job by URCO convention) is handled by
  EGS_AdvancedApplication, which calls combineResults() when
  run->finishSimulation() returns 1.

  JCF fallback (err < 0): if the filesystem doesn't support file locking,
  JCFControl returns -2.  We recover manually — same logic as before.
----------------------------------------------------------------------------*/
int EGS_ShieldApplication::finishSimulation() {
    int np = getNparallel();
    int ip = getIparallel();
    int err = EGS_AdvancedApplication::finishSimulation();

    // Fallback: if JCF failed (err < 0), force the correct path manually.
    if (err < 0 && np > 0) {
        if (ip > 0) {
            outputResults();
            err = outputData();
        } else {
            combineResults();
            err = 0;
        }
    }
    return err;
}


/*----------------------------------------------------------------------------
  Parallel combining: outputData / readData / addState / resetCounter

  After each parallel job completes, the framework writes a .egsdat file via
  outputData().  combineResults() then resets the accumulator and calls
  addState() once per job file.  Only the sufficient statistics are needed:
  sum_ratio[k], sum_ratio2[k], sum_Kt[k], sum_Kp[k], n_valid_b[k], n_completed.
----------------------------------------------------------------------------*/
int EGS_ShieldApplication::outputData() {
    int err = EGS_AdvancedApplication::outputData();
    if (err) return err;
    (*data_out) << n_completed << endl << total_cpu_time_ << endl;
    for (int k = 0; k < n_scoring; k++)
        (*data_out) << sum_ratio[k]  << " " << sum_ratio2[k] << " "
                    << sum_Kt[k]     << " " << sum_Kp[k]     << " "
                    << sum_Kp0[k]    << " " << max_ratio[k]  << " "
                    << n_valid_b[k]  << endl;

    // Top-k block, written LAST and read optionally, so that a file produced
    // here still combines with a binary that predates it and -- more usefully
    // -- so the 6400-bunch datasets already on disk remain readable.  Appending
    // fields to the per-shell lines instead would have made an old file's
    // shells misalign, which is how the sum_Kp0 addition broke compatibility on
    // 2026-08-10.
    (*data_out) << "TOPK " << N_TOPK << endl;
    for (int k = 0; k < n_scoring; k++) {
        for (int j = 0; j < N_TOPK; j++)
            (*data_out) << topk[(size_t)k * N_TOPK + j] << " ";
        (*data_out) << endl;
    }
    (*data_out) << "BUNCHR" << endl;
    for (int k = 0; k < n_scoring; k++) {
        (*data_out) << bunch_r[k].size();
        for (size_t j = 0; j < bunch_r[k].size(); j++)
            (*data_out) << " " << bunch_r[k][j];
        (*data_out) << endl;
    }
    return (*data_out) ? 0 : 99;
}

int EGS_ShieldApplication::readData() {
    int err = EGS_AdvancedApplication::readData();
    if (err) return err;
    (*data_in) >> n_completed >> total_cpu_time_;
    for (int k = 0; k < n_scoring; k++)
        (*data_in) >> sum_ratio[k]  >> sum_ratio2[k]
                   >> sum_Kt[k]     >> sum_Kp[k]
                   >> sum_Kp0[k]    >> max_ratio[k]
                   >> n_valid_b[k];
    return (*data_in) ? 0 : 99;
}

int EGS_ShieldApplication::addState(istream &data) {
    int err = EGS_AdvancedApplication::addState(data);
    if (err) return err;
    int    nc; double tcpu;
    data >> nc >> tcpu;
    n_completed     += nc;
    total_cpu_time_ += tcpu;
    for (int k = 0; k < n_scoring; k++) {
        double sr, sr2, skt, skp, skp0, mx; int nv;
        data >> sr >> sr2 >> skt >> skp >> skp0 >> mx >> nv;
        sum_ratio[k]  += sr;
        sum_ratio2[k] += sr2;
        sum_Kt[k]     += skt;
        sum_Kp[k]     += skp;
        sum_Kp0[k]    += skp0;
        if (mx > max_ratio[k]) max_ratio[k] = mx;   // max, not sum
        n_valid_b[k]  += nv;
    }
    if (!data) return 99;

    // Optional top-k block.  Absent from files written before 2026-08-13; a
    // clean EOF there is not an error, so clear the failbit and carry on with
    // an empty top-k rather than rejecting the file.
    string tag;
    if (data >> tag && tag == "TOPK") {
        int kk = 0;
        data >> kk;
        if (!data || kk < 1) return 99;
        for (int k = 0; k < n_scoring; k++)
            for (int j = 0; j < kk; j++) {
                double v = 0;
                data >> v;
                if (!data) return 99;
                insertTopK(k, v);   // merges sorted lists: exact global top-k
            }
    }
    else {
        data.clear();
        return 0;      // pre-2026-08-13 file: no TOPK, so no BUNCHR either
    }

    string tag2;
    if (data >> tag2 && tag2 == "BUNCHR") {
        for (int k = 0; k < n_scoring; k++) {
            size_t nb = 0;
            data >> nb;
            if (!data) return 99;
            for (size_t j = 0; j < nb; j++) {
                double v = 0;
                data >> v;
                if (!data) return 99;
                bunch_r[k].push_back(v);
            }
        }
    }
    else {
        data.clear();
    }
    return 0;
}

/*----------------------------------------------------------------------------
  combineResults — like the base class, but a bad .egsdat is skipped instead
  of failing the whole combine.

  EGS_Application::combineResults() returns -1 if any single file fails to
  parse, and EGS_AdvancedApplication::finishSimulation() bails on that before
  outputResults() ever runs.  One truncated file therefore costs the entire
  run's output.  That is a poor trade at 800 jobs, and more likely now that
  outputData() runs once per bunch: a job killed during a write leaves a short
  file behind.

  Failed files are rolled back and skipped.  The rollback covers this class's
  accumulators, which are the only inputs to the reported results; a partial
  read has already perturbed the base class's run/rndm/source state, which we
  cannot undo, but that affects only the cosmetic ncase and cpu columns of the
  listing below, not BUF, K or sigma.
----------------------------------------------------------------------------*/
int EGS_ShieldApplication::combineResults() {
    // A combine invoked without -P leaves n_parallel at 0.  That is a supported
    // way to run one (egs_shield -i <input> with 'calculation = combine'), and
    // EGS_Application::combineResults handles it by substituting
    // MAXIMUM_JOB_NUMBER and scanning for whatever files exist
    // (egs_application.cpp:640).  Do the same here.
    //
    // Delegating to the base class instead -- which is what this did until
    // 2026-08-13 -- produced a log that was indistinguishable from a good one:
    // the base prints the same banner and the same per-file lines, so an
    // 800-job combine listed all 800 files and simply had no batch-statistics
    // table.  Nothing indicated that a different function had done the work.
    static const int MAX_JOBS_SCAN = 8192;   // MAXIMUM_JOB_NUMBER, not exported
    int np = getNparallel();
    if (np <= 0) {
        np = MAX_JOBS_SCAN;
    }

    egsInformation(
        "\n                      Suming the following .egsdat files:\n"
        "=======================================================================\n");
    resetCounter();
    n_jobs_comb = 0;
    std::fill(job_nr2.begin(),  job_nr2.end(),  0.0);
    std::fill(job_used.begin(), job_used.end(), 0);

    // Snapshots for rollback.
    vector<double> s_ratio(n_scoring), s_ratio2(n_scoring), s_Kt(n_scoring),
           s_Kp(n_scoring), s_Kp0(n_scoring), s_max(n_scoring);
    vector<int>    s_valid(n_scoring);

    EGS_Float last_cpu = 0;
    EGS_I64   last_ncase = 0;
    int ndat = 0, nbad = 0;
    char buf[512];

    for (int j = getFirstParallel(); j < getFirstParallel() + np; j++) {
        sprintf(buf, "%s_w%d.egsdat", getFinalOutputFile().c_str(), j);
        string dfile = egsJoinPath(getAppDir(), buf);
        ifstream data(dfile.c_str());
        if (!data) {
            continue;    // job never ran, or never got far enough to write
        }

        int    save_completed = n_completed;
        double save_cpu       = total_cpu_time_;
        for (int k = 0; k < n_scoring; k++) {
            s_ratio[k]  = sum_ratio[k];
            s_ratio2[k] = sum_ratio2[k];
            s_Kt[k]     = sum_Kt[k];
            s_Kp[k]     = sum_Kp[k];
            s_Kp0[k]    = sum_Kp0[k];
            s_max[k]    = max_ratio[k];
            s_valid[k]  = n_valid_b[k];
        }

        int err = addState(data);
        if (err) {
            n_completed     = save_completed;
            total_cpu_time_ = save_cpu;
            for (int k = 0; k < n_scoring; k++) {
                sum_ratio[k]  = s_ratio[k];
                sum_ratio2[k] = s_ratio2[k];
                sum_Kt[k]     = s_Kt[k];
                sum_Kp[k]     = s_Kp[k];
                sum_Kp0[k]    = s_Kp0[k];
                max_ratio[k]  = s_max[k];
                n_valid_b[k]  = s_valid[k];
            }
            ++nbad;
            egsWarning("   %-30s SKIPPED (parse error %d)\n", buf, err);
            continue;
        }

        ++ndat;

        // Between-job batch statistics.  The snapshot taken above for rollback
        // doubles as the "before" state, so job j's own contribution is a
        // subtraction and costs nothing.  rbar_j = dS/dn is that job's mean
        // ratio; accumulating n_j * rbar_j^2 = dS^2/dn gives the weighted
        // between-job sum of squares later.  Weighting by dn matters because
        // egs-parallel gives the first (n_bunches % np) jobs one extra bunch.
        ++n_jobs_comb;
        for (int k = 0; k < n_scoring; k++) {
            int dn = n_valid_b[k] - s_valid[k];
            if (dn > 0) {
                double dS = sum_ratio[k] - s_ratio[k];
                job_nr2[k] += dS * dS / dn;
                job_used[k]++;
            }
        }

        EGS_I64   ncase = run->getNdone();
        EGS_Float cpu   = run->getCPUTime();
        egsInformation("%2d %-30s ncase=%-14lld cpu=%-11.2f\n",
                       ndat, buf, ncase - last_ncase, cpu - last_cpu);
        last_ncase = ncase;
        last_cpu   = cpu;
    }

    if (ndat > 0) {
        egsInformation(
            "=======================================================================\n");
        egsInformation("%40s%-14lld cpu=%-11.2f\n\n", "Total ncase=",
                       last_ncase, last_cpu);
    }
    if (nbad > 0) {
        egsWarning("\n  *** %d .egsdat file(s) were unreadable and have been "
                   "skipped.\n"
                   "  *** Their contribution is absent from the results below;\n"
                   "  *** the bunch count in the summary line reflects this.\n\n",
                   nbad);
    }
    // Unlike the base class, succeed as long as something was read.
    return ndat > 0 ? 0 : 1;
}


void EGS_ShieldApplication::resetCounter() {
    EGS_AdvancedApplication::resetCounter();
    n_completed     = 0;
    total_cpu_time_ = 0.0;
    for (int k = 0; k < n_scoring; k++) {
        sum_ratio[k]  = 0.0;
        sum_ratio2[k] = 0.0;
        sum_Kt[k]     = 0.0;
        sum_Kp[k]     = 0.0;
        sum_Kp0[k]    = 0.0;
        max_ratio[k]  = 0.0;
        n_valid_b[k]  = 0;
    }
    std::fill(topk.begin(), topk.end(), 0.0);
    for (size_t k = 0; k < bunch_r.size(); k++) bunch_r[k].clear();
}


/*----------------------------------------------------------------------------
  Destructor
----------------------------------------------------------------------------*/
EGS_ShieldApplication::~EGS_ShieldApplication() {
    delete E_Muen_Rho;

    delete [] scoring_reg;
    delete [] scoring_idx;
    delete [] is_scoring;
    delete [] V_shell;

    delete [] combing_reg;
    delete [] combing_idx;
    delete [] is_combing;

    if (containers) {
        for (int m = 0; m < n_combing; m++) delete containers[m];
        delete [] containers;
    }
    if (snapshots) {
        for (int m = 0; m < n_combing; m++) delete snapshots[m];
        delete [] snapshots;
    }
    delete [] sum_ratio;
    delete [] sum_ratio2;
    delete [] sum_Kt;
    delete [] sum_Kp;
    delete [] sum_Kp0;
    delete [] max_ratio;
    delete [] n_valid_b;
}


/*----------------------------------------------------------------------------
  Fortran interface: select_photon_mfp — called once per photon free path.
  Routed here via the $SELECT-PHOTON-MFP override in egs_shield.macros.
----------------------------------------------------------------------------*/
extern __extc__ void
F77_OBJ_(select_photon_mfp,SELECT_PHOTON_MFP)(EGS_Float *dpmfp) {
    EGS_Application *a = EGS_Application::activeApplication();
    EGS_ShieldApplication *app =
        dynamic_cast<EGS_ShieldApplication *>(a);
    if (!app)
        egsFatal("select_photon_mfp: active application is not "
                 "EGS_ShieldApplication\n");
    app->selectPhotonMFP(*dpmfp);
}


#ifdef BUILD_APP_LIB
    APP_LIB(EGS_ShieldApplication);
#else
    APP_MAIN(EGS_ShieldApplication);
#endif
