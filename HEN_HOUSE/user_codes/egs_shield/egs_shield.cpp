
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
          cascade_diag(false),
          n_bunches(100), n_per_bunch(1000000LL),
          current_bunch(0), n_completed(0),
          sum_ratio(nullptr), sum_ratio2(nullptr),
          sum_Kt(nullptr), sum_Kp(nullptr), sum_Kp0(nullptr),
          n_valid_b(nullptr),
          total_cpu_time_(0.0)
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
        if (ir >= 0 && ir < nreg && is_combing[ir] && latch != 0) {
            the_stack->wt[np]    = 0;
            the_epcont->idisc = -1;
        }
        return 0;
    }

    /*------------------------------------------------------------------------
      selectPhotonMFP — called by Fortran at the start of each photon free
      path via the $SELECT-PHOTON-MFP macro override in egs_shield.macros.
    ------------------------------------------------------------------------*/
    void selectPhotonMFP(EGS_Float &dpmfp) {
        bool primary = (the_stack->latch[the_stack->np - 1] == 0);
        scoreFD_all(primary);
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
    int       n_bunches;
    long long n_per_bunch;
    int       current_bunch;
    int       n_completed;     // bunches scored so far (all jobs combined)
    // Sufficient statistics for combining parallel jobs (all [n_scoring]):
    double   *sum_ratio;       // Σ (K_tot_b / K_pri_b) over valid bunches
    double   *sum_ratio2;      // Σ (K_tot_b / K_pri_b)² over valid bunches
    double   *sum_Kt;          // Σ K_tot_b
    double   *sum_Kp;          // Σ K_pri_b
    double   *sum_Kp0;         // Σ K_pri_b computed with exp(-Λ) omitted;
                               // K_pri/K_pri0 is the primary transmission, so
                               // η_eff = -ln(ΣK_pri / ΣK_pri0).  See outputResults().
    int      *n_valid_b;       // # bunches with K_pri_b > 0
    double    total_cpu_time_; // accumulated CPU time over all bunches [s]


    /*------------------------------------------------------------------------
      scoreFD_all

      Ray-trace from the current photon position along direction u,
      accumulating optical depth Lambda.

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
    int scoreFD_all(bool is_primary);

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
int EGS_ShieldApplication::scoreFD_all(bool is_primary) {
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

            if (t_shell < TSTEP_MAX) {
                EGS_Float emuen_rho = E_Muen_Rho->interpolateFast(gle);
                EGS_Float unatt = wt * emuen_rho * t_shell / V_shell[k];
                EGS_Float score = unatt * std::exp(-Lambda);
                sum_tot[k] += score;               // total = primary + scattered
                if (is_primary) {
                    sum_pri[k]  += score;
                    // Same score with the attenuation factor removed.  The
                    // ratio of the two sums is the primary transmission along
                    // the FD rays, from which outputResults() recovers the
                    // effective optical depth.  Costs one multiply-add per
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
        // Fires whenever a scattered photon reaches any combing region,
        // whether or not that region is also a scoring shell.
        // The virtual photon is placed at the combing region's EXIT so the
        // replayed particle starts in the correct medium on the far side.
        if (!is_primary && is_combing[ireg]) {
            EGS_Float t_comb = TSTEP_MAX;
            int       med_after;
            int ireg_after = geometry->howfar(ireg, x, u, t_comb, &med_after);
            // Include attenuation through the combing shell itself.
            // sigma here is for the medium just entered (same iron as surroundings).
            if (t_comb < TSTEP_MAX) Lambda += t_comb * sigma;
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

    sum_ratio  = new double[n_scoring]();
    sum_ratio2 = new double[n_scoring]();
    sum_Kt     = new double[n_scoring]();
    sum_Kp     = new double[n_scoring]();
    sum_Kp0    = new double[n_scoring]();
    n_valid_b  = new int[n_scoring]();
    n_completed = 0;

    // Enable ausgab calls needed for latch bookkeeping.
    // By default the framework only enables BeforeTransport..AfterTransport.
    // AfterCompton/Photo/Pair/Rayleigh are above AfterTransport in the enum
    // and must be turned on explicitly so that scattered photons get latch > 0.
    setAusgabCall(AfterCompton,  true);
    setAusgabCall(AfterPhoto,    true);
    setAusgabCall(AfterPair,     true);
    setAusgabCall(AfterRayleigh, true);

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

        // ---- Stage 0: source photons ----
        // Periodic combing every comb_interval photons keeps containers from
        // growing without bound.  We cap the population at comb_interval×10 so
        // that the container never grows much larger than that ceiling.
        const long long comb_interval = std::max(10000LL, n_per_bunch / 100);
        const int       comb_target   = (int)(comb_interval * 10);

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
                for (int m = 0; m < n_combing; m++)
                    containers[m]->comb(rndm, comb_target);
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
        // where Λ_min is the optical depth to the nearest combing shell.
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
        for (int iter = 0; iter < MAX_ITER; iter++) {
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

            if (w_cur <= 0 || w_cur < 1e-10 * w_initial) {
                if (diag)
                    egsInformation("  [cascade] STOP at iter %d: w_cur/w_ini = %.3g"
                                   " < 1e-10 threshold\n\n", iter,
                                   w_initial > 0 ? w_cur / w_initial : 0.0);
                break;
            }

            for (int m = 0; m < n_combing; m++) {
                std::swap(containers[m], snapshots[m]);
                containers[m]->clean();
                // Cap at current size: don't split when np < comb_target.
                // Splitting would create a fixed-cost 100K-shower replay
                // regardless of n_per_bunch, breaking time scaling.
                snapshots[m]->comb(rndm, std::min(snapshots[m]->size(), comb_target));
            }
            for (int m = 0; m < n_combing; m++)
                if (snapshots[m]->size() > 0)
                    replayContainer(snapshots[m]);
        }

        double bunch_time = timer.time();
        total_cpu_time_ += bunch_time;
        endBunch();
        egsInformation("  %.1f s\n", bunch_time);
    }

    // We replace the standard batch loop, so EGS_RunControl::finishBatch() is
    // never reached.  Call it once here, at the point the batch loop would
    // have, to do the end-of-batch bookkeeping: it sets the run control's
    // cpu_time and then calls outputData() for us.  Doing this for EVERY
    // parallel job -- the watcher included -- is what gets every job's results
    // into the combine.
    //
    // It must happen here and not in finishSimulation(): for the watcher the
    // whole combine runs inside EGS_AdvancedApplication::finishSimulation(),
    // and the partial combine inside URCO's poll loop calls resetCounter(),
    // which zeroes our accumulators.  Writing afterwards would store zeros
    // into a file nothing reads any more.
    //
    // run_dir is still set at this point, so the file lands in this job's
    // egsrun_* directory: invisible to howManyJobsDone() (which scans the
    // application directory only), so the watcher's poll loop still waits for
    // npar-1 *other* jobs.  finishRun() then moves it up to the application
    // directory before combineResults() runs.
    //
    // finishBatch() also prints a batch summary line and applies the
    // "statistical accuracy sought" early-termination test.  We do not
    // override getCurrentResult(), so that line's result/uncertainty columns
    // print as 0 and 100.00; the return value is only meaningful for loop
    // control, which we do not have, so it is discarded.
    if (getNparallel() > 0 && getIparallel() > 0) {
        run->finishBatch();
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
            n_valid_b[k]++;
        }
    }
    n_completed++;
}


/*----------------------------------------------------------------------------
  outputResults

  Effective optical depth (eta/mfp column)
  ---------------------------------------
  The FD estimator already carries the exact quantity we want.  Every primary
  score is w·exp(-Λ)·(E μen/ρ)·t/V, where Λ is the optical depth accumulated
  along the ray from the source to the shell.  Accumulating the same score
  without exp(-Λ) (sum_Kp0) makes the ratio ΣK_pri/ΣK_pri0 the primary
  transmission, so

      η_eff = -ln( ΣK_pri / ΣK_pri0 )

  For a monoenergetic point source and concentric shells every primary ray to
  a given shell is identical, so η_eff is exactly that shell's optical depth,
  to machine precision and with no user input.

  For a polyenergetic source the per-photon η varies with energy and η_eff
  becomes -ln⟨exp(-η)⟩, the transmission-weighted mean — which is the right
  characteristic depth here, because it is precisely the attenuation that
  forms the denominator of the buildup factor.  It is not the arithmetic mean
  optical depth, and for a broad spectrum it will sit below it (low-energy
  components are attenuated away and stop contributing).  The same expression
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

    // ---- Speed and efficiency summary ----
    long long total_source = (long long)n_completed * n_per_bunch;
    double speed = (total_cpu_time_ > 0)
                   ? total_source / total_cpu_time_ : 0.0;
    egsInformation("\n");
    egsInformation("  Bunches: %d of %d  |  Source photons: %lld  |  CPU: %.1f s"
                   "  |  Speed: %.0f photons/s\n",
                   n_completed, n_bunches, total_source, total_cpu_time_, speed);
    egsInformation("  eta = -ln(K_pri/K_pri_unattenuated), the transmission-weighted\n"
                   "  optical depth along the forced-detection rays.\n");
    egsInformation("  FOM = 1/(sigma_BUF^2 * T);  higher is better."
                   "  Expect ~constant with depth for this algorithm.\n");
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
                    << sum_Kp0[k]    << " "
                    << n_valid_b[k]  << endl;
    return (*data_out) ? 0 : 99;
}

int EGS_ShieldApplication::readData() {
    int err = EGS_AdvancedApplication::readData();
    if (err) return err;
    (*data_in) >> n_completed >> total_cpu_time_;
    for (int k = 0; k < n_scoring; k++)
        (*data_in) >> sum_ratio[k]  >> sum_ratio2[k]
                   >> sum_Kt[k]     >> sum_Kp[k]
                   >> sum_Kp0[k]
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
        double sr, sr2, skt, skp, skp0; int nv;
        data >> sr >> sr2 >> skt >> skp >> skp0 >> nv;
        sum_ratio[k]  += sr;
        sum_ratio2[k] += sr2;
        sum_Kt[k]     += skt;
        sum_Kp[k]     += skp;
        sum_Kp0[k]    += skp0;
        n_valid_b[k]  += nv;
    }
    return data ? 0 : 99;
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
        n_valid_b[k]  = 0;
    }
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
