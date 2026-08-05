
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
#  Algorithm (Divide et Impera, extended for backscatter):
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
#        combing regions   = 9 19 29 ...     # subset of scoring region indices
#        emuen file        = $EGS_HOME/egs_kerma/emuen_rho_air_1keV-20MeV.data
#        n bunches         = 100
#        photons per bunch = 1000000
#    :stop scoring options:
#
###############################################################################
*/

#include "egs_advanced_application.h"
#include "egs_scoring.h"
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
          kerma_tot(nullptr), kerma_pri(nullptr),
          E_Muen_Rho(nullptr),
          n_bunches(100), n_per_bunch(1000000LL),
          current_bunch(0),
          K_bunch(nullptr), Kp_bunch(nullptr),
          last_case(0)
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

    int  initScoring()   override;
    int  runSimulation() override;
    void outputResults() override;

protected:
    int startNewShower() override {
        int res = EGS_Application::startNewShower();
        if (res) return res;
        if (current_case != last_case) {
            if (kerma_tot) kerma_tot->setHistory(current_case);
            if (kerma_pri) kerma_pri->setHistory(current_case);
            last_case = current_case;
        }
        return 0;
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

    // ---- Per-bunch scoring (reset after each bunch) ----
    EGS_ScoringArray *kerma_tot;
    EGS_ScoringArray *kerma_pri;

    // ---- E*mu_en/rho interpolator (log E abscissa) ----
    EGS_Interpolator *E_Muen_Rho;

    // ---- Bunch bookkeeping ----
    int       n_bunches;
    long long n_per_bunch;
    int       current_bunch;
    double  **K_bunch;    // [n_bunches][n_scoring]
    double  **Kp_bunch;   // [n_bunches][n_scoring]

    EGS_I64 last_case;

    /*------------------------------------------------------------------------
      scoreFD_all

      Ray-trace from the current photon position along direction u,
      accumulating optical depth Lambda.

      At each scoring shell encountered:
        score += w * exp(-Lambda) * (E*muen/rho) * t_shell / V_shell[k]
        added to kerma_pri (primary) or kerma_tot (scattered).

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
                EGS_Float score = wt * std::exp(-Lambda) * emuen_rho
                                  * t_shell / V_shell[k];
                kerma_tot->score(k, score);        // total = primary + scattered
                if (is_primary) kerma_pri->score(k, score);
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

    // ---- Scoring arrays ----
    kerma_tot = new EGS_ScoringArray(n_scoring);
    kerma_pri = new EGS_ScoringArray(n_scoring);

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

    K_bunch  = new double*[n_bunches];
    Kp_bunch = new double*[n_bunches];
    for (int b = 0; b < n_bunches; b++) {
        K_bunch[b]  = new double[n_scoring]();
        Kp_bunch[b] = new double[n_scoring]();
    }

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
  runSimulation — staged + iterative backscatter-aware loop.

  APP_MAIN calls finishSimulation() which calls outputResults() after we
  return, so we do NOT call outputResults() here.
----------------------------------------------------------------------------*/
int EGS_ShieldApplication::runSimulation() {
    EGS_Timer timer;
    for (current_bunch = 0; current_bunch < n_bunches; current_bunch++) {
        timer.start();
        egsInformation("egs_shield: bunch %d / %d  (%lld photons)",
                       current_bunch + 1, n_bunches, n_per_bunch);

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

        static const int MAX_ITER = 500;
        for (int iter = 0; iter < MAX_ITER; iter++) {
            double w_cur = 0;
            for (int m = 0; m < n_combing; m++)
                for (int j = 0; j < containers[m]->size(); j++)
                    w_cur += (*containers[m])[j].wt;
            if (w_cur <= 0 || w_cur < 1e-10 * w_initial) break;

            for (int m = 0; m < n_combing; m++) {
                std::swap(containers[m], snapshots[m]);
                containers[m]->clean();
                snapshots[m]->comb(rndm, comb_target);
            }
            for (int m = 0; m < n_combing; m++)
                if (snapshots[m]->size() > 0)
                    replayContainer(snapshots[m]);
        }

        endBunch();
        egsInformation("  %.1f s\n", timer.time());
    }
    return 0;
}


/*----------------------------------------------------------------------------
  endBunch
----------------------------------------------------------------------------*/
void EGS_ShieldApplication::endBunch() {
    for (int k = 0; k < n_scoring; k++) {
        double r, dr;
        kerma_tot->currentResult(k, r, dr);
        K_bunch[current_bunch][k] = r;
        kerma_pri->currentResult(k, r, dr);
        Kp_bunch[current_bunch][k] = r;
    }
    kerma_tot->reset();
    kerma_pri->reset();
}


/*----------------------------------------------------------------------------
  outputResults
----------------------------------------------------------------------------*/
void EGS_ShieldApplication::outputResults() {
    static const double MeVtoGy = 1.6021773e-10;

    egsInformation("\n");
    egsInformation("  %-6s  %-14s  %-14s  %-12s  %-10s\n",
                   "Shell", "K_tot/[Gy]", "K_pri/[Gy]", "BUF", "sigma/%");
    egsInformation("  %s\n", string(60, '-').c_str());

    for (int k = 0; k < n_scoring; k++) {
        double sum_buf  = 0.0;
        double sum_buf2 = 0.0;
        double Kt_mean  = 0.0;
        double Kp_mean  = 0.0;
        int    n_valid  = 0;

        for (int b = 0; b < n_bunches; b++) {
            double Kp = Kp_bunch[b][k];
            if (Kp <= 0.0) continue;
            double buf  = K_bunch[b][k] / Kp;
            sum_buf    += buf;
            sum_buf2   += buf * buf;
            Kt_mean    += K_bunch[b][k];
            Kp_mean    += Kp;
            n_valid++;
        }

        if (n_valid < 2) {
            egsInformation("  %-6d  (insufficient bunches with K_pri > 0)\n", k);
            continue;
        }

        double mean = sum_buf / n_valid;
        double var  = (sum_buf2 - n_valid * mean * mean) / (n_valid - 1);
        double sem  = std::sqrt(std::max(var, 0.0) / n_valid);
        Kt_mean /= n_valid;
        Kp_mean /= n_valid;

        egsInformation("  %-6d  %-14.6g  %-14.6g  %-12.6g  %-10.4f\n",
                       k, Kt_mean * MeVtoGy, Kp_mean * MeVtoGy, mean,
                       100.0 * sem / std::max(mean, 1e-30));
    }
    egsInformation("\n");
}


/*----------------------------------------------------------------------------
  Destructor
----------------------------------------------------------------------------*/
EGS_ShieldApplication::~EGS_ShieldApplication() {
    delete kerma_tot;
    delete kerma_pri;
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
    if (K_bunch) {
        for (int b = 0; b < n_bunches; b++) delete [] K_bunch[b];
        delete [] K_bunch;
    }
    if (Kp_bunch) {
        for (int b = 0; b < n_bunches; b++) delete [] Kp_bunch[b];
        delete [] Kp_bunch;
    }
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
