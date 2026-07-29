/*
###############################################################################
#
#  EGSnrc egs++ kerma scoring ausgab object implementation
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
#  Author:          Ernesto Mainegra-Hing, 2026
#
#  Contributors:
#
###############################################################################
*/

#include "egs_kerma_scoring.h"
#include "egs_input.h"
#include "egs_functions.h"

#include <fstream>
#include <vector>
#include <cmath>

using std::vector;
using std::ifstream;
using std::log;

/* --------------------------------------------------------------------------
 * Helper: read E*muen/rho table and build log-energy interpolator.
 * Format (same as egs_kerma 'emuen file'):
 *   <ndat>
 *   <E1>  <val1>
 *   ...
 * -------------------------------------------------------------------------- */
static EGS_Interpolator *buildInterpolator(EGS_Input *inp, const char *who) {
    string fname;
    if (inp->getInput("emuen file", fname)) {
        egsFatal("%s: missing required 'emuen file' input\n", who);
        return nullptr;
    }
    fname = egsExpandPath(fname);
    ifstream ifs(fname.c_str());
    if (!ifs) {
        egsFatal("%s: cannot open emuen file '%s'\n", who, fname.c_str());
        return nullptr;
    }
    int ndat;
    ifs >> ndat;
    if (ndat < 2 || ifs.fail()) {
        egsFatal("%s: emuen file must have at least 2 data points\n", who);
        return nullptr;
    }
    vector<EGS_Float> xdat(ndat), fdat(ndat);
    for (int j = 0; j < ndat; j++) {
        ifs >> xdat[j] >> fdat[j];
    }
    if (ifs.fail()) {
        egsFatal("%s: error reading emuen data file '%s'\n", who, fname.c_str());
        return nullptr;
    }
    egsInformation(
        "\n=============== Kerma Scoring ===============\n"
        "E*muen/rho file : %s\n"
        "Data points     : %d (E range %.4g – %.4g MeV)\n"
        "=============================================\n\n",
        fname.c_str(), ndat, xdat[0], xdat[ndat-1]);
    return new EGS_Interpolator(ndat, log(xdat[0]), log(xdat[ndat-1]),
                                fdat.data());
}

/* -------------------------------------------------------------------------
 * EGS_PlanarKerma
 * ------------------------------------------------------------------------- */

EGS_PlanarKerma::EGS_PlanarKerma(const string &Name, EGS_ObjectFactory *f)
    : EGS_PlanarFluence(Name, f), E_Muen_Rho(nullptr) {}

void EGS_PlanarKerma::loadTable(EGS_Input *inp) {
    E_Muen_Rho = buildInterpolator(inp, "EGS_PlanarKerma");
}

/* -------------------------------------------------------------------------
 * EGS_VolumetricKerma
 * ------------------------------------------------------------------------- */

EGS_VolumetricKerma::EGS_VolumetricKerma(const string &Name, EGS_ObjectFactory *f)
    : EGS_VolumetricFluence(Name, f), E_Muen_Rho(nullptr) {}

void EGS_VolumetricKerma::loadTable(EGS_Input *inp) {
    E_Muen_Rho = buildInterpolator(inp, "EGS_VolumetricKerma");
}

/* -------------------------------------------------------------------------
 * EGS_SphericalKerma
 * ------------------------------------------------------------------------- */

EGS_SphericalKerma::EGS_SphericalKerma(const string &Name, EGS_ObjectFactory *f)
    : EGS_SphericalFluence(Name, f), E_Muen_Rho(nullptr) {}

void EGS_SphericalKerma::loadTable(EGS_Input *inp) {
    E_Muen_Rho = buildInterpolator(inp, "EGS_SphericalKerma");
}

/* -------------------------------------------------------------------------
 * Factory
 * -------------------------------------------------------------------------
 * Input block:
 *   :start ausgab object:
 *       library = egs_kerma_scoring
 *       type    = planar | volumetric | spherical
 *       emuen file = /path/to/emuen_table.dat
 *       # All other keys are passed through to the parent fluence AO.
 *   :stop ausgab object:
 * ------------------------------------------------------------------------- */
extern "C" {

    EGS_KERMA_SCORING_EXPORT EGS_AusgabObject *
    createAusgabObject(EGS_Input *input, EGS_ObjectFactory *f) {
        const static char *func = "createAusgabObject(kerma_scoring)";
        if (!input) {
            egsWarning("%s: null input\n", func);
            return nullptr;
        }

        string type;
        int err = input->getInput("type", type);
        if (err) {
            egsWarning("%s: missing 'type' key (planar/volumetric/spherical)\n", func);
            return nullptr;
        }

        if (input->compare("planar", type)) {
            EGS_PlanarKerma *result = new EGS_PlanarKerma("", f);
            result->setName(input);
            result->initScoring(input);
            result->loadTable(input);
            return result;
        }
        else if (input->compare("volumetric", type)) {
            EGS_VolumetricKerma *result = new EGS_VolumetricKerma("", f);
            result->setName(input);
            result->initScoring(input);
            result->loadTable(input);
            return result;
        }
        else if (input->compare("spherical", type)) {
            EGS_SphericalKerma *result = new EGS_SphericalKerma("", f);
            result->setName(input);
            result->initScoring(input);
            result->loadTable(input);
            return result;
        }
        else {
            egsWarning("%s: unknown type '%s' (expected planar/volumetric/spherical)\n",
                       func, type.c_str());
            return nullptr;
        }
    }

} // extern "C"
