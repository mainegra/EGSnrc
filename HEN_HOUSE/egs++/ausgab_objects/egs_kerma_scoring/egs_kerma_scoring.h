/*
###############################################################################
#
#  EGSnrc egs++ kerma scoring ausgab object header
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

  Ausgab objects for collision kerma scoring.

  Three concrete classes — EGS_PlanarKerma, EGS_VolumetricKerma,
  EGS_SphericalKerma — derive from the corresponding fluence AOs and
  override the energyWeight() hook to return E*muen/rho(E) from a
  user-supplied table.

  Input key: emuen file = <path>
  File format (identical to egs_kerma):
    <ndat>
    <E1>  <E1*muen/rho>
    <E2>  <E2*muen/rho>
    ...

  The table is interpolated in log(E) space. Units of E*muen/rho must
  be consistent with the simulation units (MeV, cm^2/g typically).

  Result: kerma [MeV g^-1] per source particle, using the same
  normalization as the parent fluence AO.
*/

#ifndef EGS_KERMA_SCORING_
#define EGS_KERMA_SCORING_

#include "egs_fluence_scoring.h"
#include "egs_interpolator.h"

#ifdef BUILD_KERMA_SCORING_DLL
    #define EGS_KERMA_SCORING_EXPORT EGS_EXPORT
#else
    #define EGS_KERMA_SCORING_EXPORT EGS_IMPORT
#endif

/*! \brief Planar kerma scoring AO.

  Inherits the planar fluence scorer and multiplies every scored weight
  by E*muen/rho(E) looked up from the user-supplied table.
*/
class EGS_KERMA_SCORING_EXPORT EGS_PlanarKerma : public EGS_PlanarFluence {
public:
    EGS_PlanarKerma(const string &Name="", EGS_ObjectFactory *f=nullptr);
    ~EGS_PlanarKerma() { delete E_Muen_Rho; }

    EGS_Float energyWeight(EGS_Float logE) const override {
        return E_Muen_Rho->interpolateFast(logE);
    }
    string    scoringType()   const override { return "kerma"; }
    string    columnHeader()  const override { return "K/[Gy]"; }
    string    quantityUnits() const override { return "Gy"; }
    EGS_Float outputFactor()  const override { return 1.6021773e-10; }

    void loadTable(EGS_Input *inp);

private:
    EGS_Interpolator *E_Muen_Rho;
};

/*! \brief Volumetric kerma scoring AO.

  Inherits the volumetric fluence scorer (TL + FD estimators) and
  multiplies every scored weight by E*muen/rho(E).
*/
class EGS_KERMA_SCORING_EXPORT EGS_VolumetricKerma : public EGS_VolumetricFluence {
public:
    EGS_VolumetricKerma(const string &Name="", EGS_ObjectFactory *f=nullptr);
    ~EGS_VolumetricKerma() { delete E_Muen_Rho; }

    EGS_Float energyWeight(EGS_Float logE) const override {
        return E_Muen_Rho->interpolateFast(logE);
    }
    string    scoringType()   const override { return "kerma"; }
    string    columnHeader()  const override { return "K/[Gy]"; }
    string    quantityUnits() const override { return "Gy"; }
    EGS_Float outputFactor()  const override { return 1.6021773e-10; }

    void loadTable(EGS_Input *inp);

private:
    EGS_Interpolator *E_Muen_Rho;
};

/*! \brief Spherical kerma scoring AO.

  Inherits the spherical fluence scorer (crossing + FD estimators) and
  multiplies every scored weight by E*muen/rho(E).
*/
class EGS_KERMA_SCORING_EXPORT EGS_SphericalKerma : public EGS_SphericalFluence {
public:
    EGS_SphericalKerma(const string &Name="", EGS_ObjectFactory *f=nullptr);
    ~EGS_SphericalKerma() { delete E_Muen_Rho; }

    EGS_Float energyWeight(EGS_Float logE) const override {
        return E_Muen_Rho->interpolateFast(logE);
    }
    string    scoringType()   const override { return "kerma"; }
    string    columnHeader()  const override { return "K/[Gy]"; }
    string    quantityUnits() const override { return "Gy"; }
    EGS_Float outputFactor()  const override { return 1.6021773e-10; }

    void loadTable(EGS_Input *inp);

private:
    EGS_Interpolator *E_Muen_Rho;
};

#endif
