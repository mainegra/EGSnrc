#!/usr/bin/env python3

###############################################################################
#
#  EGSnrc egs_shield input generator
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

"""
buf_inputs.py -- generate egs_kerma and egs_shield input files for
deep-penetration buildup-factor calculations.

Emits complete, self-consistent .egsinp files for buildup-factor calculations
in a homogeneous sphere with a central photon source, for either code.

Why this exists
---------------
Every (material, energy) point needs a new mean free path, 38 shell radii at
eta*mfp, shell volumes, and either an importance map (egs_kerma) or combing
shells (egs_shield) -- with region indices that shift whenever anything is
inserted.  Doing that by hand produced, in this project alone:

  * importance transitions sitting on detectors, worth 1.2% at 10 mfp
  * combing shells coincident with scoring shells, a 2x overcount
  * a combing shell with no iron buffer, leaving the backscatter loop empty
    and the last shell 10% low

All three are index/placement errors that a generator cannot make.

Geometry conventions (both codes)
---------------------------------
Region 0 is the central core; region i spans radii[i-1] -> radii[i].
Scoring shells are thin air; everything else is the shield material.

egs_kerma layout, per the validated Fe/1 MeV inputs:
    fine section : core + 20 alternating [air, material] pairs -> regions 0..40
    coarse       : groups of THREE regions [material, air, material], so the
                   importance transition at each group boundary falls midway
                   between detectors, never on one.
    importance   : exp(eta) evaluated at the group's inner edge.

egs_shield layout:
    same scoring shells, plus thin combing shells of material placed between
    detector positions (offset, never coincident).

Obtaining the mean free path
----------------------------
mfp = 1 / (mu/rho * rho).  Take mu/rho from NIST XCOM for the material at the
photon energy; rho is the mass density.  Supply it explicitly -- deriving it
from EGSnrc's own data would make the geometry depend on the cross-section set
being tested.  A few values are tabulated in MFP_TABLE below.

Usage
-----
    python3 buf_inputs.py --list
    python3 buf_inputs.py --material iron --energy 1.0 --code both
    python3 buf_inputs.py --material lead --energy 0.5 --mfp 0.42 --code kerma
"""

import argparse
import math
import sys

# ---------------------------------------------------------------------------
# Material data.  mfp in cm, keyed by (material, energy in MeV).
# density correction file names are EGSnrc's, under $HEN_HOUSE/pegs4/density_corrections/
# ---------------------------------------------------------------------------
MATERIALS = {
    "iron":     dict(dcfile="iron",     rho=7.874),
    "lead":     dict(dcfile="lead",     rho=11.350),
    "tungsten": dict(dcfile="tungsten", rho=19.300),
    "water":    dict(dcfile="water_liquid", rho=1.000),
    "concrete": dict(dcfile="concrete_portland", rho=2.300),
}

# mfp = 1/(mu_total/rho * rho), NIST XCOM total attenuation with coherent.
MFP_TABLE = {
    ("iron", 1.0): 2.12007,      # 5.990370e-2 cm2/g * 7.874 g/cm3 -- validated
}

DETECTOR_MEDIUM = "air"
DETECTOR_DC     = "air_dry_nearsealevel"
SHELL_T         = 1.0e-4         # air scoring shell thickness, cm (1 um)
COMB_T          = 1.0e-3         # combing shell thickness, cm


def default_schedule():
    """38 shells: 0.5-10 mfp in 0.5 steps, then 15-100 mfp in 5s."""
    return [0.5 * k for k in range(1, 21)] + [15.0 + 5.0 * j for j in range(18)]


# ---------------------------------------------------------------------------
# Geometry construction
# ---------------------------------------------------------------------------
def shell_volume(eta, mfp, t=SHELL_T):
    r = eta * mfp
    return 4.0 * math.pi * r * r * t


def build_kerma_geometry(etas, mfp, n_fine):
    """
    Returns (radii, scoring_regions, groups) where groups is a list of
    (r_first, r_last, importance) covering every region.

    Fine section: core + alternating pairs.  Coarse: groups of three with the
    detector in the middle, so each importance transition lands midway between
    detectors.
    """
    fine, coarse = etas[:n_fine], etas[n_fine:]
    radii, scoring = [], []

    for e in fine:
        r = e * mfp
        radii += [r - SHELL_T / 2, r + SHELL_T / 2]
        scoring.append(len(radii) - 1)          # air shell region index

    n_fine_reg = len(radii)                     # last fine region index
    groups = []
    prev = fine[-1]
    for e in coarse:
        trans = 0.5 * (prev + e) * mfp          # midway between detectors
        r = e * mfp
        radii += [trans, r - SHELL_T / 2, r + SHELL_T / 2]
        # region i spans radii[i-1]->radii[i]; the three appended radii create
        # regions [material, AIR, material], the detector being the last index.
        scoring.append(len(radii) - 1)          # middle region of the group
        prev = e

    radii.append(1.10 * etas[-1] * mfp)         # outer absorber

    # importance groups: everything up to the first transition is one group
    groups.append((0, n_fine_reg, math.exp(fine[-1])))
    first = n_fine_reg + 1
    prev = fine[-1]
    for j, e in enumerate(coarse):
        inner = 0.5 * (prev + e)
        groups.append((first + 3 * j, first + 3 * j + 2, math.exp(inner)))
        prev = e
    return radii, scoring, groups


def build_shield_geometry(etas, mfp, comb_spacing, comb_start=None):
    """
    Same scoring shells plus thin combing shells of material.

    Combing surfaces are placed at MIDPOINTS between consecutive detectors, so
    they are structurally incapable of coinciding with one, and every combing
    shell has a material buffer on both sides (required: a replayed virtual
    photon must start in the shield, not in a detector).  Every k-th midpoint
    is taken, k chosen to give the requested spacing.

    comb_start defaults to the first midpoint beyond the fine section, matching
    the validated Fe/1 MeV geometry.  The fine section needs no population
    control: forced detection alone is already efficient at shallow depth, and
    combing there costs CPU for no variance reduction.
    """
    mids = [0.5 * (a + b) for a, b in zip(etas, etas[1:])]
    if comb_start is None:
        step = etas[-1] - etas[-2]
        comb_start = next((m for m, (a, b) in zip(mids, zip(etas, etas[1:]))
                           if b - a >= step - 1e-9), mids[0])
    combs, last = [], -1e30
    for m in mids:
        if m >= comb_start - 1e-9 and m >= last + comb_spacing - 1e-9:
            combs.append(m); last = m

    feats = [(x, "s") for x in etas] + [(x, "c") for x in combs]
    feats.sort()
    radii, scoring, combing = [], [], []
    for x, kind in feats:
        r = x * mfp
        if kind == "s":
            radii += [r - SHELL_T / 2, r + SHELL_T / 2]
            scoring.append(len(radii) - 1)
        else:
            radii += [r, r + COMB_T]
            combing.append(len(radii) - 1)
    radii.append(1.10 * etas[-1] * mfp)
    return radii, scoring, combing, combs


# ---------------------------------------------------------------------------
# Invariants -- checked on every generation, not just in tests.
# Each corresponds to a real error made by hand in this project.
# ---------------------------------------------------------------------------
def check_common(radii, scoring, mfp, etas):
    assert all(b > a for a, b in zip(radii, radii[1:])), \
        "radii not strictly increasing"
    assert len(scoring) == len(etas), \
        f"{len(scoring)} detectors for {len(etas)} requested depths"
    for reg, e in zip(scoring, etas):
        mid = 0.5 * (radii[reg - 1] + radii[reg])
        assert abs(mid / mfp - e) < 1e-6, \
            f"region {reg} sits at eta={mid/mfp:.6f}, expected {e}"


def check_kerma(radii, scoring, groups, mfp, etas):
    check_common(radii, scoring, mfp, etas)
    covered = set()
    for a, b, _ in groups:
        assert a <= b, f"malformed importance group {a}-{b}"
        assert not (set(range(a, b + 1)) & covered), f"group {a}-{b} overlaps"
        covered |= set(range(a, b + 1))
    assert covered == set(range(len(radii))), "importance groups do not tile the geometry"
    # a transition ON a detector corrupts that detector's score (1.2% at 10 mfp)
    edges = {a for a, _, _ in groups} | {b for _, b, _ in groups}
    bad = [r for r in scoring if r in edges]
    assert not bad, f"detector(s) {bad} sit on an importance-group boundary"
    imps = [i for _, _, i in groups]
    assert all(b >= a for a, b in zip(imps, imps[1:])), "importance not monotonic"


def check_shield(radii, scoring, combing, mfp, etas):
    check_common(radii, scoring, mfp, etas)
    assert not (set(scoring) & set(combing)), \
        "combing shell coincides with a detector (2x overcount)"
    # a replayed virtual photon must start in the shield, never in a detector
    for c in combing:
        for nb in (c - 1, c + 1):
            assert nb not in scoring, \
                f"combing shell {c} is adjacent to detector {nb}: no material buffer"


# ---------------------------------------------------------------------------
# Emission
# ---------------------------------------------------------------------------
def wrap(items, per, indent, width):
    return " \\\n".join(
        indent + " ".join(f"{x:>{width}}" for x in items[i:i + per])
        for i in range(0, len(items), per))


def media_blocks(mat, radii, scoring, geom_name):
    dc = MATERIALS[mat]["dcfile"]
    setmed = "\n".join(f"            set medium = {r:5d} 1" for r in scoring)
    geom = f"""    :start geometry:
        name    = {geom_name}
        library = egs_spheres
        radii   = \\
{wrap([f"{r:.6f}" for r in radii], 6, "           ", 13)}

        :start media input:
            media = {mat} {DETECTOR_MEDIUM}
            # air scoring shells; all other regions default to {mat}
{setmed}
        :stop media input:
    :stop geometry:"""
    media = f""":start media definition:

    ae = 0.521
    ue = 2.511
    ap = 0.010
    up = 2.0

    :start {mat}:
        density correction file = {dc}
    :stop {mat}:

    :start {DETECTOR_MEDIUM}:
        density correction file = {DETECTOR_DC}
    :stop {DETECTOR_MEDIUM}:

:stop media definition:"""
    return geom, media


TRANSPORT = """:start MC transport parameter:
    Global ECUT                  = 0.521
    Global PCUT                  = 0.010
    Rayleigh scattering          = On
    Photoelectron angular sampling = Off
    Brems angular sampling       = Simple
    Bound Compton scattering     = norej
    Pair angular sampling        = Off
    ESTEPE                       = 0.25
    XIMAX                        = 0.5
    Skin depth for BCA           = 3
    Boundary crossing algorithm  = EXACT
    Electron-step algorithm      = PRESTA-II
:stop MC transport parameter:"""


def source_block(energy, seeds):
    return f""":start source definition:

    :start source:
        name      = point_source
        library   = egs_isotropic_source
        charge    = 0
        :start spectrum:
            type   = monoenergetic
            energy = {energy}
        :stop spectrum:
        :start shape:
            type     = point
            position = 0 0 0
        :stop shape:
    :stop source:

    simulation source = point_source

:stop source definition:

:start rng definition:
    type = xoshiro256++
    initial seeds = {seeds[0]} {seeds[1]}
:stop rng definition:"""


def emit_kerma(mat, energy, mfp, etas, n_fine, ncase, emuen, seeds):
    radii, scoring, groups = build_kerma_geometry(etas, mfp, n_fine)
    check_kerma(radii, scoring, groups, mfp, etas)
    gname = f"{mat}_sphere"
    geom, media = media_blocks(mat, radii, scoring, gname)
    vols = [f"{shell_volume(e, mfp):.4E}" for e in etas]
    imp = " \\\n".join(f"        {a:5d} {b:5d}  {i:.5E}" for a, b, i in groups)
    hdr = f"""##############################################################################
# egs_kerma buildup factors: {energy} MeV photons, {mat} sphere, to {etas[-1]:.0f} mfp
#
# GENERATED by buf_inputs.py -- do not hand-edit region indices.
#
# mfp = {mfp} cm      {len(etas)} scoring shells ({n_fine} fine + {len(etas)-n_fine} coarse)
# geometry: {len(radii)} regions (0 = core, {len(radii)-1} = outer absorber)
#
# Importance map: I = exp(eta) at each group's INNER edge.  The coarse section
# is grouped in threes [material, detector, material] so that every importance
# transition falls midway between detectors and never on one -- a transition
# coincident with a detector corrupts that detector's score (worth 1.2% at
# 10 mfp when it was first found).
#
# Transitions at eta = {", ".join(f"{0.5*(etas[n_fine-1+j]+etas[n_fine+j]):.1f}" for j in range(min(4,len(etas)-n_fine)))} ... mfp
#
# REQUIRES the IS lost-FD-score fix (branch fix-egs_kerma-is-lost-fd-score,
# commit 409e460b).  Without it every result is biased low in proportion to
# the splitting rate -- 10% at 15 mfp for a e^2.5 map.
##############################################################################

:start run control:
    ncase = {ncase}
    nbatch = 1
    rco type = uniform
    interval wait time  = 60000
    number of intervals = 180
:stop run control:

:start geometry definition:

{geom}

    simulation geometry = {gname}

:stop geometry definition:

{media}

{source_block(energy, seeds)}

:start scoring options:

    score primaries = yes
    verbose = yes

    :start calculation geometry:
        geometry name = {gname}
        scoring regions = \\
{wrap([str(r) for r in scoring], 10, "            ", 5)}

        scoring region volumes = \\
{wrap(vols, 5, "            ", 11)}

        importance region ranges = \\
{imp}
    :stop calculation geometry:

    emuen file = {emuen}

:stop scoring options:

:start variance reduction:
:stop variance reduction:

{TRANSPORT}
"""
    return hdr, dict(regions=len(radii), scoring=scoring, groups=groups)


def emit_shield(mat, energy, mfp, etas, nb, npb, emuen, seeds,
                comb_spacing, comb_start, cutoff, comb_target):
    radii, scoring, combing, combs = build_shield_geometry(
        etas, mfp, comb_spacing, comb_start)
    check_shield(radii, scoring, combing, mfp, etas)
    gname = f"{mat}_sphere"
    geom, media = media_blocks(mat, radii, scoring, gname)
    vols = [f"{shell_volume(e, mfp):.4E}" for e in etas]
    hdr = f"""##############################################################################
# egs_shield buildup factors: {energy} MeV photons, {mat} sphere, to {etas[-1]:.0f} mfp
#
# GENERATED by buf_inputs.py -- do not hand-edit region indices.
#
# mfp = {mfp} cm      {len(etas)} scoring shells, {len(combing)} combing shells
# combing spacing {comb_spacing} mfp, placed at midpoints between detectors
# geometry: {len(radii)} regions (0 = core, {len(radii)-1} = outer absorber)
#
# Combing spacing: 10 mfp is at or near the optimum for a HOMOGENEOUS sphere.
# Finer staging was measured (Fe, 1 MeV, 50 bunches each, matched sample count)
# and is worse on every axis -- cost per bunch 281 / 427 / 1643 s at 10 / 4 /
# 2 mfp, and at 25 mfp the figure of merit is 1.5e-2 / 1.3e-4 / 2.9e-4.
#
# The reason: combing is itself a stochastic re-weighting (roulette below the
# mean weight, splitting above), so it pays only when the arriving weight
# spread is large enough to justify the noise it injects.  At 10 mfp spacing
# the spread is ~e^10 and combing earns its keep; at 2 mfp it is ~7x and does
# not -- 49 doses of combing noise instead of 9, for less useful work each.
#
# Do not "improve" this by adding surfaces without re-measuring.  Note this
# conclusion is specific to a HOMOGENEOUS shield: fine staging pays off when
# there are low-attenuation paths for the population to concentrate along.
#
# 'cascade weight cutoff' must let the cascade frontier reach EVERY container.
# The frontier advances ~1 container per iteration while the weight falls by a
# roughly constant factor, so N containers need ~N+6 iterations.  Too high a
# cutoff truncates the cascade and drives the deep shells to BUF = 1 exactly.
##############################################################################

:start run control:
    ncase = {nb * npb}
    nbatch = 1
    rco type = uniform
    interval wait time  = 60000
    number of intervals = 180
:stop run control:

:start geometry definition:

{geom}

    simulation geometry = {gname}

:stop geometry definition:

{media}

{source_block(energy, seeds)}

:start scoring options:

    geometry name = {gname}

    scoring regions = \\
{wrap([str(r) for r in scoring], 10, "        ", 5)}

    scoring volumes = \\
{wrap(vols, 5, "        ", 11)}

    combing regions = \\
{wrap([str(r) for r in combing], 10, "        ", 5)}

    emuen file = {emuen}

    n bunches         = {nb}
    photons per bunch = {npb}

    cascade weight cutoff = {cutoff}
    comb target           = {comb_target}

    bunch statistics   = yes
    cascade diagnostic = yes

:stop scoring options:

:start variance reduction:
:stop variance reduction:

{TRANSPORT}
"""
    return hdr, dict(regions=len(radii), scoring=scoring,
                     combing=combing, comb_etas=combs)


# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--material", default="iron", choices=sorted(MATERIALS))
    p.add_argument("--energy", type=float, default=1.0, help="MeV")
    p.add_argument("--mfp", type=float, default=None,
                   help="cm; required unless the (material,energy) pair is tabulated")
    p.add_argument("--code", default="both", choices=["kerma", "shield", "both"])
    p.add_argument("--max-mfp", type=float, default=100.0)
    p.add_argument("--ncase", default="384e9")
    p.add_argument("--bunches", type=int, default=6400)
    p.add_argument("--per-bunch", type=int, default=1000000)
    p.add_argument("--comb-spacing", type=float, default=10.0)
    p.add_argument("--comb-start", type=float, default=None,
                   help="first combing eta; default = first coarse midpoint")
    p.add_argument("--cutoff", default="1e-50")
    p.add_argument("--comb-target", type=int, default=0)
    p.add_argument("--emuen", default="$EGS_HOME/egs_kerma/emuen_rho_air_1keV-20MeV.data")
    p.add_argument("--seeds", nargs=2, type=int, default=[13579, 24680])
    p.add_argument("--prefix", default=None)
    p.add_argument("--list", action="store_true", help="show tabulated mfp values and exit")
    a = p.parse_args()

    if a.list:
        print("tabulated mfp values (cm):")
        for (m, e), v in sorted(MFP_TABLE.items()):
            print(f"  {m:10s} {e:6.3f} MeV   {v}")
        print("\nmaterials:", ", ".join(sorted(MATERIALS)))
        return 0

    mfp = a.mfp or MFP_TABLE.get((a.material, a.energy))
    if mfp is None:
        sys.exit(f"error: no tabulated mfp for ({a.material}, {a.energy} MeV).\n"
                 f"       supply --mfp; compute as 1/(mu_total/rho * rho) with\n"
                 f"       mu/rho from NIST XCOM (total attenuation with coherent),\n"
                 f"       rho = {MATERIALS[a.material]['rho']} g/cm3.")

    etas = [e for e in default_schedule() if e <= a.max_mfp + 1e-9]
    n_fine = sum(1 for e in etas if e <= 10.0 + 1e-9)
    tag = a.prefix or f"BUF_{a.material}_{a.energy:g}MeV"

    if a.code in ("kerma", "both"):
        txt, info = emit_kerma(a.material, a.energy, mfp, etas, n_fine,
                               a.ncase, a.emuen, a.seeds)
        fn = f"{tag}_kerma.egsinp"
        open(fn, "w").write(txt)
        print(f"{fn}: {info['regions']} regions, {len(info['scoring'])} scoring, "
              f"{len(info['groups'])} importance groups")
    if a.code in ("shield", "both"):
        txt, info = emit_shield(a.material, a.energy, mfp, etas, a.bunches,
                                a.per_bunch, a.emuen, a.seeds, a.comb_spacing,
                                a.comb_start, a.cutoff, a.comb_target)
        fn = f"{tag}_shield.egsinp"
        open(fn, "w").write(txt)
        print(f"{fn}: {info['regions']} regions, {len(info['scoring'])} scoring, "
              f"{len(info['combing'])} combing at eta = "
              f"{', '.join(f'{e:g}' for e in info['comb_etas'][:4])} ...")
    return 0


if __name__ == "__main__":
    sys.exit(main())
