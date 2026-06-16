"""
tests/test_wannier_export.py
============================
Round-trip test for ``structure/wannier_export.py``.

Builds a minimal 1-band simple-cubic tight-binding model, calls
``write_wannier90``, then re-reads the three output files using the *exact*
parsing logic from ``Wannier90Calculation._readWout / _readNnkp / _readHr``
and verifies:

1. All three files exist and are non-empty.
2. ``.wout``:  Length Unit = Ang detected.
3. ``.nnkp``:  real_lattice matches ``tb.rvec``.
4. ``.nnkp``:  recip_lattice × 2π matches ``tb.kvec``.
5. ``.nnkp``:  k-mesh count and point set match ``tb.kpoints``.
6. ``.nnkp``:  nproj == ``tb.energyBandMax``.
7. ``_hr.dat``: nwann, nrp, all mult == 1.
8. ``_hr.dat``: every R-point in ``tb.rpoints`` is present.
9. ``_hr.dat``: H(R) values match ``tb.hr`` element-by-element.
10. Round-trip dispersion: e(k) = Σ_R H(R) exp(i·2π·k·R) reproduces
    ``tb.energies[0]`` at every k-point.

Run from the linretrace package root::

    python tests/test_wannier_export.py
"""

from __future__ import annotations

import math
import os
import sys
import tempfile
import textwrap
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Make sure the package root is on sys.path
# ---------------------------------------------------------------------------
_root = Path(__file__).resolve().parent
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))

from structure.tb             import TightBinding
from structure.wannier_export import write_wannier90


# ---------------------------------------------------------------------------
# Parsers that mirror Wannier90Calculation exactly
# ---------------------------------------------------------------------------

def _parse_wout(fwout: str) -> str:
    """Return 'Ang' or 'Bohr' — mirrors _readWout."""
    with open(fwout) as f:
        for line in f:
            if 'Length Unit' in line:
                return 'Ang' if 'Ang' in line else 'Bohr'
    raise AssertionError("Length Unit not found in " + fwout)


def _parse_nnkp(fnnkp: str, lengthscale: float):
    """Parse all four sections — mirrors _readNnkp field-for-field."""

    # real_lattice
    rvec = []
    with open(fnnkp) as f:
        while not f.readline().startswith('begin real_lattice'):
            pass
        for _ in range(3):
            line = f.readline()
            rvec.append([float(line[:12]), float(line[12:24]), float(line[24:36])])
        assert f.readline().startswith('end real_lattice'), \
            "end real_lattice not found"
    rvec = np.array(rvec) * lengthscale

    # recip_lattice
    kvec = []
    with open(fnnkp) as f:
        while not f.readline().startswith('begin recip_lattice'):
            pass
        for _ in range(3):
            line = f.readline()
            kvec.append([float(line[:12]), float(line[12:24]), float(line[24:36])])
        assert f.readline().startswith('end recip_lattice'), \
            "end recip_lattice not found"
    kvec = np.array(kvec) / lengthscale   # matches _readNnkp: kvec /= lengthscale

    # kpoints
    kpoints = []
    with open(fnnkp) as f:
        while not f.readline().startswith('begin kpoints'):
            pass
        nkp = int(f.readline())
        for _ in range(nkp):
            line = f.readline()
            kpoints.append([float(line[:14]), float(line[14:28]), float(line[28:42])])
        assert f.readline().startswith('end kpoints'), "end kpoints not found"
    kpoints = np.array(kpoints)

    # projections
    plist = []
    with open(fnnkp) as f:
        while not f.readline().startswith('begin projections'):
            pass
        nproj = int(f.readline())
        for _ in range(nproj):
            line = f.readline()
            plist.append([float(line[:10]), float(line[10:21]), float(line[21:32])])
            f.readline()   # skip auxiliary line
        assert f.readline().startswith('end projections'), \
            "end projections not found"
    plist = np.array(plist)

    return rvec, kvec, nkp, kpoints, nproj, plist


def _parse_hr(fhr: str):
    """Parse _hr.dat — mirrors _readHr field-for-field."""
    with open(fhr) as f:
        f.readline()                       # date comment
        nproj = int(f.readline())
        nrp   = int(f.readline())

        # multiplicities
        rmult = []
        while len(rmult) < nrp:
            rmult += f.readline().split()
        rmult = np.array(rmult, dtype=float)

        rlist  = []
        matrix = np.zeros((nrp, nproj, nproj), dtype=complex)
        for ir in range(nrp):
            for _ in range(nproj ** 2):
                line = f.readline()
                rx, ry, rz, p1, p2 = [int(line[i*5:(i+1)*5]) for i in range(5)]
                tr = float(line[25:37])
                ti = float(line[37:])
                matrix[ir, p1 - 1, p2 - 1] = tr + 1j * ti
            rlist.append([rx, ry, rz])     # R from last line of each block
        rlist = np.array(rlist, dtype=int)

        assert f.readline() == '', \
            "_hr.dat has trailing content after last R-block"

    return nproj, nrp, rlist, rmult, matrix


# ---------------------------------------------------------------------------
# Tight-binding model for the test
# ---------------------------------------------------------------------------

# Simple-cubic, 1 band.  On-site energy 0.5 eV, NN hopping t = 0.25 eV.
# Convention in ltb: on-site lines carry E_0 (not -E_0); the sign flip in
# _readTb makes hr[R=0, 0, 0] = +0.5.  All other hr entries are -t = -0.25.
TBDATA = textwrap.dedent("""\
    begin atoms
      1  0.0  0.0  0.0
    end atoms

    begin real_lattice
      5.0  0.0  0.0
      0.0  5.0  0.0
      0.0  0.0  5.0
    end real_lattice

    begin hopping
    # R_x R_y R_z  orb1  orb2  t_real
      0    0   0    1     1     0.5
      1    0   0    1     1     0.25
     -1    0   0    1     1     0.25
      0    1   0    1     1     0.25
      0   -1   0    1     1     0.25
      0    0   1    1     1     0.25
      0    0  -1    1     1     0.25
    end hopping
    """)


# ---------------------------------------------------------------------------
# Test runner
# ---------------------------------------------------------------------------

def run_tests():
    OK = "\033[32m✓\033[0m"   # green tick

    with tempfile.TemporaryDirectory() as tmp:
        tbfile = os.path.join(tmp, 'model.tb')
        with open(tbfile, 'w') as f:
            f.write(TBDATA)

        stem = os.path.join(tmp, 'wann')
        NKX, NKY, NKZ = 4, 4, 4

        # Build TB model on a reducible grid so tb.kpoints is immediately
        # available as the full NKX×NKY×NKZ mesh.
        tb = TightBinding(nkx=NKX, nky=NKY, nkz=NKZ,
                          irreducible=False, kshift=False)
        tb.computeData(tbfile=tbfile, charge=1.0)

        # Export
        write_wannier90(tb, stem)

        fhr   = stem + '_hr.dat'
        fnnkp = stem + '.nnkp'
        fwout = stem + '.wout'

        # ── 1. Files exist and are non-empty ─────────────────────────────
        for fp in [fhr, fnnkp, fwout]:
            assert os.path.isfile(fp),       "Missing: " + fp
            assert os.stat(fp).st_size > 0,  "Empty: "   + fp
        print(OK, "All three files exist and are non-empty.")

        # ── 2. .wout: Length Unit ─────────────────────────────────────────
        unit = _parse_wout(fwout)
        assert unit == 'Ang', f"Expected Ang, got {unit}"
        lengthscale = 1.0
        print(OK, ".wout: Length Unit = Ang")

        # ── 3. .nnkp: real_lattice ────────────────────────────────────────
        rvec, kvec, nkp, kpoints, nproj, plist = _parse_nnkp(fnnkp, lengthscale)
        assert np.allclose(rvec, tb.rvec, atol=1e-5), \
            f"real_lattice mismatch:\n{rvec}\nvs tb.rvec:\n{tb.rvec}"
        print(OK, ".nnkp: real_lattice matches tb.rvec")

        # ── 4. .nnkp: recip_lattice ───────────────────────────────────────
        # _readNnkp does kvec /= lengthscale (=1), so kvec_in_memory = file value.
        # File stores b/(2π); tb.kvec includes 2π.
        # Check: kvec_from_file × 2π ≈ tb.kvec
        assert np.allclose(kvec * 2 * math.pi, tb.kvec, atol=1e-5), \
            f"recip_lattice mismatch (×2π):\n{kvec*2*math.pi}\nvs tb.kvec:\n{tb.kvec}"
        print(OK, ".nnkp: recip_lattice × 2π matches tb.kvec")

        # ── 5. .nnkp: k-mesh ─────────────────────────────────────────────
        assert nkp == NKX * NKY * NKZ, \
            f"nkp {nkp} != {NKX*NKY*NKZ}"
        assert kpoints.shape == (nkp, 3)
        # Every k-point in tb.kpoints must appear in kpoints
        for k in tb.kpoints:
            diffs = np.linalg.norm(kpoints - k[None, :], axis=1)
            assert np.any(diffs < 1e-7), f"k-point {k} not found in written mesh"
        print(OK, f".nnkp: k-mesh has {nkp} points; all tb.kpoints are present")

        # ── 6. .nnkp: nproj ──────────────────────────────────────────────
        assert nproj == tb.energyBandMax, \
            f"nproj {nproj} != energyBandMax {tb.energyBandMax}"
        print(OK, f".nnkp: nproj = {nproj}")

        # ── 7. _hr.dat: header ────────────────────────────────────────────
        nproj_hr, nrp_hr, rlist_hr, rmult_hr, matrix_hr = _parse_hr(fhr)
        assert nproj_hr == tb.energyBandMax, \
            f"nwann in _hr.dat {nproj_hr} != {tb.energyBandMax}"
        assert nrp_hr == tb.nrp, \
            f"nrp in _hr.dat {nrp_hr} != tb.nrp {tb.nrp}"
        assert np.all(rmult_hr == 1.0), \
            "Not all multiplicities are 1"
        print(OK, f"_hr.dat: nwann={nproj_hr}, nrp={nrp_hr}, all mult=1")

        # ── 8. _hr.dat: R-points ─────────────────────────────────────────
        for r in tb.rpoints:
            diffs = np.max(np.abs(rlist_hr - r[None, :]), axis=1)
            assert np.any(diffs == 0), f"R-point {r} missing from _hr.dat"
        print(OK, "_hr.dat: all tb.rpoints are present")

        # ── 9. _hr.dat: H(R) values ──────────────────────────────────────
        for ir_tb, r in enumerate(tb.rpoints):
            diffs = np.max(np.abs(rlist_hr - r[None, :]), axis=1)
            ir_file = int(np.argmin(diffs))
            assert np.allclose(matrix_hr[ir_file], tb.hr[ir_tb], atol=1e-8), \
                "H(R) mismatch at R={}:\nfile:\n{}\ntb:\n{}".format(
                    r, matrix_hr[ir_file], tb.hr[ir_tb])
        print(OK, "_hr.dat: H(R) values match tb.hr exactly")

        # ── 10. Round-trip dispersion ─────────────────────────────────────
        # e(k) = Σ_R H(R) exp(i·2π·k·R) / mult(R)
        # (mult=1 everywhere, so the division is a no-op)
        ek_rt = np.zeros((nkp, nproj_hr, nproj_hr), dtype=complex)
        for ir_file in range(nrp_hr):
            R      = rlist_hr[ir_file].astype(float)
            phases = np.exp(2j * math.pi * kpoints.dot(R))   # (nkp,)
            ek_rt += phases[:, None, None] * matrix_hr[ir_file][None, :, :]

        # For 1-band the Hamiltonian is already diagonal; eigenvalue = H[0,0]
        ek_rt_diag = ek_rt[:, 0, 0].real   # (nkp,)

        for ik_tb, k in enumerate(tb.kpoints):
            # Find the matching k in the written mesh
            diffs    = np.linalg.norm(kpoints - k[None, :], axis=1)
            ik_file  = int(np.argmin(diffs))
            expected = tb.energies[0][ik_tb, 0]
            got      = ek_rt_diag[ik_file]
            assert abs(expected - got) < 1e-5, \
                "Dispersion mismatch at k={}: ltb={:.8f}, from_hr={:.8f}".format(
                    k, expected, got)
        print(OK, "Round-trip dispersion matches ltb energies at all k-points")

    print("\nAll tests passed.\n")


if __name__ == '__main__':
    run_tests()
