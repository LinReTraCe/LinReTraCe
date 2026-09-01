#!/usr/bin/env python3
"""
testsuite/tests/test_refined_mesh_flags.py

Regression tests for KI-R11: `.kmesh/irreducible` must describe the k-MESH
(regular grid, weights reconstructible from multiplicity/(nkx*nky*nkz)) and
NOT the symmetry treatment of the optical elements.  Getting the two confused
made `linretrace` discard `/.kmesh/weights` on every refined irreducible mesh
and assign `weightsum` to every k-point.

Run with:  python3 -m pytest testsuite/tests/test_refined_mesh_flags.py
or simply: python3 testsuite/tests/test_refined_mesh_flags.py
"""

from __future__ import annotations

import os
import sys
import tempfile

import numpy as np

_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if _root not in sys.path:
    sys.path.insert(0, _root)

import h5py

from structure.inout import h5output
from structure.tb import TightBinding
from structure.generators.ltb_gen import LtbGenerator
from scripts.kmesh_refinement import read_source_dims, read_source_symmetry


GRAPHENE = """
begin hopping
   0  0  0     1 2 1.0
   0  0  0     2 1 1.0
  +1  0  0     2 1 1.0
   0 +1  0     2 1 1.0
  -1  0  0     1 2 1.0
   0 -1  0     1 2 1.0
end hopping

begin atoms
  1 0.33333333 0.33333333 0
  1 0.66666667 0.66666667 0
end atoms

begin orbitals
  1 0.33333333 0.33333333 0
  2 0.66666667 0.66666667 0
end orbitals

begin real_lattice
  1   0            0
  0.5 0.8660254038 0
  0   0            1
end real_lattice
"""


def _coarse(tmpdir, nk, irreducible):
    """Write a coarse graphene HDF5 and return its path."""
    tbfile = os.path.join(tmpdir, 'graphene.tbdata')
    if not os.path.exists(tbfile):
        with open(tbfile, 'w') as fh:
            fh.write(GRAPHENE)
    out = os.path.join(tmpdir, 'coarse_{}_{}.hdf5'.format(nk, int(irreducible)))
    tb = TightBinding(nkx=nk, nky=nk, nkz=1, irreducible=irreducible)
    tb.computeData(tbfile=tbfile, charge=2.0)
    h5output(out, tb, tb)
    return tbfile, out


def _refine_once(tbfile, source, out, factor=3):
    """Subdivide every k-point of *source* by *factor* along the active axes.

    Deliberately simpler than the adaptive loop: the point of these tests is
    the metadata written for a custom mesh, not the hotspot selection.
    """
    with h5py.File(source, 'r') as f:
        pts = f['/.kmesh/points'][:]
        wts = f['/.kmesh/weights'][:]
        nkx = int(f['/.kmesh/nkx'][()])
        nky = int(f['/.kmesh/nky'][()])
        dims = np.asarray(f['/.unitcell/dims'][()], dtype=bool)

    deltas = np.array([1.0 / nkx, 1.0 / nky, 0.0])
    offs = (np.arange(factor) - (factor - 1) / 2.0) / factor
    axes = [pts[:, a, None] + offs[None, :] * deltas[a] if dims[a]
            else pts[:, a, None] for a in range(3)]
    n = [ax.shape[1] for ax in axes]
    nchild = n[0] * n[1] * n[2]

    child = np.empty((pts.shape[0], nchild, 3))
    grid = np.meshgrid(np.arange(n[0]), np.arange(n[1]), np.arange(n[2]),
                       indexing='ij')
    idx = [g.ravel() for g in grid]
    for a in range(3):
        child[:, :, a] = axes[a][:, idx[a]]

    new_pts = (child.reshape(-1, 3)) % 1.0
    new_wts = np.repeat(wts / nchild, nchild)

    source_dims = read_source_dims(source)
    _, symop = read_source_symmetry(source)
    gen = LtbGenerator(tb_file=tbfile, filling=2.0, mu=0.0,
                       dims=source_dims, symop=symop)
    gen.generate(new_pts, new_wts, out)
    return out, symop


# ---------------------------------------------------------------------------


def test_refined_irreducible_is_not_flagged_as_regular_grid():
    """KI-R11: a refined wedge must NOT carry .kmesh/irreducible = True.

    linretrace reads that flag to mean "weights are reconstructible from
    multiplicity/(nkx*nky*nkz)" and skips /.kmesh/weights entirely.  With the
    custom-mesh convention (nkx=nky=nkz=1, multiplicity=1) that gives every
    k-point the full weightsum.
    """
    with tempfile.TemporaryDirectory() as tmp:
        tbfile, coarse = _coarse(tmp, 6, irreducible=True)
        refined, symop = _refine_once(tbfile, coarse,
                                      os.path.join(tmp, 'refined.hdf5'))
        assert symop is not None, 'coarse irreducible mesh must carry symop'

        with h5py.File(refined, 'r') as f:
            assert not bool(f['/.kmesh/irreducible'][()]), (
                'refined custom mesh flagged as a regular irreducible grid')
            assert bool(f['/.kmesh/symmetrized'][()]), (
                'symmetrisation marker missing; a continuation run would drop it')
            # the condition linretrace now checks
            mult = f['/.kmesh/multiplicity'][:]
            nkred = (int(f['/.kmesh/nkx'][()]) * int(f['/.kmesh/nky'][()])
                     * int(f['/.kmesh/nkz'][()]))
            assert int(np.sum(mult)) != nkred, 'test premise: custom mesh'
            assert np.isclose(np.sum(f['/.kmesh/weights'][:]),
                              float(f['/.kmesh/weightsum'][()]))


def test_refined_reducible_carries_no_symmetrisation_marker():
    with tempfile.TemporaryDirectory() as tmp:
        tbfile, coarse = _coarse(tmp, 6, irreducible=False)
        refined, symop = _refine_once(tbfile, coarse,
                                      os.path.join(tmp, 'refined.hdf5'))
        assert symop is None
        with h5py.File(refined, 'r') as f:
            assert not bool(f['/.kmesh/irreducible'][()])
            assert not bool(f['/.kmesh/symmetrized'][()])


def _moment_sum(path):
    with h5py.File(path, 'r') as f:
        w = f['/.kmesh/weights'][:]
        m = f['/momentsDiagonal'][:]
    return np.einsum('k,kbd->d', w, m)


def test_coarse_wedge_matches_coarse_reducible_exactly():
    """On a REGULAR grid the wedge sum is exact, so demand machine precision.

    This is the anchor for the looser refined-mesh test below: it pins the
    star-average machinery itself, with no subdivision in the way.
    """
    with tempfile.TemporaryDirectory() as tmp:
        _, coarse_i = _coarse(tmp, 8, irreducible=True)
        _, coarse_r = _coarse(tmp, 8, irreducible=False)
        si, sr = _moment_sum(coarse_i), _moment_sum(coarse_r)
        assert np.allclose(si[:3], sr[:3], rtol=1e-12), (
            'wedge {} != reducible {}'.format(si[:3], sr[:3]))


def test_symmetrised_wedge_preserves_xx_yy_equality():
    """KI-R3: graphene's hexagonal symmetry mixes x and y, so a raw
    (unsymmetrised) wedge sum breaks the equality of the xx and yy Onsager
    elements while the band energies stay correct.  That equality is the
    property the symmetrising path exists to deliver, and unlike the value
    itself it is exact on a refined wedge too.

    The grid division is deliberately NOT a multiple of 3, so that the Dirac
    points K/K' stay off the mesh: at an exact degeneracy |v_nn|^2 is gauge
    dependent and the comparison would be meaningless (KI-11).
    """
    with tempfile.TemporaryDirectory() as tmp:
        tbfile, coarse_i = _coarse(tmp, 8, irreducible=True)
        ref_i, _ = _refine_once(tbfile, coarse_i, os.path.join(tmp, 'ri.hdf5'))
        si = _moment_sum(ref_i)
        assert np.isclose(si[0], si[1], rtol=1e-10), (
            'symmetrised wedge broke the xx/yy equality: {}'.format(si[:2]))


def test_refined_wedge_and_refined_reducible_agree_to_quadrature_error():
    """The two refinement routes are two *different* quadratures of the same
    integral, and must agree only to the discretisation error -- NOT bitwise.

    The exactness argument in ``_setCustomSymmetries`` (each image cell
    S.C_p is tiled by the images of the children of k_p) needs the 3x3
    sub-tiling to be invariant under the point group.  That holds for cubic /
    orthorhombic lattices -- which is what KI-R3 was validated on -- but not
    for a hexagonal one, where the operations mix the fractional axes and map
    the parallelogram patch of children onto a sheared one.  Both quadratures
    carry correct weights, exact group-averaged moments and conserved total
    weight, and both converge to the same limit; they simply sample different
    points.

    Measured for graphene, 8x8x1 subdivided by 3:

        refined reducible  = 1.167906   (identical to the plain 24x24 grid)
        refined wedge      = 1.169197
        converged (192^2)  = 1.172925

    so the two differ by ~0.11%, an order of magnitude below their common
    distance from the converged value.  The tolerance below is set to catch a
    weighting error (which would show up as a factor, not a per-mille shift)
    while tolerating the quadrature difference.
    """
    with tempfile.TemporaryDirectory() as tmp:
        tbfile, coarse_i = _coarse(tmp, 8, irreducible=True)
        _,      coarse_r = _coarse(tmp, 8, irreducible=False)
        ref_i, _ = _refine_once(tbfile, coarse_i, os.path.join(tmp, 'ri.hdf5'))
        ref_r, _ = _refine_once(tbfile, coarse_r, os.path.join(tmp, 'rr.hdf5'))

        si, sr = _moment_sum(ref_i), _moment_sum(ref_r)
        assert np.allclose(si[:3], sr[:3], rtol=1e-2, atol=1e-12), (
            'wedge {} and reducible {} differ by more than the quadrature '
            'error; this is a weighting bug, not a sampling difference'
            .format(si[:3], sr[:3]))

        # total weight is the quantity that must be exact in both routes
        for path in (ref_i, ref_r):
            with h5py.File(path, 'r') as f:
                assert np.isclose(np.sum(f['/.kmesh/weights'][:]),
                                  float(f['/.kmesh/weightsum'][()]), rtol=1e-12)


def test_h5output_rejects_inconsistent_irreducible_flag():
    """The guard that makes KI-R11 impossible to reintroduce silently."""
    with tempfile.TemporaryDirectory() as tmp:
        tbfile, coarse = _coarse(tmp, 6, irreducible=True)
        tb = TightBinding(nkx=6, nky=6, nkz=1, irreducible=True)
        tb.computeData(tbfile=tbfile, charge=2.0)
        # forge the custom-mesh convention onto an irreducible-flagged object
        tb.nkx = tb.nky = tb.nkz = 1
        try:
            h5output(os.path.join(tmp, 'bad.hdf5'), tb, tb)
        except IOError as exc:
            assert 'sum(multiplicity)' in str(exc)
        else:
            raise AssertionError('h5output accepted an inconsistent k-mesh')


if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.WARNING)
    failures = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith('test_') or not callable(fn):
            continue
        try:
            fn()
            print('PASS  {}'.format(name))
        except Exception as exc:                       # noqa: BLE001
            failures += 1
            print('FAIL  {}: {}'.format(name, exc))
    sys.exit(1 if failures else 0)
