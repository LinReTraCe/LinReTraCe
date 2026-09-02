#!/usr/bin/env python3
"""Tests for the k-mesh anisotropy estimator and the parent-mesh qualifier.

Uses analytic tight-binding dispersions rather than generator output, so the
tests need no HDF5 files and run in well under a second.  Run directly:

    python3 testsuite/tests/test_mesh_anisotropy.py
"""

from __future__ import annotations

import sys
from pathlib import Path

_root = Path(__file__).resolve().parents[2]
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))

import numpy as np

from structure.meshrefine import (axis_anisotropy, compute_error, build_band_axis,
                                  decimate_mask, estimate_metric_scaling,
                                  metric_components,
                                  infer_mesh_multiple, recommend_density,
                                  suggest_mesh_ratio, _symmetry_axis_orbits)

FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print("PASS  %s" % name)
    else:
        print("FAIL  %s   %s" % (name, detail))
        FAILURES.append(name)


def square_model(nk=32, tx=0.25, ty=0.05, tz=None):
    """Analytic single-band model with per-axis hoppings.

    E(f) = -2 tx cos(2 pi fx) - 2 ty cos(2 pi fy) [- 2 tz cos(2 pi fz)]

    Returns (energies, momentsDiagonal, weights, kvec, dims) in the layout
    that axis_anisotropy expects.  Note the energy variation per FRACTIONAL
    step is set by the hoppings alone and is independent of the lattice
    constants, which is exactly why the geometric rule nk ~ |b| is wrong for
    tight-binding models.
    """
    dims = np.array([True, True, tz is not None])
    axes = [np.arange(nk) / nk, np.arange(nk) / nk,
            np.arange(nk) / nk if tz is not None else np.zeros(1)]
    FX, FY, FZ = np.meshgrid(*axes, indexing='ij')
    E = -2 * tx * np.cos(2 * np.pi * FX) - 2 * ty * np.cos(2 * np.pi * FY)
    if tz is not None:
        E = E - 2 * tz * np.cos(2 * np.pi * FZ)
    # v_a = dE/dk_a = (1/L_a) dE/df_a ; use L_a = 1 so k and f coincide
    vx = 4 * np.pi * tx * np.sin(2 * np.pi * FX)
    vy = 4 * np.pi * ty * np.sin(2 * np.pi * FY)
    vz = (4 * np.pi * tz * np.sin(2 * np.pi * FZ)
          if tz is not None else np.zeros_like(FX))
    n = E.size
    energies = E.reshape(n, 1)
    moments = np.stack([vx.ravel() ** 2, vy.ravel() ** 2, vz.ravel() ** 2],
                       axis=-1).reshape(n, 1, 3)
    weights = np.full(n, 1.0 / n)
    kvec = 2 * np.pi * np.eye(3)
    return energies, moments, weights, kvec, dims


def test_anisotropy_recovers_the_hopping_ratio():
    """t_x/t_y = 5 must show up as an axis ratio of roughly 5, not 1."""
    d = axis_anisotropy(*square_model(tx=0.25, ty=0.05), 300.0, 0.01)
    raw, _ = suggest_mesh_ratio(d, np.array([True, True, False]))
    check("test_anisotropy_recovers_the_hopping_ratio",
          4.0 < raw[0] / raw[1] < 12.0 and abs(raw[1] - 1.0) < 1e-12,
          "raw=%s" % raw)


def test_isotropic_model_returns_unity():
    d = axis_anisotropy(*square_model(tx=0.25, ty=0.25), 300.0, 0.01)
    raw, sug = suggest_mesh_ratio(d, np.array([True, True, False]))
    check("test_isotropic_model_returns_unity",
          np.allclose(raw[:2], 1.0, atol=1e-9) and np.allclose(sug, 1.0),
          "raw=%s sug=%s" % (raw, sug))


def test_estimate_is_independent_of_the_mesh_it_is_measured_on():
    """The estimator must be usable on a cheap trial mesh.

    Measured convergence on this model: nk = 16 gives 6.09, nk = 32 gives 6.97
    and nk = 48 gives 6.91, so it is settled by about nk = 32 but a very coarse
    trial mesh under-estimates the ratio by roughly 12%.  Since the estimate is
    used as a lower bound that gets rounded up, erring low on a coarse mesh is
    the harmless direction.
    """
    ratios = []
    for nk in (32, 48, 64):
        d = axis_anisotropy(*square_model(nk=nk, tx=0.25, ty=0.05), 300.0, 0.01)
        ratios.append(d[0] / d[1])
    check("test_estimate_is_independent_of_the_mesh_it_is_measured_on",
          np.ptp(ratios) / np.mean(ratios) < 0.05, "ratios=%s" % ratios)


def test_inactive_dimension_is_reported_as_zero():
    d = axis_anisotropy(*square_model(tx=0.25, ty=0.05), 300.0, 0.01)
    check("test_inactive_dimension_is_reported_as_zero", d[2] == 0.0, repr(d))


def test_three_dimensional_model_ranks_all_three_axes():
    d = axis_anisotropy(*square_model(nk=16, tx=0.25, ty=0.05, tz=0.125),
                        300.0, 0.01)
    order = np.argsort(d)[::-1]
    check("test_three_dimensional_model_ranks_all_three_axes",
          tuple(order) == (0, 2, 1) and np.all(d > 0), "d=%s" % d)


def test_suggested_ratio_rounds_up_and_keeps_the_base_at_one():
    """A raw 6.97 must suggest 7, not 4.

    Forcing each ratio to be even and then renormalising halved it; parity
    belongs on the division counts, not on the ratio.
    """
    raw, sug = suggest_mesh_ratio(np.array([6.972, 1.0, 0.0]),
                                  np.array([True, True, False]))
    check("test_suggested_ratio_rounds_up_and_keeps_the_base_at_one",
          sug[0] == 7.0 and sug[1] == 1.0 and sug[2] == 1.0,
          "raw=%s sug=%s" % (raw, sug))


def test_symmetry_orbits_merge_only_coupled_axes():
    c4 = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    mm2 = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
    dims = np.array([True, True, False])
    coupled = _symmetry_axis_orbits(np.stack([np.eye(3, dtype=int), c4]), dims)
    apart = _symmetry_axis_orbits(np.stack([np.eye(3, dtype=int), mm2]), dims)
    check("test_symmetry_orbits_merge_only_coupled_axes",
          any(len(o) > 1 and 0 in o and 1 in o for o in coupled)
          and all(not (0 in o and 1 in o) for o in apart),
          "coupled=%s apart=%s" % (coupled, apart))


def test_density_recommendation_grows_only_active_axes():
    rec, ok, _ = recommend_density((32, 32, 1), 0.0178, 5e-3,
                                   np.array([True, True, False]), multiple=2)
    check("test_density_recommendation_grows_only_active_axes",
          ok and rec[2] == 1 and rec[0] == rec[1] and rec[0] > 32
          and rec[0] % 2 == 0, repr(rec))


def test_saturation_floor_is_reported_as_unreachable():
    """Below the metric's floor no uniform mesh can help, and saying 'grow the
    mesh' would be a lie."""
    rec, ok, msg = recommend_density((16, 16, 16), 0.0050, 1e-3,
                                     np.array([True, True, True]),
                                     floor_error=0.0043)
    check("test_saturation_floor_is_reported_as_unreachable",
          not ok and "saturates" in msg and rec == (16, 16, 16), msg)


def test_recommendation_is_a_noop_when_already_converged():
    rec, ok, _ = recommend_density((48, 48, 1), 1e-4, 5e-3,
                                   np.array([True, True, False]))
    check("test_recommendation_is_a_noop_when_already_converged",
          ok and rec == (48, 48, 1), repr(rec))


def test_commensurability_is_inherited_from_the_starting_mesh():
    """Graphene needs n divisible by 3 to put K on the mesh; a parent without
    it stalls the refinement completely."""
    dims2 = np.array([True, True, False])
    got = (infer_mesh_multiple((48, 48, 1), dims2),
           infer_mesh_multiple((16, 16, 16), np.array([True] * 3)),
           infer_mesh_multiple((49, 49, 1), dims2))
    check("test_commensurability_is_inherited_from_the_starting_mesh",
          got[0] == 6 and got[1] == 4 and got[2] == 1, repr(got))


def graphene_grid(n):
    """Regular graphene mesh: fractional points and the two band energies."""
    i = np.arange(n) / n
    KX, KY = np.meshgrid(i, i, indexing='ij')
    pts = np.stack([KX.ravel(), KY.ravel(), np.zeros(n * n)], axis=-1)
    E = np.abs(1 + np.exp(2j * np.pi * KX) + np.exp(2j * np.pi * KY)).ravel()
    return pts, np.stack([-E, E], axis=-1)


def test_decimation_reproduces_the_coarser_mesh_exactly():
    """Decimating a regular mesh must give the coarse mesh's ENERGY SET.

    Symmetry operations are integer matrices in fractional coordinates, so the
    star of a coarse-grid point stays on the coarse grid; this is what makes
    the trick valid on an irreducible wedge too.
    """
    pts, E = graphene_grid(48)
    ok = True
    for m in (2, 3, 4):
        mask = decimate_mask(pts, (48, 48, 1), m)
        got = compute_error(build_band_axis(E[mask]), 300.0, 1e-3)
        _, Ec = graphene_grid(48 // m)
        want = compute_error(build_band_axis(Ec), 300.0, 1e-3)
        ok = ok and abs(got - want) < 1e-9
    check("test_decimation_reproduces_the_coarser_mesh_exactly", ok)


def test_exponent_is_the_active_dimension_count():
    """Assuming 1/n asked for 768x768 where 420x420 reaches the target.

    Measured local slopes are 2.3 in 2D and 3.1 in 3D, so p = d slightly
    under-estimates and therefore errs towards a denser mesh.
    """
    pts, E = graphene_grid(48)
    p2, raw, _fit, _n = estimate_metric_scaling(pts, E, (48, 48, 1), 300.0, 1e-3)
    pts3 = np.stack(np.meshgrid(*(np.arange(8) / 8,) * 3, indexing='ij'), -1).reshape(-1, 3)
    E3 = (-0.5 * np.sum(np.cos(2 * np.pi * pts3), axis=1)).reshape(-1, 1)
    p3, _, _, _ = estimate_metric_scaling(pts3, E3, (8, 8, 8), 300.0, 1e-2)
    # raw is now the REFINABLE component, not the total (floor 0.000123 here)
    check("test_exponent_is_the_active_dimension_count",
          p2 == 2.0 and p3 == 3.0 and abs(raw - 0.302175) < 1e-4,
          "p2=%s p3=%s raw=%.6f" % (p2, p3, raw))


def test_metric_splits_into_a_refinable_part_and_a_fixed_floor():
    """The floor is kernel weight outside the sampled band range and does not
    move with the mesh, so convergence must be judged on the refinable part."""
    floors, refinables = [], []
    for n in (96, 192):
        _pts, E = graphene_grid(n)
        tot, ref, floor = metric_components(build_band_axis(E), 300.0, 1e-3)
        floors.append(floor)
        refinables.append(ref)
        if abs(tot - (ref + floor)) > 1e-12:
            floors.append(float('nan'))
    # The floor settles once the band edges are resolved (from n = 96 here);
    # the refinable part keeps falling.
    check("test_metric_splits_into_a_refinable_part_and_a_fixed_floor",
          len(floors) == 2
          and abs(floors[0] - floors[1]) < 0.05 * max(floors)
          and refinables[1] < 0.6 * refinables[0],
          "floors=%s refinable=%s" % (floors, refinables))


def test_estimate_survives_a_non_monotone_metric():
    """The metric is not monotone in density: graphene at (500 K, 1 meV) scores
    0.0414 at n=48 but 0.0738 at n=66, so an estimate anchored on one lucky
    mesh swung between 132 and 480 for a true answer near 264.  Fitting the
    amplitude at a fixed exponent must keep two very different starting meshes
    in rough agreement.
    """
    recs = []
    for n in (48, 96):
        pts, E = graphene_grid(n)
        p, _raw, fit, _k = estimate_metric_scaling(pts, E, (n, n, 1), 500.0, 1e-3)
        rec, ok, _ = recommend_density((n, n, 1), fit, 5e-3,
                                       np.array([True, True, False]),
                                       multiple=6, exponent=p)
        recs.append(rec[0])
    ratio = max(recs) / min(recs)
    check("test_estimate_survives_a_non_monotone_metric",
          ratio < 1.6 and all(200 <= r <= 420 for r in recs),
          "recs=%s ratio=%.2f" % (recs, ratio))


def test_recommendation_lands_on_the_target():
    """Graphene 48x48 at (300 K, 1 meV) must suggest about 420, not 768."""
    pts, E = graphene_grid(48)
    p, _raw, fit, _n = estimate_metric_scaling(pts, E, (48, 48, 1), 300.0, 1e-3)
    rec, ok, _ = recommend_density((48, 48, 1), fit, 5e-3,
                                   np.array([True, True, False]),
                                   multiple=6, exponent=p)
    check("test_recommendation_lands_on_the_target",
          ok and 350 <= rec[0] <= 550 and rec[2] == 1, repr(rec))


def main():
    test_anisotropy_recovers_the_hopping_ratio()
    test_isotropic_model_returns_unity()
    test_estimate_is_independent_of_the_mesh_it_is_measured_on()
    test_inactive_dimension_is_reported_as_zero()
    test_three_dimensional_model_ranks_all_three_axes()
    test_suggested_ratio_rounds_up_and_keeps_the_base_at_one()
    test_symmetry_orbits_merge_only_coupled_axes()
    test_density_recommendation_grows_only_active_axes()
    test_saturation_floor_is_reported_as_unreachable()
    test_recommendation_is_a_noop_when_already_converged()
    test_commensurability_is_inherited_from_the_starting_mesh()
    test_decimation_reproduces_the_coarser_mesh_exactly()
    test_exponent_is_the_active_dimension_count()
    test_metric_splits_into_a_refinable_part_and_a_fixed_floor()
    test_estimate_survives_a_non_monotone_metric()
    test_recommendation_lands_on_the_target()
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
