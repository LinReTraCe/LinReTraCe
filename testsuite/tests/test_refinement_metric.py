#!/usr/bin/env python3
"""
testsuite/tests/test_refinement_metric.py

Regression tests for the df/da refinement metric and hotspot selection
(KI-R14).  Two properties are pinned:

  * ``compute_error`` is the SUM OF ABSOLUTE per-panel quadrature defects, not
    the absolute value of their signed sum.  The signed version could dip on
    one iterate through accidental cancellation between the peak and the tails
    and rise again on the next, strictly finer one.

  * ``select_refinement_panels`` marks panels by the defect they carry, and
    counts a panel as in-window when it OVERLAPS the window rather than when
    its centre lies inside.  On a coarse mesh the panel straddling mu can be
    much wider than the window.

Run with:  python3 -m pytest testsuite/tests/test_refinement_metric.py
or simply: python3 testsuite/tests/test_refinement_metric.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if _root not in sys.path:
    sys.path.insert(0, _root)

from structure.meshrefine import (          # noqa: E402
    compute_error,
    df_da,
    hotspot_mask_from_panels,
    panel_defects,
    select_refinement_panels,
)

T_K   = 1.0
GAMMA = 1.0e-3
MU    = 0.0


def _axis(n, emin=-3.0, emax=3.0, cluster=None):
    """Band axis: uniform over the range, optionally clustered around mu."""
    a = np.linspace(emin, emax, n)
    if cluster:
        off = np.geomspace(1e-6, cluster, 200)
        a = np.concatenate((a, MU + off, MU - off, [MU]))
    return np.sort(np.unique(a))


def test_metric_is_the_L1_defect_not_the_signed_sum():
    """The two definitions must differ, and the L1 one must be the larger.

    |sum d_i| <= sum |d_i| always, so the new metric is conservative: nothing
    that was converged under the old one becomes unconverged.
    """
    band = _axis(4000, cluster=1e-2)
    _, defects, exact_total = panel_defects(band, T_K, GAMMA)

    signed = abs(1.0 - float(np.trapezoid(df_da(band, T_K, GAMMA), band)))
    l1     = compute_error(band, T_K, GAMMA)

    assert l1 >= signed - 1e-12, 'L1 metric fell below the signed metric'
    assert np.sum(np.abs(defects)) > abs(np.sum(defects)), (
        'test premise: peak and tail defects must have opposite signs')
    # the reported value is discretisation + truncation of the sampled range
    assert np.isclose(l1, np.sum(np.abs(defects)) + abs(1.0 - exact_total),
                      rtol=1e-12)


def test_metric_decreases_under_pure_refinement():
    """Adding energies must never increase the metric.

    This is the property the old signed metric lacked and the reason the
    graphene run reported a -29% 'regression' on a strictly finer mesh.
    """
    band = _axis(400)
    prev = compute_error(band, T_K, GAMMA)
    for _ in range(6):
        # bisect every panel: the new axis is a strict superset
        band = np.sort(np.unique(np.concatenate((band, 0.5 * (band[:-1] + band[1:])))))
        cur  = compute_error(band, T_K, GAMMA)
        assert cur <= prev * (1.0 + 1e-9), (
            'metric increased under pure refinement: {} -> {}'.format(prev, cur))
        prev = cur


def test_wide_panel_straddling_mu_is_not_discarded():
    """A panel far wider than the window must still be selectable.

    Graphene 24x24x1 has its first energy above the Dirac point near 0.25 eV,
    so the single panel holding essentially the whole defect has its centre at
    -0.125 eV -- outside a 0.1 eV window.  A centre-based test reported
    'converged' at an error of 81.9.
    """
    band = np.array([-3.0, -0.25, 0.0, 0.25, 3.0])
    marked, info = select_refinement_panels(band, T_K, GAMMA, MU,
                                            energy_window=0.1)
    assert marked.any(), 'no panel marked despite a large defect'
    assert info['total'] > 1.0, 'test premise: the defect must be large here'
    # the two panels adjacent to mu are the ones that matter
    mid = 0.5 * (band[:-1] + band[1:])
    assert marked[np.argmin(np.abs(mid - MU))]


def test_selection_captures_the_requested_defect_fraction():
    band = _axis(2000, cluster=1e-2)
    _, defects, _ = panel_defects(band, T_K, GAMMA)
    for frac in (0.5, 0.9, 0.99):
        marked, info = select_refinement_panels(band, T_K, GAMMA, MU,
                                                energy_window=0.1,
                                                defect_fraction=frac)
        assert info['captured'] >= frac * info['total'] - 1e-12
        assert marked.sum() == info['n_marked'] > 0
    # a stricter fraction can never mark fewer panels
    n = [select_refinement_panels(band, T_K, GAMMA, MU, energy_window=0.1,
                                  defect_fraction=f)[1]['n_marked']
         for f in (0.5, 0.9, 0.99)]
    assert n[0] <= n[1] <= n[2]


def test_hotspot_mask_selects_the_bracketing_kpoints():
    """A marked panel must select exactly the k-points owning its endpoints."""
    band = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
    # 4 k-points, 2 bands, energies drawn from the axis
    energies = np.array([[-1.0, 1.0],
                         [-0.5, 0.5],
                         [0.0, 0.5],
                         [-1.0, -0.5]])
    marked = np.array([False, True, False, False])   # panel [-0.5, 0.0]
    mask = hotspot_mask_from_panels(energies, band, marked)
    # k-points carrying -0.5 or 0.0: indices 1, 2, 3
    assert mask.tolist() == [False, True, True, True]

    none = hotspot_mask_from_panels(energies, band, np.zeros(4, dtype=bool))
    assert not none.any()


if __name__ == '__main__':
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
