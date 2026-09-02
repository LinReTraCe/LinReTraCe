#!/usr/bin/env python3
"""Tests for the refinement cascade (patch 20).

Covers the kernel-width ladder, the per-stage tolerance policy and the
output-size budget helpers.  Run directly:

    python3 testsuite/tests/test_kernel_cascade.py
"""

from __future__ import annotations

import sys
from pathlib import Path

_root = Path(__file__).resolve().parents[2]
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))

import numpy as np

from structure.meshrefine import kernel_ladder, kernel_width, stage_tolerances
from scripts.kmesh_refinement import format_bytes, parse_size, predict_output_bytes
from structure.units import kB_eV


FAILURES = []


def check(name, condition, detail=""):
    if condition:
        print("PASS  %s" % name)
    else:
        print("FAIL  %s   %s" % (name, detail))
        FAILURES.append(name)


def test_default_ladder_is_a_single_stage():
    """No maxima given -> exactly the pre-cascade behaviour.

    This is the backwards-compatibility guarantee: every existing invocation
    of ltb/lwann refine must keep running one loop at (T_min, gamma_min).
    """
    ladder = kernel_ladder(1.0, 1e-3)
    check("test_default_ladder_is_a_single_stage",
          ladder == [(1.0, 1e-3)], repr(ladder))


def test_ladder_is_widest_first_and_ends_exactly_at_the_target():
    """Ordering is the whole point: the wide stage places the parents.

    Taking max() over the ladder inside one loop would NOT do this, since on a
    coarse mesh the narrow-kernel defect dominates and the maximum is always
    the narrow entry.
    """
    ladder = kernel_ladder(1.0, 1e-3, T_max=300.0, gamma_max=1e-2)
    widths = [kernel_width(t, g) for t, g in ladder]
    check("test_ladder_is_widest_first_and_ends_exactly_at_the_target",
          len(ladder) > 1
          and all(a > b for a, b in zip(widths, widths[1:]))
          and ladder[-1] == (1.0, 1e-3)
          and abs(widths[0] - max(kB_eV * 300.0, 1e-2)) < 1e-12,
          "widths=%s" % widths)


def test_ladder_spacing_is_logarithmic_at_the_requested_ratio():
    """Successive widths differ by a constant factor <= ratio."""
    ladder = kernel_ladder(1.0, 1e-3, T_max=300.0, gamma_max=1e-2, ratio=3.0)
    widths = np.array([kernel_width(t, g) for t, g in ladder])
    steps = widths[:-1] / widths[1:]
    check("test_ladder_spacing_is_logarithmic_at_the_requested_ratio",
          np.all(steps <= 3.0 + 1e-9) and np.ptp(steps) < 1e-9,
          "steps=%s" % steps)


def test_ladder_respects_the_stage_cap():
    ladder = kernel_ladder(1.0, 1e-6, T_max=1000.0, gamma_max=1.0,
                           ratio=2.0, max_stages=4)
    check("test_ladder_respects_the_stage_cap",
          len(ladder) == 4 and ladder[-1] == (1.0, 1e-6), repr(len(ladder)))


def test_maxima_below_the_minima_are_raised_not_honoured():
    """The cascade must always END at the target corner."""
    ladder = kernel_ladder(10.0, 1e-2, T_max=1.0, gamma_max=1e-4)
    check("test_maxima_below_the_minima_are_raised_not_honoured",
          ladder == [(10.0, 1e-2)], repr(ladder))


def test_gamma_only_ladder_leaves_temperature_alone():
    ladder = kernel_ladder(1.0, 1e-3, gamma_max=1e-1)
    check("test_gamma_only_ladder_leaves_temperature_alone",
          len(ladder) > 1 and all(t == 1.0 for t, _ in ladder), repr(ladder))


def test_zero_gamma_min_does_not_break_the_interpolation():
    """gamma_min = 0 is legal: the kernel is then purely thermal."""
    try:
        ladder = kernel_ladder(1.0, 0.0, T_max=300.0, gamma_max=1e-2)
        ok = ladder[-1] == (1.0, 0.0) and all(np.isfinite(g) for _, g in ladder)
    except Exception as exc:                            # pragma: no cover
        ok, ladder = False, exc
    check("test_zero_gamma_min_does_not_break_the_interpolation", ok, repr(ladder))


def test_default_stage_tolerance_is_uniform_across_widths():
    """exponent = 0 holds h/W fixed, i.e. grades the mesh self-similarly."""
    ladder = kernel_ladder(1.0, 1e-3, T_max=300.0, gamma_max=1e-2)
    tols = stage_tolerances(ladder, 5e-3)
    check("test_default_stage_tolerance_is_uniform_across_widths",
          len(tols) == len(ladder) and all(t == 5e-3 for t in tols), repr(tols))


def test_positive_exponent_loosens_only_the_wide_stages():
    ladder = kernel_ladder(1.0, 1e-3, T_max=300.0, gamma_max=1e-2)
    tols = stage_tolerances(ladder, 5e-3, exponent=1.0)
    widths = [kernel_width(t, g) for t, g in ladder]
    expected_first = 5e-3 * widths[0] / widths[-1]
    check("test_positive_exponent_loosens_only_the_wide_stages",
          abs(tols[-1] - 5e-3) < 1e-15
          and abs(tols[0] - expected_first) < 1e-12
          and all(a > b for a, b in zip(tols, tols[1:])),
          repr(tols))


def test_size_strings_parse_with_decimal_and_binary_prefixes():
    cases = [("2048", 2048), ("500MB", 500 * 1000**2), ("4GB", 4 * 1000**3),
             ("2.5GiB", int(2.5 * 1024**3)), ("1e9", 10**9), (" 512 kiB", 512 * 1024)]
    bad = [(s, parse_size(s), want) for s, want in cases if parse_size(s) != want]
    rejected = True
    for s in ("", "-5MB", "10 parsecs", "MB"):
        try:
            parse_size(s)
            rejected = False
        except ValueError:
            pass
    check("test_size_strings_parse_with_decimal_and_binary_prefixes",
          not bad and rejected, "bad=%s rejected=%s" % (bad, rejected))


def test_output_size_prediction_is_linear_in_the_kpoint_count(tmpdir):
    """The generator writes one moments dataset per k-point (KI-02), so the
    file size is linear in the mesh up to a small fixed header."""
    probe = tmpdir / "probe.bin"
    probe.write_bytes(b"x" * 1000)
    got = predict_output_bytes(probe, 100, 350)
    check("test_output_size_prediction_is_linear_in_the_kpoint_count",
          abs(got - 3500.0) < 1e-6 and predict_output_bytes(probe, 0, 350) == 0.0,
          repr(got))


def test_format_bytes_is_readable():
    check("test_format_bytes_is_readable",
          format_bytes(512) == "512 B"
          and format_bytes(1024).startswith("1.0 kiB")
          and format_bytes(5 * 1024**3).startswith("5.0 GiB"),
          format_bytes(1024))


def main():
    import tempfile
    test_default_ladder_is_a_single_stage()
    test_ladder_is_widest_first_and_ends_exactly_at_the_target()
    test_ladder_spacing_is_logarithmic_at_the_requested_ratio()
    test_ladder_respects_the_stage_cap()
    test_maxima_below_the_minima_are_raised_not_honoured()
    test_gamma_only_ladder_leaves_temperature_alone()
    test_zero_gamma_min_does_not_break_the_interpolation()
    test_default_stage_tolerance_is_uniform_across_widths()
    test_positive_exponent_loosens_only_the_wide_stages()
    test_size_strings_parse_with_decimal_and_binary_prefixes()
    with tempfile.TemporaryDirectory() as d:
        test_output_size_prediction_is_linear_in_the_kpoint_count(Path(d))
    test_format_bytes_is_readable()
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
