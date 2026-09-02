"""
scripts/kmesh_refinement.py

Generator-agnostic iterative k-mesh refinement loop.

This module contains no knowledge of tight-binding, Wannier, or any other
physical model.  It drives the refinement purely through the MeshGenerator
protocol defined in structure/meshgenerator.py.

Entry points for specific generators (ltb, lwann, ...) live in the
corresponding sub-packages under scripts/ and call run_refinement() after
constructing the appropriate generator object.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import h5py
import numpy as np

from structure.meshgenerator import MeshGenerator
from structure.meshrefine import (
    BandData,
    MissingDependencyError,
    build_band_axis,
    axis_anisotropy,
    estimate_metric_scaling,
    infer_mesh_multiple,
    metric_components,
    kernel_ladder,
    kernel_width,
    load_band_data,
    recommend_density,
    refine_kmesh,
    stage_tolerances,
    suggest_mesh_ratio,
    write_custom_mesh,
    select_refinement_panels,
    hotspot_mask_from_panels,
)
from structure.units import kB_eV

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Parameters dataclass
# ---------------------------------------------------------------------------

@dataclass
class RefinementParams:
    """
    All tuning knobs for the refinement loop, collected in one place so that
    generator-specific entry points can build this object from their own
    argparse namespace without duplicating field names.

    Fields
    ------
    initial_hdf5 : Path
        HDF5 file produced by the first (coarse) generator run.
    chemical_potential : float
        mu used for hotspot detection and passed to the generator.
    gamma_min : float
        Minimum scattering rate (eV) for the df/da error metric.
    T_min : float
        Minimum temperature (K) for the df/da error metric.
    error_tol : float
        Stop when the df/da integration error falls below this value.
    max_iter : int
        Hard cap on refinement iterations.
    refinement_factor : int
        Number of subdivisions per k-axis for hotspot regions; should be odd.
    energy_window : float
        CEILING (eV) on the distance from mu at which a band-axis panel may
        still be marked for refinement.  It is not the refinement region:
        panels are selected by the quadrature defect they carry, and this
        only excludes far-away van Hove defects at the band edges, which
        carry real defect but no transport weight.
    workdir : Path
        Directory where intermediate mesh and HDF5 files are written.
    keep_intermediate : bool
        If True, do not delete intermediate files after each iteration.
    plateau_tol : float
        Error-plateau detection (fraction, default 0.05 = 5%): stop the
        refinement if an iteration reduced the error by less than this
        fraction relative to the previous iteration.  The check only becomes
        active from the fourth refinement step onward (iteration index >= 3),
        since the error may equilibrate non-monotonically in the early
        stages.  The persisting plateau value shrinks with a finer initial
        mesh, so on plateau the loop stops with a warning suggesting one.
        Set to 0 (or negative) to disable the check.
    defect_fraction : float
        Fraction (default 0.9) of the total in-window panel defect that the
        marked panels must account for.  Lower values refine fewer k-points
        per iteration and converge more slowly; higher values approach
        "refine everything that is imperfect" and grow the mesh faster.
    T_max, gamma_max : float or None
        Wide corner of the CASCADE.  ``None`` (the default) means "same as
        the minimum", which gives a one-stage ladder and reproduces the
        pre-cascade behaviour exactly.  When either is larger than its
        minimum, the refinement loop is run once per stage of
        :func:`structure.meshrefine.kernel_ladder`, WIDEST KERNEL FIRST, so
        that the wide stage places k-points across the whole +-W_max shell
        while it is still cheap and those points become the parents of the
        narrow stages.  T_min/gamma_min alone say how FINE the mesh must be
        near mu and nothing about how FAR from mu it must be sampled; a mesh
        refined for a 1 meV kernel is otherwise certified converged while
        being useless at 26 meV.
    ladder_ratio : float
        Ratio between the kernel widths of successive cascade stages.
        Default 3.0, matching the default refinement_factor.
    ladder_max_stages : int
        Cap on the number of cascade stages.
    stage_tol_exponent : float
        ``tol_j = error_tol * (W_j/W_min)**exponent``.  0 (default) applies
        the same tolerance at every width, which yields a self-similar mesh
        whose local energy spacing is proportional to its distance from mu.
        See :func:`structure.meshrefine.stage_tolerances`.
    parent_tol : float or None
        Target for the sum-rule metric that the STARTING mesh must already
        meet at the widest kernel, before any refinement.  ``None`` means
        "use error_tol".  Refinement cannot repair a parent that is too
        coarse at the wide corner, so this check refuses rather than warns.
        Kept separate from ``error_tol`` because the two answer different
        questions and because the metric saturates at a system-dependent
        floor (about 0.0043 for the cubic model, 0.002 for graphene).
    skip_parent_check : bool
        Bypass the parent qualification entirely.
    max_output_bytes : int or None
        Refuse to generate a mesh whose HDF5 file is predicted to exceed this
        size, estimated from the bytes-per-k-point of the current file.  The
        run then stops and keeps the last mesh that fitted.  None disables
        the check.
    """
    initial_hdf5:        Path
    chemical_potential:  float
    gamma_min:           float
    T_min:               float
    error_tol:           float = 5e-3
    max_iter:            int   = 10
    refinement_factor:   int   = 3
    energy_window:       float = 0.1
    workdir:             Path  = Path('.')
    keep_intermediate:   bool  = False
    plateau_tol:         float = 0.05
    defect_fraction:     float = 0.9
    T_max:               Optional[float] = None
    gamma_max:           Optional[float] = None
    ladder_ratio:        float = 3.0
    ladder_max_stages:   int   = 6
    stage_tol_exponent:  float = 0.0
    max_output_bytes:    Optional[int] = None
    parent_tol:          Optional[float] = None
    skip_parent_check:   bool = False


# first iteration index (0-based) at which the plateau check is active,
# i.e. the check can first trigger on the Nth refinement step
PLATEAU_MIN_ITER = 4


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_ITER_FILENAME_RE = re.compile(r"refined_iter_(\d+)\.hdf5$")


def _patch_refinement_metadata(
    output_path: Path,
    mesh_path: Path,
    iteration_index: int,
    parent: Path,
) -> None:
    """Copy /.kmesh/cell_deltas from *mesh_path* into *output_path* and stamp
    refinement provenance attributes.

    Background
    ----------
    ``write_custom_mesh`` writes ``cell_deltas`` to a lightweight mesh-only
    HDF5 (``custom_mesh_iter_N.hdf5``).  The generator then reads that file
    and produces a full output HDF5 (``refined_iter_N.hdf5``) via
    ``h5output``, which knows nothing about ``cell_deltas`` and therefore
    does not write it.

    The next iteration's ``load_band_data`` reads the generator output, finds
    no ``cell_deltas`` dataset, and would fall back to ``nkx/nky/nkz`` --
    which are always 1 for custom-mesh runs.  That would give zero cell
    widths, silently turning all subsequent subdivisions into no-ops
    (load_band_data now raises a hard error on that combination instead).

    The fix is to patch ``cell_deltas`` into the generator output immediately
    after it is written, before ``mesh_path`` is deleted.  Opening in append
    mode (``"a"``) leaves all other datasets untouched.

    Provenance attributes
    ---------------------
    ``refinement_iteration`` (int)
        Global 1-based index of the refinement iteration that produced this
        file.  Read back by :func:`_detect_refinement_state` so that a later
        ``ltb refine`` run starting from this file continues the numbering
        instead of restarting at 1 (and thereby colliding with -- or even
        overwriting -- earlier outputs, including the input file itself).
    ``refinement_parent`` (str)
        Path of the HDF5 file this iteration refined.  Informational.
    """
    with h5py.File(mesh_path,   "r") as src, \
         h5py.File(output_path, "a") as dst:
        cell_deltas = src["/.kmesh/cell_deltas"][:]
        kmesh = dst["/.kmesh"]
        if "cell_deltas" in kmesh:
            del kmesh["cell_deltas"]
        kmesh.create_dataset("cell_deltas", data=cell_deltas)
        dst.attrs["refinement_iteration"] = int(iteration_index)
        dst.attrs["refinement_parent"]    = str(parent)


def _detect_refinement_state(path: Path):
    """Classify *path* as fresh coarse mesh or (partially) refined mesh.

    Returns
    -------
    (is_continuation, last_iteration) : tuple[bool, Optional[int]]
        is_continuation : True if the file is an already-refined mesh, i.e.
            refinement should *continue* from it rather than start fresh.
        last_iteration : the 1-based index of the iteration that produced
            the file, if known (from the ``refinement_iteration`` attribute,
            or -- for files produced before that attribute existed -- parsed
            from a ``refined_iter_N.hdf5`` filename); None if unknown.

    Detection is hierarchical:
      1. ``refinement_iteration`` HDF5 attribute (robust against renames).
      2. Presence of ``/.kmesh/cell_deltas`` -> continuation; index possibly
         recovered from the filename.
      3. Otherwise: fresh coarse mesh.
    """
    with h5py.File(path, "r") as f:
        attr = f.attrs.get("refinement_iteration", None)
        has_deltas = "cell_deltas" in f["/.kmesh"]

    if attr is not None:
        return True, int(attr)

    if has_deltas:
        m = _ITER_FILENAME_RE.search(path.name)
        return True, (int(m.group(1)) if m else None)

    return False, None


def _max_existing_index(workdir: Path) -> int:
    """Highest N among refined_iter_N.hdf5 / custom_mesh_iter_N.hdf5 files
    already present in *workdir*; 0 if none.  Used so that a new run in a
    working directory containing outputs of earlier runs numbers its own
    outputs strictly after them, instead of interleaving into gaps.
    """
    pattern = re.compile(r"(?:refined|custom_mesh)_iter_(\d+)\.hdf5$")
    indices = [
        int(m.group(1))
        for f in workdir.iterdir()
        if (m := pattern.match(f.name))
    ]
    return max(indices, default=0)


def _reserve_iteration_paths(
    workdir: Path,
    index: int,
    forbidden: "set[Path]",
):
    """Return (mesh_path, output_path, index) for the next iteration such
    that neither file exists and neither resolves into *forbidden* (the
    input files of this run).  Existing files -- e.g. outputs of a previous
    refinement run in the same working directory -- are never overwritten;
    the index is bumped past them with a warning instead.
    """
    bumped_from = index
    while True:
        mesh_path   = workdir / f"custom_mesh_iter_{index}.hdf5"
        output_path = workdir / f"refined_iter_{index}.hdf5"
        collision = (
            mesh_path.exists()
            or output_path.exists()
            or mesh_path.resolve()   in forbidden
            or output_path.resolve() in forbidden
        )
        if not collision:
            break
        index += 1
    if index != bumped_from:
        logger.warning(
            "Output index bumped %d -> %d: refined_iter_%d.hdf5 / "
            "custom_mesh_iter_%d.hdf5 already exist in %s (or coincide with "
            "the input file) and will not be overwritten.",
            bumped_from, index, bumped_from, bumped_from, workdir,
        )
    return mesh_path, output_path, index


# ---------------------------------------------------------------------------
# Shared CLI surface for the cascade (used by every generator entry point)
# ---------------------------------------------------------------------------

_SIZE_UNITS = {
    "":    1,
    "B":   1,
    "K":   1000,     "KB":  1000,     "KIB": 1024,
    "M":   1000**2,  "MB":  1000**2,  "MIB": 1024**2,
    "G":   1000**3,  "GB":  1000**3,  "GIB": 1024**3,
    "T":   1000**4,  "TB":  1000**4,  "TIB": 1024**4,
}


def parse_size(text: str) -> int:
    """Parse '500MB', '4 GiB', '1e9' or '2048' into a byte count.

    Decimal prefixes are powers of 1000 and binary ones (kiB, MiB, ...) powers
    of 1024, following IEC, so that a budget written as the number read off a
    filesystem quota means what the user expects.
    """
    s = str(text).strip().replace(" ", "")
    m = re.fullmatch(r"([0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)([A-Za-z]*)", s)
    if not m:
        raise ValueError("cannot parse size %r; use e.g. 500MB, 4GiB or 2048" % text)
    value, unit = float(m.group(1)), m.group(2).upper()
    if unit not in _SIZE_UNITS:
        raise ValueError("unknown size unit %r in %r; use B, kB, MB, GB, TB or "
                         "their binary forms kiB, MiB, GiB, TiB" % (m.group(2), text))
    n = int(value * _SIZE_UNITS[unit])
    if n <= 0:
        raise ValueError("size must be positive, got %r" % text)
    return n


def add_cascade_arguments(parser) -> None:
    """Add the cascade / size-budget options to a generator's argparse parser.

    Kept here rather than in each entry point so that ltb and lwann cannot
    drift apart in defaults or wording.
    """
    parser.add_argument("--T_max", type=float, default=None, metavar="K",
                        help="Widest temperature [K] of the refinement CASCADE. "
                             "gamma_min/T_min only say how FINE the mesh must be near "
                             "mu; they say nothing about how FAR from mu it has to be "
                             "sampled, so a mesh refined for a 1 meV kernel is "
                             "certified converged while being useless at 26 meV. Give "
                             "the widest corner of the planned sweep here and the loop "
                             "runs once per kernel width, widest first. "
                             "Default: T_min (single stage, previous behaviour).")
    parser.add_argument("--gamma_max", type=float, default=None, metavar="EV",
                        help="Widest scattering rate [eV] of the refinement cascade. "
                             "See --T_max. Default: gamma_min.")
    parser.add_argument("--ladder_ratio", type=float, default=3.0,
                        help="Ratio between the kernel widths of successive cascade "
                             "stages, i.e. the step of the logarithmic ladder from "
                             "W_max down to W_min. Default: 3 (one stage step = one "
                             "refinement_factor subdivision).")
    parser.add_argument("--ladder_max_stages", type=int, default=6,
                        help="Cap on the number of cascade stages. Default: 6.")
    parser.add_argument("--stage_tol_exponent", type=float, default=0.0,
                        help="Loosen the wide cascade stages: the stage tolerance is "
                             "error_tol * (W_stage/W_min)**EXPONENT. 0 (default) uses "
                             "the same tolerance at every width, which grades the mesh "
                             "self-similarly (local energy spacing proportional to the "
                             "distance from mu) and is what the metric's (h/W)**2 "
                             "scaling asks for. Positive values trade accuracy at the "
                             "wide end for speed; the shortfall is NOT recoverable "
                             "later, since at W_min those panels carry no selectable "
                             "defect.")
    parser.add_argument("--parent_tol", type=float, default=None, metavar="TOL",
                        help="Sum-rule metric the STARTING mesh must already meet at "
                             "the widest kernel, before any refinement. Refinement "
                             "cannot repair a parent that is too coarse there -- the "
                             "wide-corner accuracy is set by the parent and points not "
                             "placed there are never recovered -- so this check "
                             "REFUSES rather than warns. Default: error_tol.")
    parser.add_argument("--skip_parent_check", action="store_true",
                        help="Bypass the starting-mesh qualification. Expert flag: the "
                             "wide-corner error is then unverified.")
    parser.add_argument("--max_output_size", type=str, default=None, metavar="SIZE",
                        help="Upper limit on the HDF5 file the refinement may produce, "
                             "e.g. '500MB', '4GB', '2.5GiB' or a plain byte count. The "
                             "size of the NEXT output is predicted from the current "
                             "file's bytes-per-k-point before the generator runs, so "
                             "nothing oversized is ever written; the run stops and "
                             "keeps the last mesh that fitted. Default: no limit.")


def validate_cascade_args(args) -> None:
    """Range-check the cascade options and set ``args.max_output_bytes``."""
    if args.ladder_ratio <= 1.0:
        raise ValueError("ladder_ratio must be greater than 1.")
    if args.ladder_max_stages < 1:
        raise ValueError("ladder_max_stages must be at least 1.")
    if args.stage_tol_exponent < 0.0:
        raise ValueError("stage_tol_exponent must be non-negative: a negative "
                         "value would demand MORE accuracy at the wide kernel "
                         "than at the target one.")
    if args.T_max is not None and args.T_max < args.T_min:
        raise ValueError("T_max must not be below T_min; the cascade ends at "
                         "the target corner.")
    if args.gamma_max is not None and args.gamma_max < args.gamma_min:
        raise ValueError("gamma_max must not be below gamma_min; the cascade "
                         "ends at the target corner.")
    if args.parent_tol is not None and args.parent_tol <= 0:
        raise ValueError("parent_tol must be positive.")
    args.max_output_bytes = (
        None if args.max_output_size is None else parse_size(args.max_output_size)
    )


def inspect_mesh(
    path:        Path,
    temperature: float,
    gamma:       float,
    target:      "float | None" = None,
) -> dict:
    """Report what a mesh can and cannot resolve at one kernel width.

    Answers two independent questions, which need two different criteria
    because they are not the same question:

    * DENSITY -- is the mesh fine enough?  Judged by the sum-rule metric at
      the given corner.  On a uniform mesh this is a good proxy: it tracks
      the transport error monotonically (cubic 0.2806/0.0326/0.0080/0.0050
      against +145/+21.6/+8.2/+1.8%).  On an ADAPTIVE mesh it is not, which
      is why the qualifier only trusts it for regular grids.
    * ASPECT RATIO -- are the axes divided in the right proportion?  The
      metric is blind to this, ranking a 32x32 mesh above a 64x8 one that is
      250x more accurate, because mesh anisotropy is a k-space property and
      the metric lives on the energy axis.  Judged instead by
      :func:`structure.meshrefine.axis_anisotropy`.

    Returns a dict with the raw numbers; formatting is the caller's job.
    """

    with h5py.File(path, 'r') as h5file:
        data = load_band_data(h5file)
        kvec  = np.asarray(h5file['.unitcell/kvec'][()]).reshape(3, 3)
        dims  = np.asarray(h5file['.unitcell/dims'][()]).astype(bool).ravel()[:3]
        symop = (np.asarray(h5file['.unitcell/symop'][()])
                 if '.unitcell/symop' in h5file else None)
        nk = tuple(int(np.asarray(h5file['.kmesh/%s' % k][()]).ravel()[0])
                   for k in ('nkx', 'nky', 'nkz'))
        uniform = bool(np.asarray(h5file['.kmesh/irreducible'][()]).ravel()[0]) \
            if '.kmesh/irreducible' in h5file else False
        moments = (np.real(np.asarray(h5file['momentsDiagonal'][()]))
                   if 'momentsDiagonal' in h5file else None)

    band_axis = build_band_axis(data.energies)
    # Judge convergence on the refinable component: the floor is the kernel
    # weight outside the sampled band range and no mesh can reduce it.
    total, error, floor = metric_components(band_axis, temperature, gamma)

    # The coarse end of the metric-versus-density curve is available for free
    # by decimating this mesh, so the extrapolation exponent is measured rather
    # than assumed.  Only meaningful on a regular grid.
    exponent, fitted, n_fit = 1.0, error, 0
    if uniform:
        exponent, _raw, fitted, n_fit = estimate_metric_scaling(
            data.k_points, data.energies, nk, temperature, gamma,
        )

    result = {
        'path': str(path), 'nk': nk, 'nkp': int(data.k_points.shape[0]),
        'uniform_grid': uniform, 'dims': dims,
        'width': kernel_width(temperature, gamma),
        'temperature': temperature, 'gamma': gamma,
        'error': error, 'total': total, 'floor': floor,
        'target': target, 'exponent': exponent,
        'fitted': fitted, 'n_fit': n_fit,
        'raw_ratio': None, 'suggested_ratio': None, 'anisotropy': None,
    }

    if moments is not None:
        d = axis_anisotropy(data.energies, moments, data.weights, kvec, dims,
                            temperature, gamma, symop=symop)
        raw, sug = suggest_mesh_ratio(d, dims)
        result['anisotropy'] = d
        result['raw_ratio'] = raw
        result['suggested_ratio'] = sug
    else:
        logger.warning(
            "%s has no momentsDiagonal dataset, so the aspect ratio cannot be "
            "estimated; only the density check is reported.", path,
        )

    if target is not None:
        multiple = infer_mesh_multiple(nk, dims)
        # A second, coarser evaluation is not available here, so the
        # saturation floor is left to the caller when it has one.
        # Extrapolate from the fitted trend, not from this mesh's raw metric:
        # the metric is not monotone in density, so a lucky mesh would
        # otherwise pull the recommendation far too low.
        rec, message = recommend_density(
            nk, fitted, target, dims, multiple=multiple, exponent=exponent,
        )
        result.update({'recommended_nk': rec, 'message': message,
                       'multiple': multiple, 'restructured_nk': None})

        # Second recommendation at the suggested proportions, when the mesh is
        # badly proportioned.  Same k-point count -- the metric is inversely
        # proportional to the total, however it is shared out -- but far better
        # transport accuracy, which is what the ratio actually buys.
        sug = result.get('suggested_ratio')
        if sug is not None:
            act = [i for i in range(3) if dims[i] and nk[i] > 1]
            if act:
                have = np.array([float(nk[i]) for i in act]); have /= have.min()
                want = np.array([float(sug[i]) for i in act]); want /= want.min()
                if np.any(np.abs(np.log(want / np.maximum(have, 1e-30)))
                          > np.log(1.5)):
                    alt, _msg = recommend_density(
                        nk, fitted, target, dims, multiple=multiple,
                        exponent=exponent, ratio=sug,
                    )
                    if tuple(alt) != tuple(rec):
                        result['restructured_nk'] = tuple(alt)
    return result


_AXIS_NAMES = ("kx", "ky", "kz")


def _floor_flag(info: dict) -> str:
    """Inline marker for a floor large enough that the user should know."""
    floor, target = info['floor'], info.get('target')
    if target is not None and floor >= float(target):
        return "   <-- EXCEEDS the target"  # total can never reach it
    if target is not None and floor >= 0.5 * float(target):
        return "   <-- comparable to the target"
    if floor > info['error']:
        return "   <-- dominates the metric"
    return ""


def _floor_advice(info: dict) -> "list[str]":
    """Explain a substantial floor and name the knob that reduces it.

    The comparison that matters is floor against TARGET, not floor against the
    refinable part: a floor of 0.0107 against a 5e-3 target is worth saying
    even when the refinable part is momentarily larger, because the user will
    refine that away and then find the total metric parked at 0.0107 with no
    explanation ever offered.
    """
    floor, target = info['floor'], info.get('target')
    if target is None:
        if floor <= info['error']:
            return []
        threshold_msg = "is the larger part of the metric on this mesh"
    elif floor >= float(target):
        threshold_msg = ("is %.1fx the target, so the TOTAL metric cannot go "
                         "below it however many k-points are added"
                         % (floor / float(target)))
    elif floor >= 0.5 * float(target):
        threshold_msg = ("is %.0f%% of the target, so the TOTAL metric "
                         "saturates just under it and a finer mesh buys little"
                         % (100.0 * floor / float(target)))
    else:
        return []

    return [
        "  The mesh-independent floor %s." % threshold_msg,
        "     It is the kernel weight lying beyond the band range: df/da "
        "integrates to one over infinite energy, but the bands span a finite "
        "interval. NO setting reduces it -- it is fixed by the bandwidth "
        "relative to the kernel width (a 1.2 eV band gives 0.0107 where "
        "graphene's 6 eV gives 0.0002). --energy_window does not affect it; "
        "that knob enters hotspot selection only. This is why convergence "
        "above is judged on the refinable part alone.",
    ]


def format_inspection(info: dict) -> str:
    """Human-readable rendering of :func:`inspect_mesh`.

    Inactive dimensions are named rather than shown as a bare ``1``: a "1 : 1 :
    1" ratio on a 2D system invites the reading that kz carries a meaningful
    ratio, when it simply has a single division.  The density recommendation
    is the actionable line, so it goes last.
    """

    nk   = info['nk']
    dims = np.asarray(info['dims']).astype(bool).ravel()[:3]
    active = [i for i in range(3) if dims[i]]
    idle   = [i for i in range(3) if not dims[i]]
    idle_note = ("   (%s inactive: 1 division)"
                 % ", ".join(_AXIS_NAMES[i] for i in idle)) if idle else ""

    lines = [
        "Mesh:            %s" % info['path'],
        "  divisions:     %d x %d x %d  (%d k-points in file, %s)"
        % (nk[0], nk[1], nk[2], info['nkp'],
           "regular grid" if info['uniform_grid'] else "not a regular grid"),
        "  kernel:        T = %.6g K, gamma = %.4e eV  ->  W = %.4e eV"
        % (info['temperature'], info['gamma'], info['width']),
        "  sum-rule metric at this width: %.6f" % info['total'],
        "     refinable by a finer mesh: %.6f%s"
        % (info['error'],
           "   (trend %.6f, from n^-%.0f fitted over %d densities)"
           % (info['fitted'], info['exponent'], info['n_fit'])
           if info.get('n_fit') else ""),
        "     mesh-independent floor:    %.6f%s"
        % (info['floor'], _floor_flag(info)),
    ]

    if info.get('raw_ratio') is not None:
        raw, sug = info['raw_ratio'], info['suggested_ratio']
        aniso = info['anisotropy']
        isotropic = np.allclose([sug[i] for i in active], sug[active[0]]) \
            if active else True
        lines += [
            "  axis anisotropy (RMS dE per fractional step):  %s%s"
            % ("   ".join("%s %.4f" % (_AXIS_NAMES[i], aniso[i]) for i in active),
               idle_note),
            "  raw axis ratio:       %s"
            % " : ".join("%s %.2f" % (_AXIS_NAMES[i], raw[i]) for i in active),
            "  suggested (rounded):  %s%s"
            % (" : ".join("%s %d" % (_AXIS_NAMES[i], round(sug[i])) for i in active),
               "   (the axes are equivalent; nothing to change)"
               if isotropic else ""),
        ]
        # The rationale is only worth printing when there is something to act
        # on; on an isotropic mesh it is noise.
        if not isotropic:
            lines.append(
                "  Rounding up is the cheap direction to err in: on an "
                "orthorhombic test model the raw ratio 6.97 sat right at the "
                "onset of the accuracy plateau (4.6:1 gave +0.50%, 8:1 gave "
                "+0.004%), and the plateau extended to at least 32:1, while "
                "1:1 gave +12.3% at the same cost."
            )

    lines += _floor_advice(info)

    # The verdict and the recommendation go last: they are what the user acts on.
    if info.get('target') is not None:
        lines.append("")
        if info['error'] <= info['target']:
            lines.append("  target %.4e -> OK" % info['target'])
        else:
            rec = info['recommended_nk']
            lines += [
                "  target %.4e -> TOO COARSE" % info['target'],
                "  ==> suggested divisions:  %d x %d x %d"
                % (rec[0], rec[1], rec[2]),
                "      %s" % info['message'],
            ]
            if info.get('multiple', 1) > 1:
                lines.append("      kept commensurate with a multiple of %d, "
                             "inherited from the current mesh" % info['multiple'])
            # The first recommendation scales the CURRENT proportions, since
            # inspect-mesh will not silently restructure a mesh the user chose.
            # When the proportions are off, offer the restructured mesh too.
            alt = info.get('restructured_nk')
            if alt is not None:
                sug = info['suggested_ratio']
                have = np.array([float(nk[i]) for i in active]); have /= have.min()
                lines += [
                    "      (this keeps the current %s proportions, %s)"
                    % (":".join(_AXIS_NAMES[i] for i in active),
                       " : ".join("%.2g" % v for v in have)),
                    "  ==> or, at the suggested %s ratio:  %d x %d x %d"
                    % (" : ".join("%d" % round(sug[i]) for i in active),
                       alt[0], alt[1], alt[2]),
                    "      same total k-point count, but proportioned to the "
                    "band structure. The metric target needs the same number of "
                    "k-points either way; what the ratio buys is transport "
                    "accuracy, which on a test model was 250x better at matched "
                    "cost.",
                ]
    return "\n".join(lines)


def qualify_parent_mesh(params: RefinementParams, path: Path) -> bool:
    """Check a starting mesh at the WIDE corner before refining it.

    Refinement cannot repair a parent that is too coarse for the widest
    kernel of the planned sweep: the accuracy at that corner is set by the
    parent and is not recoverable later, because at W_min those panels carry
    no selectable defect.  Measured on graphene, tightening error_tol by 100x
    moved the 300 K result by 0.02% while changing the parent from 48 to 300
    moved it from -12.3% to -1.0%.  On the cubic model refinement made the
    observable WORSE than its own coarse parent (+21.6% -> -29.9%).

    Returns True when refinement should proceed.
    """

    if params.T_max is None and params.gamma_max is None:
        return True                      # no wide corner declared, nothing to check

    T_w = params.T_max if params.T_max is not None else params.T_min
    g_w = params.gamma_max if params.gamma_max is not None else params.gamma_min
    target = params.parent_tol if params.parent_tol is not None else params.error_tol

    info = inspect_mesh(path, T_w, g_w, target=target)

    logger.info("Parent mesh qualification at the widest kernel:\n%s",
                format_inspection(info))

    # ── Aspect ratio: warn loudly, never refuse ───────────────────────────
    sug = info.get('suggested_ratio')
    raw = info.get('raw_ratio')
    if sug is not None:
        nk = np.array(info['nk'], dtype=float)
        dims = info['dims']
        active = dims & (nk > 0)
        if np.any(active):
            have = nk[active] / np.min(nk[active])
            want = np.asarray(sug)[active]
            if np.any(want / np.maximum(have, 1e-30) > 1.5) or \
               np.any(have / np.maximum(want, 1e-30) > 1.5):
                logger.warning(
                    "ASPECT RATIO: this mesh is divided %s but the band "
                    "structure suggests %s (raw %s). The sum-rule metric "
                    "CANNOT see this -- it ranked a 32x32 mesh above a 64x8 "
                    "one that was 250x more accurate -- so the density check "
                    "below may pass while the mesh is still badly "
                    "proportioned. Refinement will proceed, but run "
                    "'inspect-mesh' and consider regenerating the starting "
                    "mesh with the suggested proportions.",
                    " : ".join("%d" % v for v in np.array(info['nk'])),
                    " : ".join("%d" % round(v) for v in sug),
                    " : ".join("%.2f" % v for v in raw),
                )

    # ── Density: refuse ───────────────────────────────────────────────────
    if info['error'] <= target:
        return True

    if not info['uniform_grid']:
        logger.warning(
            "The starting mesh is not a regular grid, and the sum-rule metric "
            "is only a reliable density proxy on regular grids, so the wide "
            "corner cannot be qualified (metric %.6f against target %.4e). "
            "Proceeding, but the result at the wide corner is unverified.",
            info['error'], target,
        )
        return True

    rec = info.get('recommended_nk')
    logger.error(
        "STARTING MESH TOO COARSE at the widest kernel: metric %.6f "
        "against the target %.4e at W = %.4e eV. Refinement cannot repair "
        "this -- the accuracy at the wide corner is set by the parent, and "
        "points not placed there are never recovered, because at the "
        "narrow corner those panels carry no selectable defect. "
        "Regenerate the starting mesh with about %d x %d x %d divisions "
        "and try again, or raise --parent_tol if you accept the error.",
        info['error'], target, info['width'], rec[0], rec[1], rec[2],
    )
    return False


def format_bytes(n: float) -> str:
    """Human-readable byte count (binary prefixes)."""
    for unit in ("B", "kiB", "MiB", "GiB", "TiB"):
        if abs(n) < 1024.0 or unit == "TiB":
            return f"{n:.1f} {unit}" if unit != "B" else f"{n:.0f} B"
        n /= 1024.0
    return f"{n:.1f} TiB"                              # pragma: no cover


def predict_output_bytes(current_hdf5: Path, n_before: int, n_after: int) -> float:
    """Estimate the size of the next generator output.

    The generator writes one ``moments`` and one ``momentsBfield`` dataset per
    k-point (KI-02), so the file size is very nearly linear in the k-point
    count with a small fixed header.  Scaling the current file's
    bytes-per-k-point is therefore accurate to a few percent, which is all a
    budget check needs.  Returns 0.0 if the current size cannot be read.
    """
    try:
        size = float(current_hdf5.stat().st_size)
    except OSError as exc:                             # pragma: no cover
        logger.debug("Cannot stat %s for the size estimate: %s", current_hdf5, exc)
        return 0.0
    if n_before <= 0:
        return 0.0
    return size * float(n_after) / float(n_before)


# ---------------------------------------------------------------------------
# Shared refinement loop
# ---------------------------------------------------------------------------

def read_source_dims(path: Path):
    """Return ``.unitcell/dims`` of *path* as a length-3 bool array, or None.

    Generators pass this on to setCustomKmesh so that a refined mesh inherits
    the dimensionality of the calculation it was derived from instead of
    re-deriving it.  Re-derivation is not always possible: a layered system
    can have a three-dimensional k-mesh while the Wannier projection is
    effectively two-dimensional.
    """
    try:
        with h5py.File(path, "r") as f:
            if "/.unitcell/dims" not in f:
                return None
            dims = np.asarray(f["/.unitcell/dims"][()], dtype=bool).ravel()
        return dims if dims.shape == (3,) else None
    except Exception as exc:                      # pragma: no cover
        logger.warning("Could not read /.unitcell/dims from %s: %s", path, exc)
        return None


def read_source_symmetry(path: Path):
    """Return ``(wedge, symop)`` of *path*.

    ``wedge`` is True when the stored k-points cover only an irreducible part
    of the Brillouin zone, which is the case for

      * a regular irreducible grid from ``ltb run``   -> ``.kmesh/irreducible``
      * a refined custom mesh derived from one        -> ``.kmesh/symmetrized``

    The second marker is needed because a custom mesh must be written with
    ``.kmesh/irreducible = False`` (its weights are non-uniform and cannot be
    reconstructed from ``multiplicity/(nkx*nky*nkz)``, which is the only thing
    linretrace uses that flag for).  Without it a continuation run
    (``ltb refine refined_iter_N.hdf5 ...``) would silently drop the
    symmetrisation of the optical elements.

    ``symop`` is the ``.unitcell/symop`` array, returned only when the mesh is
    a wedge and carries more than the identity; otherwise None.  Generators
    pass it on to setCustomKmesh, which then makes the evaluation symmetrise
    the optical / B-field moments over the star of every k-point.
    """
    try:
        with h5py.File(path, "r") as f:
            irr = bool(f["/.kmesh/irreducible"][()]) if "/.kmesh/irreducible" in f else False
            sym = bool(f["/.kmesh/symmetrized"][()]) if "/.kmesh/symmetrized" in f else False
            wedge = irr or sym
            if not wedge or "/.unitcell/symop" not in f:
                return wedge, None
            symop = np.asarray(f["/.unitcell/symop"][()], dtype=float)
            irr = wedge
    except Exception as exc:                      # pragma: no cover
        logger.warning("Could not read symmetry information from %s: %s", path, exc)
        return False, None

    if symop.ndim != 3 or symop.shape[1:] != (3, 3) or symop.shape[0] < 1:
        logger.warning("%s has a malformed /.unitcell/symop (shape %s); "
                       "treating the mesh as reducible.", path, symop.shape)
        return irr, None
    return irr, symop


def _reject_unsymmetrised_irreducible(path: Path, generator) -> bool:
    """True (and logs why) if *path* is irreducible but *generator* cannot
    symmetrise.

    An irreducible mesh covers only a wedge.  Band energies are symmetry
    invariant and come out right either way, but the optical matrix elements
    v_i v_j do not: summed raw over the wedge they break the symmetry of the
    Onsager tensor (cubic single-band test model, 8x8x8: the wedge sum gives
    (0.2148, 0.3203, 0.2148) against the correct (0.25, 0.25, 0.25)).  The
    failure is silent, because the refinement error metric only sees energies.

    Generators that accept a ``symop`` argument take the symmetrising path and
    are fine.  Anything else must refuse.
    """
    with h5py.File(path, "r") as f:
        irr = "/.kmesh/irreducible" in f and bool(f["/.kmesh/irreducible"][()])
        sym = "/.kmesh/symmetrized" in f and bool(f["/.kmesh/symmetrized"][()])
        if not (irr or sym):
            return False

    if getattr(generator, "symop", None) is not None:
        return False

    logger.error(
        "%s holds an IRREDUCIBLE k-mesh, but %s was constructed without "
        "symmetry operations, so the refined mesh would be evaluated through "
        "the reducible code path. That path does not symmetrise the optical "
        "matrix elements over the star of each k-point, and the resulting "
        "transport tensor would be wrong while the band energies -- and hence "
        "the refinement error metric -- still looked correct. Either supply "
        "/.unitcell/symop (read_source_symmetry) or regenerate the coarse "
        "mesh as a reducible grid.",
        path, type(generator).__name__,
    )
    return True


def run_refinement(params: RefinementParams, generator: MeshGenerator) -> int:
    """
    Iteratively refine the k-mesh stored in *params.initial_hdf5*.

    On each iteration:
      1. Load band data (energies, k-points, weights, cell_deltas) from the
         current HDF5 file.
      2. Evaluate the df/da integration error.
      3. If the error is below *params.error_tol*, stop.
      4. Detect hotspot k-points near *params.chemical_potential*.
      5. Refine those k-points (subdivide) while preserving total weight.
         Cell widths are taken from the stored cell_deltas array, so the
         subdivision is correct regardless of whether the starting mesh is
         reducible or irreducible and regardless of which iteration created
         each point.
      6. Call *generator.generate()* to produce the next HDF5 file.
      7. Patch cell_deltas and provenance metadata from the mesh file into
         the generator output so that they are available on the next
         iteration or a later continuation run (_patch_refinement_metadata).
      8. Optionally clean up intermediate files (only files created by this
         run are ever deleted; the input file is never touched).

    Continuation
    ------------
    *params.initial_hdf5* may itself be the output of a previous refinement
    run (e.g. ``refined_iter_3.hdf5``).  This is detected automatically from
    the ``refinement_iteration`` attribute / the ``cell_deltas`` dataset:
    per-point cell widths are then read from the file (NOT re-derived from
    coarse grid steps), and output numbering resumes after the input's
    iteration index.  Output files never overwrite existing files: if a name
    is taken, the index is bumped past it with a warning.

    Parameters
    ----------
    params : RefinementParams
    generator : MeshGenerator

    Returns
    -------
    int
        0 on success, 1 on error (suitable as a process exit code).
    """
    assert isinstance(generator, MeshGenerator), (
        f"generator must implement the MeshGenerator protocol, got {type(generator)}"
    )

    workdir = params.workdir.resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    # Make the energy scale the metric actually resolves explicit: the df/da
    # peak at mu has a half-width of roughly max(gamma, kB*T).  The mesh must
    # provide several distinct band energies inside that width, otherwise the
    # sum rule cannot be satisfied no matter how the hotspots are chosen.
    thermal_width = kB_eV * params.T_min
    logger.info(
        "df/da metric: T_min = %.6g K (kB*T = %.4e eV), gamma_min = %.4e eV, "
        "beta = %.4e eV^-1 -> peak half-width ~ %.4e eV.",
        params.T_min, thermal_width, params.gamma_min,
        1.0 / thermal_width if thermal_width > 0 else float('inf'),
        max(thermal_width, params.gamma_min),
    )
    # T_min/gamma_min say how FINE the mesh must be near mu, and nothing about
    # how FAR from mu it must be sampled; without a wide corner the loop
    # certifies a mesh that is useless at any larger kernel width.
    if params.T_max is None and params.gamma_max is None:
        logger.info(
            "No cascade requested: the mesh will be certified at W = %.4e eV "
            "only. If the planned sweep reaches higher temperatures or "
            "scattering rates, give --T_max / --gamma_max so that the refinement "
            "also samples the wider kernel.",
            max(thermal_width, params.gamma_min),
        )

    # ── Qualify the starting mesh at the WIDE corner, before refining ─────
    if params.skip_parent_check:
        if params.T_max is not None or params.gamma_max is not None:
            logger.warning(
                "Parent qualification skipped at the user's request. The "
                "accuracy at the wide corner is set by the starting mesh and "
                "refinement cannot repair it.",
            )
    elif not qualify_parent_mesh(params, params.initial_hdf5.resolve()):
        return 1

    current_hdf5: Path            = params.initial_hdf5.resolve()
    final_error:  Optional[float] = None
    previous_error: Optional[float] = None  # error of the preceding iteration (plateau check)

    # files this run created and may therefore delete; the input file is
    # never in this set and is consequently never deleted or overwritten
    created_by_run: set = set()

    # ── Continuation detection ────────────────────────────────────────────
    try:
        is_continuation, last_iter = _detect_refinement_state(current_hdf5)
    except Exception as exc:
        logger.error("Cannot inspect %s: %s", current_hdf5, exc)
        return 1

    # Both fresh coarse input and a continuation input may be an irreducible
    # wedge (the latter carries .kmesh/symmetrized), so the check applies to
    # both: a wedge evaluated through the non-symmetrising path gives a wrong
    # transport tensor while the band energies still look correct.
    try:
        if _reject_unsymmetrised_irreducible(current_hdf5, generator):
            return 1
    except Exception as exc:
        logger.error("Cannot inspect %s: %s", current_hdf5, exc)
        return 1
    # resume after the input's iteration index, and never number below
    # refinement files already present in the working directory
    next_index = max(
        (last_iter + 1) if last_iter is not None else 1,
        _max_existing_index(workdir) + 1,
    )

    # ── Cascade over kernel widths, WIDEST FIRST ──────────────────────────
    ladder = kernel_ladder(
        params.T_min, params.gamma_min, params.T_max, params.gamma_max,
        ratio=params.ladder_ratio, max_stages=params.ladder_max_stages,
    )
    tolerances = stage_tolerances(ladder, params.error_tol,
                                  params.stage_tol_exponent)

    if len(ladder) > 1:
        logger.info(
            "Refinement cascade: %d stages from W = %.4e eV down to %.4e eV "
            "(ratio %.3g). Widest first, so the wide stage places the parents "
            "that the narrow stages subdivide.",
            len(ladder), kernel_width(*ladder[0]), kernel_width(*ladder[-1]),
            params.ladder_ratio,
        )
        for j, ((t, g), tol) in enumerate(zip(ladder, tolerances)):
            logger.info(
                "  stage %d/%d: T = %.6g K, gamma = %.4e eV -> W = %.4e eV, "
                "error_tol = %.4e", j + 1, len(ladder), t, g,
                kernel_width(t, g), tol,
            )

    first_stage = True
    stop_cascade = False

    for stage_index, ((stage_T, stage_gamma), stage_tol) in enumerate(
        zip(ladder, tolerances)
    ):
        if stop_cascade:
            break

        stage_width = kernel_width(stage_T, stage_gamma)
        if len(ladder) > 1:
            logger.info(
                "=== Cascade stage %d/%d: T = %.6g K, gamma = %.4e eV "
                "(W = %.4e eV), error_tol = %.4e ===",
                stage_index + 1, len(ladder), stage_T, stage_gamma,
                stage_width, stage_tol,
            )

        # the plateau test compares consecutive iterations of the SAME stage
        previous_error = None

        for iteration in range(params.max_iter):
            logger.info("--- Iteration %d ---", iteration)

            # ── 1. Load band data ─────────────────────────────────────────────
            try:
                with h5py.File(current_hdf5, 'r') as h5file:
                    data: BandData = load_band_data(h5file)
                    band_axis = build_band_axis(data.energies)

                    if iteration == 0 and first_stage:
                        first_stage = False
                        if is_continuation:
                            assert data.cell_deltas_source == "file", (
                                "continuation input must provide /.kmesh/cell_deltas"
                            )
                            logger.info(
                                "Initial mesh is an already-refined custom mesh "
                                "(%d k-points; produced by refinement iteration %s). "
                                "Continuing refinement: per-point cell widths read "
                                "from /.kmesh/cell_deltas, output numbering resumes "
                                "at %d.",
                                data.k_points.shape[0],
                                str(last_iter) if last_iter is not None else "unknown",
                                next_index,
                            )
                        else:
                            irreducible = bool(h5file['/.kmesh/irreducible'][()])
                            mesh_type   = "irreducible" if irreducible else "reducible"
                            logger.info(
                                "Initial mesh is %s (%d k-points). "
                                "Cell widths initialised from coarse grid steps 1/nk_i.",
                                mesh_type, data.k_points.shape[0],
                            )
            except ValueError as exc:
                # e.g. custom-mesh file without cell_deltas (see load_band_data)
                logger.error("%s", exc)
                return 1

            # ── 2. Evaluate error ─────────────────────────────────────────────
            # Convergence is judged on the REFINABLE component.  The rest of
            # the metric is kernel weight lying beyond the band range, fixed by
            # the bandwidth relative to the kernel width, and no number of
            # k-points reduces it -- so testing the total against the tolerance
            # asks the loop to move something it cannot.  It is also what
            # select_refinement_panels already acts on, so this makes the
            # stopping criterion agree with the selection criterion.  On an
            # orthorhombic model with a 1.2 eV bandwidth the floor is 0.0107
            # against a 5e-3 tolerance, and a mesh whose refinable defect had
            # fallen to 0.0011 -- transport error +0.19%, i.e. converged -- was
            # reported as a failure.
            try:
                total_error, final_error, floor_error = metric_components(
                    band_axis, stage_T, stage_gamma)
            except MissingDependencyError as exc:
                logger.error("%s", exc)
                return 1

            floor_note = ("" if floor_error < 0.1 * stage_tol
                          else "  [+ %.6f mesh-independent floor]" % floor_error)
            if len(ladder) > 1:
                logger.info("Stage %d/%d iteration %d: error = %.6f (W = %.4e eV)%s",
                            stage_index + 1, len(ladder), iteration, final_error,
                            stage_width, floor_note)
            else:
                logger.info("Iteration %d: error = %.6f%s",
                            iteration, final_error, floor_note)

            # ── 3. Convergence check ──────────────────────────────────────────
            if final_error <= stage_tol:
                logger.info(
                    "Target error reached (%.6f <= %.6f)%s.%s", final_error, stage_tol,
                    " at W = %.4e eV" % stage_width if len(ladder) > 1 else "",
                    ("  The total metric is %.6f; the difference is kernel weight "
                     "beyond the band range, which no mesh can reduce."
                     % total_error) if floor_note else "",
                )
                break

            # ── 3b. Error-plateau check ───────────────────────────────────────
            # The refinement error typically plateaus after a few steps at a
            # value dictated by the resolution of the initial mesh.  Stop once an
            # iteration reduced the error by less than plateau_tol (relative);
            # active from the 4th refinement step onward only, since the error
            # may equilibrate non-monotonically in the early stages.
            #
            # Two distinct events land here and they must not be reported alike:
            #
            #   stalled  -- 0 <= improvement < plateau_tol.  The mesh is being
            #               refined but the metric no longer responds.
            #   backward -- improvement < 0, i.e. the metric went UP.  This does
            #               NOT mean the mesh got worse: refinement only ever ADDS
            #               k-points, so each iterate is a strict superset of its
            #               predecessor and is strictly better as a BZ sampling.
            #               The sum-rule metric is |1 - int df/da|, a SIGNED sum of
            #               per-panel trapezoid defects, and those defects carry
            #               opposite signs in the peak and in the tails.  An
            #               iterate whose peak defect happens to cancel part of the
            #               tail defect scores better than the finer mesh that
            #               follows it.  Reporting that as a regression -- or worse,
            #               keeping the coarser mesh -- would be wrong.  The final
            #               mesh written is always the last (finest) one.
            if (
                params.plateau_tol > 0
                and iteration >= PLATEAU_MIN_ITER
                and previous_error is not None
                and previous_error > 0
                and (previous_error - final_error) < params.plateau_tol * previous_error
            ):
                improvement = 100.0 * (previous_error - final_error) / previous_error
                if final_error > previous_error:
                    logger.warning(
                        "Sum-rule metric increased at iteration %d (%.6f -> %.6f). "
                        "The mesh itself did not get worse -- refinement only adds "
                        "k-points, so this iterate contains every point of the "
                        "previous one. The metric is a signed sum of per-panel "
                        "quadrature defects whose peak and tail contributions have "
                        "opposite signs, so an accidental cancellation at the "
                        "previous iterate can score better than a strictly finer "
                        "mesh. Stopping here: the metric can no longer steer the "
                        "refinement.",
                        iteration, previous_error, final_error,
                    )
                else:
                    logger.warning(
                        "Error plateau detected: iteration %d reduced the error by "
                        "only %.2f%% (%.6f -> %.6f), less than the plateau threshold "
                        "of %.2f%%.",
                        iteration, improvement, previous_error, final_error,
                        100.0 * params.plateau_tol,
                    )
                if stage_index + 1 < len(ladder):
                    # An intermediate stage is only there to place parents for
                    # the narrower ones; a plateau at a wide kernel is not a
                    # failure of the run, so hand over instead of alarming.
                    logger.warning(
                        "Stage %d/%d (W = %.4e eV) plateaued at %.6f above its "
                        "tolerance %.6f; handing the mesh to the next, narrower "
                        "stage.",
                        stage_index + 1, len(ladder), stage_width, final_error,
                        stage_tol,
                    )
                else:
                    logger.warning(
                        "Target precision (%.6f) NOT reached; final metric %.6f on %d "
                        "k-points. The residual is dominated by the trapezoid error in "
                        "the kernel TAILS (|e-mu| well beyond the %.2e eV kernel width), "
                        "which the hotspot window does not target. Options: start from a "
                        "finer initial mesh, widen --energy_window, relax --error_tol, or "
                        "disable this check with --plateau_tol 0.",
                        stage_tol, final_error, band_axis.shape[0],
                        max(kB_eV * stage_T, stage_gamma),
                    )
                break
            previous_error = final_error

            # ── 4. Detect hotspots ────────────────────────────────────────────
            # Selection and metric are the same quantity: the per-panel quadrature
            # defect of the df/da sum rule.  See select_refinement_panels for why
            # the previous pointwise criterion could not reach the kernel tails.
            try:
                marked, panel_info = select_refinement_panels(
                    band_axis,
                    stage_T,
                    stage_gamma,
                    params.chemical_potential,
                    energy_window=params.energy_window,
                    defect_fraction=params.defect_fraction,
                )
            except MissingDependencyError as exc:
                logger.error("%s", exc)
                return 1

            if not marked.any():
                logger.warning(
                    "No band-axis panel inside the energy window still carries a "
                    "significant quadrature defect, so additional k-points cannot "
                    "reduce the error further (residual %.6f is set by the band "
                    "range / initial mesh). Stopping.", final_error,
                )
                break

            hotspot_mask = hotspot_mask_from_panels(data.energies, band_axis, marked)

            logger.info(
                "Hotspot panels: %d of %d carry %.1f%% of the defect (%.4e of "
                "%.4e); reach %.4e eV from mu (ceiling energy_window=%.4e, kernel "
                "half-width=%.4e); %d of %d k-points selected.",
                panel_info["n_marked"], panel_info["n_panels"],
                100.0 * panel_info["captured"] / panel_info["total"]
                if panel_info["total"] > 0 else 0.0,
                panel_info["captured"], panel_info["total"], panel_info["reach"],
                params.energy_window,
                max(kB_eV * stage_T, stage_gamma),
                int(hotspot_mask.sum()), hotspot_mask.size,
            )

            # ── 5. Refine mesh ────────────────────────────────────────────────
            refined_points, refined_weights, refined_deltas = refine_kmesh(
                data,
                target_energy=params.chemical_potential,
                tolerance=0.0,
                refinement_factor=params.refinement_factor,
                hotspot_mask=hotspot_mask,
            )

            n_before = data.k_points.shape[0]
            n_after  = refined_points.shape[0]

            if n_after == n_before:
                logger.warning(
                    "Refinement did not modify the mesh; stopping to avoid "
                    "infinite loop."
                )
                break

            predicted = predict_output_bytes(current_hdf5, n_before, n_after)
            logger.info("Mesh size: %d → %d k-points (next output ~%s).",
                        n_before, n_after, format_bytes(predicted))

            # ── 5b. Output-size budget ────────────────────────────────────────
            # The generator writes a moments dataset per k-point (KI-02), so the
            # file grows linearly in the mesh and a refinement chasing an
            # extended Fermi surface can reach tens of GB without warning.  The
            # check happens BEFORE the generator runs, so nothing oversized is
            # ever written; the last mesh that fitted is kept and reported.
            if (
                params.max_output_bytes is not None
                and predicted > float(params.max_output_bytes)
            ):
                logger.warning(
                    "Refining to %d k-points would produce an HDF5 of about %s, "
                    "over the budget of %s set by --max_output_size. Stopping "
                    "here and keeping the current %d-k-point mesh (metric %.6f "
                    "at W = %.4e eV). Raise the budget, relax --error_tol, "
                    "lower --defect_fraction to grow the mesh more slowly, or "
                    "narrow the cascade with --T_max / --gamma_max.",
                    n_after, format_bytes(predicted),
                    format_bytes(float(params.max_output_bytes)),
                    n_before, final_error, stage_width,
                )
                stop_cascade = True
                break

            before_weight = np.sum(data.weights)
            after_weight  = np.sum(refined_weights)
            if not np.allclose(before_weight, after_weight, rtol=1e-10, atol=1e-12):
                logger.warning(
                    "Weight conservation check failed: before=%s after=%s",
                    before_weight, after_weight,
                )

            # ── 6. Write mesh file and call generator ─────────────────────────
            # Collision-free naming: never overwrite existing files (outputs of
            # earlier runs, or the input file of a continuation run in the same
            # working directory).  next_index resumes after the input's iteration
            # for continuation runs.
            mesh_path, output_path, next_index = _reserve_iteration_paths(
                workdir, next_index, forbidden={params.initial_hdf5.resolve()},
            )

            write_custom_mesh(str(mesh_path), refined_points, refined_weights,
                              refined_deltas)
            created_by_run.add(mesh_path)

            try:
                generator.generate(refined_points, refined_weights, str(output_path))
            except Exception as exc:
                logger.error("Generator failed on iteration %d: %s", iteration, exc)
                return 1
            created_by_run.add(output_path)

            # ── 7. Patch cell_deltas and provenance into the generator output ─
            # generator.generate() calls h5output which does not know about
            # cell_deltas.  Copy it from mesh_path into output_path so that the
            # next iteration's (or a later continuation run's) load_band_data
            # finds it directly and never falls back to the nkx/nky/nkz path
            # (which always gives 1 for custom meshes).  Also stamp the global
            # iteration index so continuation runs resume the numbering.
            try:
                _patch_refinement_metadata(output_path, mesh_path,
                                           iteration_index=next_index,
                                           parent=current_hdf5)
            except Exception as exc:
                logger.error("Failed to patch refinement metadata into %s: %s",
                             output_path, exc)
                return 1

            # ── 8. Clean up intermediate files ────────────────────────────────
            # Only files created by this run are ever deleted; the input file
            # (params.initial_hdf5) is never in created_by_run.
            if not params.keep_intermediate:
                if mesh_path.exists():
                    mesh_path.unlink()
                if current_hdf5 in created_by_run and current_hdf5.exists():
                    # output of the previous iteration of this run, now superseded
                    current_hdf5.unlink()

            current_hdf5 = output_path
            next_index  += 1

        else:
            logger.warning(
                "Maximum iterations (%d) reached%s. Final error = %.6f. "
                "max_iter is counted PER CASCADE STAGE.",
                params.max_iter,
                " in stage %d/%d (W = %.4e eV)"
                % (stage_index + 1, len(ladder), stage_width)
                if len(ladder) > 1 else "",
                final_error if final_error is not None else float('nan'),
            )

    # ── Cascade summary: the metric at every rung, on the final mesh ──────
    # The single number reported by the last stage only certifies the
    # NARROWEST kernel.  Printing the whole ladder makes the remaining
    # inadequacy at wider kernels visible instead of implicit.
    if len(ladder) > 1:
        try:
            with h5py.File(current_hdf5, 'r') as h5file:
                final_axis = build_band_axis(load_band_data(h5file).energies)
            logger.info("Final mesh evaluated at every rung of the ladder:")
            for (t, g), tol in zip(ladder, tolerances):
                _tot, err, flr = metric_components(final_axis, t, g)
                logger.info(
                    "  W = %.4e eV (T = %.6g K, gamma = %.4e eV): error = %.6f "
                    "%s tol %.4e%s", kernel_width(t, g), t, g, err,
                    "<=" if err <= tol else " >", tol,
                    "  [+ %.6f floor]" % flr if flr >= 0.1 * tol else "",
                )
        except Exception as exc:                       # pragma: no cover
            logger.debug("Could not build the cascade summary: %s", exc)

    logger.info("Refinement completed. Final mesh: %s", current_hdf5)
    return 0
