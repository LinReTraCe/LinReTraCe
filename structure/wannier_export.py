"""
structure/wannier_export.py
===========================
Export a :class:`~structure.tb.TightBinding` model produced by ``ltb run``
as a Wannier90-compatible file set so that ``lwann`` /
:class:`~structure.wannier.Wannier90Calculation` can read it without any
modification.

Three files are written:

``<stem>_hr.dat``
    Real-space Hamiltonian H(R) in the standard Wannier90 format, read by
    ``Wannier90Calculation._readHr``.

``<stem>.nnkp``
    Minimal .nnkp file containing the four sections that
    ``Wannier90Calculation._readNnkp`` reads:
    ``real_lattice``, ``recip_lattice``, ``kpoints``, ``projections``.
    Two mandatory stub sections (``nnkpts``, ``exclude_bands``) are appended
    so that the file looks valid to Wannier90 itself.

``<stem>.wout``
    One-line stub supplying the ``Length Unit`` string that
    ``Wannier90Calculation._readWout`` searches for.

Public API
----------
::

    from structure.wannier_export import write_wannier90
    write_wannier90(tb, "wann")   # → wann_hr.dat, wann.nnkp, wann.wout

The three files can be placed in a directory and passed straight to ``lwann``::

    ltb run model.tb 10 10 10 2.0 --write-hr outdir/wann
    lwann outdir/ --charge 2.0

Sign-convention note
--------------------
``TightBinding._readTb`` transforms the user's hopping parameters into
``tb.hr`` as follows:

* **On-site** (R = 0, orb1 == orb2): the stored value is
  ``hr[ir, m, m] = +E_onsite`` (two sign flips cancel).
* **All other entries**: ``hr[ir, orb1-1, orb2-1] = −t``.

``Wannier90Calculation.computeHamiltonian`` evaluates::

    H(k) = Σ_R  hr[R] · exp(i·2π·k·R) / rmultiplicity[R]

which, with all multiplicities = 1, is exactly the Wannier90 convention::

    H(k) = Σ_R  H(R) · exp(i·k·R)        (k in reduced coordinates, 2π included)

so ``tb.hr`` can be written verbatim to ``_hr.dat``.

Field-width reference (must match reader slice indices)
-------------------------------------------------------
``_hr.dat`` data lines::

    chars  0– 4  : rx  (int, %5d)
    chars  5– 9  : ry  (int, %5d)
    chars 10–14  : rz  (int, %5d)
    chars 15–19  : p1  (int, %5d, 1-based row index)
    chars 20–24  : p2  (int, %5d, 1-based col index)
    chars 25–36  : Re  (float, %12.6f)
    chars 37–    : Im  (float, %12.6f + newline)

``_readHr`` reads::

    rx,ry,rz,p1,p2 = [int(line[i*5:i*5+5]) for i in range(5)]
    tr, tj = float(line[25:37]), float(line[37:])
    matrix[ir, p1-1, p2-1] = tr + 1j*tj

``_readHr`` determines the R-point of each block from the (rx,ry,rz) of
its *last* line, so all nproj² lines in a block must carry the same R.

``_readNnkp`` reads::

    real_lattice  : float(line[:12]), float(line[12:24]), float(line[24:36])
    recip_lattice : same slices
    kpoints       : float(line[:14]), float(line[14:28]), float(line[28:42])
    projections   : float(line[:10]), float(line[10:21]), float(line[21:32])

``recip_lattice`` is stored WITHOUT the 2π factor (Wannier90 crystallographic
convention b_i / 2π in Å⁻¹); ``_readNnkp`` then divides by ``lengthscale``
(= 1.0 for Ang).  ``tb.kvec`` contains the factor 2π, so we write
``tb.kvec / (2π * lengthscale)``.
"""

from __future__ import annotations

import datetime
import logging
import math
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _datestamp() -> str:
    return datetime.datetime.now().strftime("%e%b%Y at %H:%M:%S")


def _reducible_kpoints(tb) -> np.ndarray:
    """
    Return the full reducible k-mesh as an (nkx*nky*nkz, 3) array of
    fractional coordinates, outer-x / inner-z loop order (same as
    ``TightBinding._setupKmesh`` for the reducible case).

    ``lwann`` / ``Wannier90Calculation.readData`` requires a full reducible
    grid and checks ``nkx * nky * nkz == nkp``.  If ``ltb`` was run with
    ``--red`` the grid is already reducible; otherwise we reconstruct it.
    """
    if not tb.irreducible:
        # Already a reducible grid stored in tb.kpoints
        return tb.kpoints

    # Reconstruct the plain reducible grid from the grid dimensions and the
    # kshift flag, replicating the logic in TightBinding._setupKmesh.
    nkx, nky, nkz = tb.nkx, tb.nky, tb.nkz
    kshift = getattr(tb, 'kshift', False)

    if kshift:
        dims = getattr(tb, 'dims',
                       np.array([nkx > 1, nky > 1, nkz > 1], dtype=bool))
        is_shift = np.array([int(d) for d in dims], dtype=np.float64)
    else:
        is_shift = np.zeros(3, dtype=np.float64)

    kgrid = np.array([nkx, nky, nkz], dtype=np.float64)
    kpoints = []
    for ikx in np.linspace(0, 1, nkx, endpoint=False):
        for iky in np.linspace(0, 1, nky, endpoint=False):
            for ikz in np.linspace(0, 1, nkz, endpoint=False):
                kpoints.append([ikx, iky, ikz])
    kpoints = np.array(kpoints, dtype=np.float64)
    if kshift:
        kpoints += is_shift / 2.0 / kgrid[None, :]

    return kpoints


# ---------------------------------------------------------------------------
# _hr.dat writer
# ---------------------------------------------------------------------------

def _write_hr(tb, path: Path) -> None:
    """
    Write ``<stem>_hr.dat``.

    Loop order: outer = R-point (in tb.rpoints order), middle = p1 (row,
    1-based), inner = p2 (col, 1-based).  This reproduces the standard
    Wannier90 column-major-within-R layout and matches what ``_readHr``
    expects (it reads nproj² consecutive lines per R-block and derives
    matrix[ir, p1-1, p2-1]).
    """
    nrp   = tb.nrp
    nband = tb.energyBandMax
    rpts  = tb.rpoints    # (nrp, 3) int
    hr    = tb.hr         # (nrp, nband, nband) complex128

    # All ltb R-points have multiplicity 1: they come directly from the
    # unique R-vectors listed in the hopping file with no Wigner-Seitz
    # averaging.  Wannier90Calculation.computeHamiltonian divides by
    # rmultiplicity, so writing 1 everywhere is correct.
    per_row = 15
    n_mult_rows = (nrp + per_row - 1) // per_row

    logger.info("Writing {}".format(path))
    with open(str(path), 'w') as f:
        # Line 1: date comment
        f.write(" written on {}\n".format(_datestamp()))
        # Line 2: number of Wannier functions
        f.write("{:12d}\n".format(nband))
        # Line 3: number of R-points
        f.write("{:12d}\n".format(nrp))
        # Multiplicities: 15 per row, right-aligned in 5-char fields
        for row in range(n_mult_rows):
            start = row * per_row
            end   = min(start + per_row, nrp)
            f.write("".join("{:5d}".format(1) for _ in range(end - start)) + "\n")

        # H(R) matrix elements.
        # _readHr slices: rx=line[0:5] ry=line[5:10] rz=line[10:15]
        #                 p1=line[15:20] p2=line[20:25]
        #                 tr=line[25:37] ti=line[37:]
        # _readHr stores: matrix[ir, p1-1, p2-1] = tr + 1j*ti
        # tb.hr[ir, n, m] = H_{n,m}(R)  (0-based)
        # We write p1 = n+1 (row), p2 = m+1 (col).
        # Loop: outer p1 (row), inner p2 (col) — standard Wannier90 order.
        for ir in range(nrp):
            rx = int(rpts[ir, 0])
            ry = int(rpts[ir, 1])
            rz = int(rpts[ir, 2])
            for p1 in range(1, nband + 1):        # row index (1-based)
                for p2 in range(1, nband + 1):    # col index (1-based)
                    val = hr[ir, p1 - 1, p2 - 1]
                    f.write(
                        "{:5d}{:5d}{:5d}{:5d}{:5d}{:12.6f}{:12.6f}\n".format(
                            rx, ry, rz, p1, p2, val.real, val.imag
                        )
                    )

    logger.info("  {:d} R-points, {:d} Wannier function(s).".format(nrp, nband))


# ---------------------------------------------------------------------------
# .nnkp writer
# ---------------------------------------------------------------------------

def _write_nnkp(tb, path: Path) -> None:
    """
    Write ``<stem>.nnkp``.

    ``Wannier90Calculation._readNnkp`` reads four sections sequentially,
    reopening the file for each one.  We write them in order:
    real_lattice → recip_lattice → kpoints → projections, followed by
    mandatory stubs for nnkpts and exclude_bands.

    reciprocal lattice convention
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Wannier90 stores b_i / 2π in Å⁻¹ (i.e. without the 2π factor).
    ``_readNnkp`` then does ``self.kvec /= self.lengthscale`` (= 1.0 for Ang)
    so the in-memory ``kvec`` equals the file value in Å⁻¹.
    ``computeHamiltonian`` uses ``kvec`` only via ``prefactor_r``:
    ``prefactor_r = rvec.T @ rpoints`` (real-space prefactors in Å),
    and ``rdotk = 2π * kpoints @ rpoints``.  The reciprocal lattice is used
    only to derive ``vol`` and ``ortho`` in ``_readNnkp``, where the missing
    2π does not matter.  Therefore we write ``tb.kvec / (2π)``.
    """
    rvec      = tb.rvec          # (3,3) Å, rows = vectors
    kvec_file = tb.kvec / (2.0 * math.pi)   # store without 2π
    nband     = tb.energyBandMax
    kpoints   = _reducible_kpoints(tb)
    nkp       = len(kpoints)

    # Projection centres in fractional coordinates.
    # tb.orbitals is either None or (nband, 3) with 0-based ordering.
    if tb.orbitals is not None:
        plist = tb.orbitals       # (nband, 3)
    else:
        plist = np.zeros((nband, 3), dtype=np.float64)

    logger.info("Writing {}".format(path))
    with open(str(path), 'w') as f:
        f.write(" File written on {}\n".format(_datestamp()))
        f.write("\n")
        f.write(" calc_only_A  :  F\n")
        f.write("\n")

        # ── real_lattice ──────────────────────────────────────────────────
        # Reader: float(line[:12]), float(line[12:24]), float(line[24:36])
        # Format: 3 × %12.7f  → each field exactly 12 chars wide
        f.write("begin real_lattice\n")
        for i in range(3):
            f.write("{:12.7f}{:12.7f}{:12.7f}\n".format(
                rvec[i, 0], rvec[i, 1], rvec[i, 2]))
        f.write("end real_lattice\n")
        f.write("\n")

        # ── recip_lattice (without 2π, in Å⁻¹) ──────────────────────────
        # Same slice widths as real_lattice.
        f.write("begin recip_lattice\n")
        for i in range(3):
            f.write("{:12.7f}{:12.7f}{:12.7f}\n".format(
                kvec_file[i, 0], kvec_file[i, 1], kvec_file[i, 2]))
        f.write("end recip_lattice\n")
        f.write("\n")

        # ── kpoints ───────────────────────────────────────────────────────
        # Reader: float(line[:14]), float(line[14:28]), float(line[28:42])
        # Format: 3 × %14.8f  → each field exactly 14 chars wide
        f.write("begin kpoints\n")
        f.write("{:6d}\n".format(nkp))
        for k in kpoints:
            f.write("{:14.8f}{:14.8f}{:14.8f}\n".format(k[0], k[1], k[2]))
        f.write("end kpoints\n")
        f.write("\n")

        # ── projections ───────────────────────────────────────────────────
        # Reader (position line):
        #   float(line[:10]), float(line[10:21]), float(line[21:32])
        # Format: %10.5f %11.5f %11.5f  → fields of width 10, 11, 11
        # The reader then skips the next (auxiliary) line unconditionally.
        f.write("begin projections\n")
        f.write("{:6d}\n".format(nband))
        for i in range(nband):
            px, py, pz = plist[i, 0], plist[i, 1], plist[i, 2]
            # Position line — field widths must match reader slice indices.
            f.write("{:10.5f}{:11.5f}{:11.5f}{:8d}{:5d}{:5d}\n".format(
                px, py, pz, 0, 1, 1))
            # Auxiliary line (Euler angles + zona) — skipped by _readNnkp.
            f.write(
                "{:12.7f}{:12.7f}{:12.7f}{:12.7f}"
                "{:12.7f}{:12.7f}{:8.2f}\n".format(
                    0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0))
        f.write("end projections\n")
        f.write("\n")

        # ── nnkpts stub (not read by lwann) ──────────────────────────────
        f.write("begin nnkpts\n")
        f.write("{:6d}\n".format(0))
        f.write("end nnkpts\n")
        f.write("\n")

        # ── exclude_bands stub ────────────────────────────────────────────
        f.write("begin exclude_bands\n")
        f.write("{:6d}\n".format(0))
        f.write("end exclude_bands\n")
        f.write("\n")

    logger.info("  {:d} k-point(s), {:d} projection(s).".format(nkp, nband))


# ---------------------------------------------------------------------------
# .wout stub writer
# ---------------------------------------------------------------------------

def _write_wout(path: Path) -> None:
    """
    Write ``<stem>.wout``.

    ``Wannier90Calculation._readWout`` iterates lines looking for one that
    contains both ``'Length Unit'`` and either ``'Ang'`` (→ lengthscale 1.0)
    or anything else (→ lengthscale = bohr2angstrom).

    ``readData`` wraps the ``_readWout`` call in a bare ``except`` and
    defaults to ``lengthscale = 1.0`` on any failure, so the wout file is
    technically optional — but writing a proper stub avoids the log warning
    and makes the workflow cleaner.

    ``tb.rvec`` is in Ångström, so we always write ``Length Unit = Ang``.
    """
    logger.info("Writing {}".format(path))
    with open(str(path), 'w') as f:
        f.write(
            " ltb Wannier90 stub -- generated by linretrace ltb --write-hr\n"
        )
        f.write(" Generated on {}\n".format(_datestamp()))
        f.write("\n")
        # The one line that _readWout actually needs:
        f.write("                 Length Unit = Ang\n")
        f.write("\n")


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def write_wannier90(tb, stem: str) -> None:
    """
    Write a Wannier90-compatible file set from a :class:`TightBinding`
    instance so that ``lwann`` can read it directly.

    Parameters
    ----------
    tb : TightBinding
        A fully initialised instance — ``computeData`` must have been called
        so that ``tb.hr``, ``tb.rpoints``, ``tb.nrp``, ``tb.rvec``,
        ``tb.kvec``, ``tb.nkx/y/z``, ``tb.kpoints``, ``tb.orbitals``, and
        ``tb.irreducible`` are all available.
    stem : str
        Filename prefix (may include a directory path).  Three files are
        created:

        * ``<stem>_hr.dat``
        * ``<stem>.nnkp``
        * ``<stem>.wout``

    Raises
    ------
    OSError
        If any of the three output files cannot be opened for writing.
    """
    stem  = str(stem)
    fhr   = Path(stem + '_hr.dat')
    fnnkp = Path(stem + '.nnkp')
    fwout = Path(stem + '.wout')

    logger.info(
        "\nWriting Wannier90-compatible output (stem={!r})".format(stem)
    )
    _write_hr(tb, fhr)
    _write_nnkp(tb, fnnkp)
    _write_wout(fwout)
    logger.info(
        "Wannier90 export complete:\n  {}\n  {}\n  {}".format(
            fhr, fnnkp, fwout)
    )
