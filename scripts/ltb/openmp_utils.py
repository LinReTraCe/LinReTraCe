"""
scripts/ltb/openmp_utils.py
---------------------------
Shared helper for setting OpenMP / BLAS thread-count environment variables
BEFORE numpy (or any BLAS-linked library) is imported.

Usage pattern in an ltb entry-point
------------------------------------
At the very top of the file, before *any* scientific imports:

    from scripts.ltb.openmp_utils import preparse_openmp, add_openmp_argument

    # Act on --openmp immediately, before numpy is imported.
    _ncores = preparse_openmp()   # reads sys.argv, sets env vars, returns int

Then inside parse_args():

    add_openmp_argument(parser)   # registers the flag for --help / validation

Then inside main(), after normal arg-parsing:

    logger.info(openmp_info_line(_ncores))   # optional: log the chosen value
"""

from __future__ import annotations

import logging
import os
import sys

logger = logging.getLogger(__name__)

# Every variable that controls thread counts in common BLAS / OpenMP stacks.
_THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)


def set_openmp_threads(ncores: int) -> None:
    """Export thread-count variables to *ncores* in the current process environment."""
    value = str(ncores)
    for var in _THREAD_ENV_VARS:
        os.environ[var] = value


def preparse_openmp() -> int | None:
    """
    Scan *sys.argv* for ``--openmp [N]`` **without** a full argparse pass so
    that the environment variables are set before numpy is imported.

    Returns
    -------
    int | None
        The number of threads that was applied, or *None* if ``--openmp`` was
        not present on the command line (environment left unchanged).

    Accepted forms
    --------------
    ``--openmp``      →  use all logical CPUs (os.cpu_count())
    ``--openmp N``    →  use exactly N threads  (N must be a positive integer)
    """
    argv = sys.argv[1:]
    ncores: int | None = None

    for i, token in enumerate(argv):
        if token == "--openmp":
            # Check whether the next token is a positive integer.
            if i + 1 < len(argv) and argv[i + 1].isdigit() and int(argv[i + 1]) > 0:
                ncores = int(argv[i + 1])
            else:
                ncores = os.cpu_count() or 1
            break
        if token.startswith("--openmp="):
            rest = token.split("=", 1)[1]
            if rest.isdigit() and int(rest) > 0:
                ncores = int(rest)
            else:
                ncores = os.cpu_count() or 1
            break

    if ncores is not None:
        set_openmp_threads(ncores)

    return ncores


def add_openmp_argument(parser) -> None:  # type: ignore[no-untyped-def]
    """
    Register ``--openmp`` with *parser* so it appears in ``--help`` output and
    is accepted during the normal argparse pass.

    The flag is purely informational at parse time — the actual env-var export
    was already done by :func:`preparse_openmp` before numpy was imported.
    """
    parser.add_argument(
        "--openmp",
        metavar="N",
        nargs="?",          # --openmp alone is valid (means "all cores")
        const=None,         # placeholder; real value was set in preparse_openmp
        default=None,
        type=int,
        help=(
            "Number of OpenMP / BLAS threads (sets OMP_NUM_THREADS and related\n"
            "variables before numpy is loaded).\n"
            "Omit N to use all logical CPU cores (default: 1 thread, i.e. this\n"
            "flag must be given explicitly to enable multi-threading)."
        ),
    )


def openmp_info_line(ncores: int | None) -> str:
    """Return a short log string suitable for logger.info()."""
    if ncores is None:
        return "OpenMP: threading not requested (single-threaded BLAS)."
    return (
        f"OpenMP: OMP_NUM_THREADS={ncores} "
        f"(and OPENBLAS_NUM_THREADS, MKL_NUM_THREADS, VECLIB_MAXIMUM_THREADS, "
        f"NUMEXPR_NUM_THREADS)."
    )
