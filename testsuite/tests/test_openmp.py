#!/usr/bin/env python3
"""
testsuite/tests/test_openmp.py
------------------------------
Checks whether np.linalg.eigh benefits from multithreading via the
OpenBLAS / MKL thread pool, and warns if the setup is suboptimal.

Exit codes
----------
  0  Threading is effective (speedup >= SPEEDUP_WARN at any tested size), OR
     only 1 logical CPU is available (nothing to parallelise).
  1  Threading gave no meaningful speedup despite multiple CPUs being present.
     This is a WARNING, not a hard failure: the test suite continues, but the
     user is advised to check their BLAS installation.

Invoked by
----------
  make test-openmp
  python3 testsuite/tests/test_openmp.py          # standalone
  python3 testsuite/tests/test_openmp.py --quick  # fewer sizes, faster CI
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

# ---------------------------------------------------------------------------
# Thresholds
# ---------------------------------------------------------------------------
SPEEDUP_WARN  = 1.4   # below this on the largest tested size → warning
SIZES_FULL    = [100, 200, 400, 800]
SIZES_QUICK   = [100, 300]
BATCH         = 6     # eigh calls per timing sample
REPEATS       = 3     # best-of-N repeats


# ---------------------------------------------------------------------------
# Subprocess worker
# ---------------------------------------------------------------------------
# The worker is a self-contained re-invocation of this script.  Setting
# OMP_NUM_THREADS must happen before numpy is imported, so we cannot do it
# inside a function in the same process; a fresh subprocess is the only
# portable guarantee.

def _worker_main(n: int, nthreads: int, batch: int) -> None:
    """Run `batch` eigh calls on n×n matrices; print JSON timing to stdout."""
    import numpy as np

    rng = np.random.default_rng(0)

    def make_hermitian(n: int) -> np.ndarray:
        A = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
        return (A + A.conj().T) / 2.0

    matrices = [make_hermitian(n) for _ in range(batch)]

    # threadpoolctl gives a reliable runtime knob when available
    try:
        from threadpoolctl import threadpool_limits
        ctx = threadpool_limits(limits=nthreads)
    except ImportError:
        from contextlib import nullcontext
        ctx = nullcontext()  # type: ignore[assignment]

    with ctx:
        t0 = time.perf_counter()
        for M in matrices:
            np.linalg.eigh(M)
        elapsed = time.perf_counter() - t0

    print(json.dumps({"n": n, "threads": nthreads, "elapsed": elapsed}))


def _run_subprocess(n: int, nthreads: int, batch: int) -> float:
    """Spawn a clean Python process; return elapsed seconds."""
    env = os.environ.copy()
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env[var] = str(nthreads)

    cmd = [sys.executable, __file__,
           "--_worker",
           "--_worker_n",       str(n),
           "--_worker_threads", str(nthreads),
           "--_worker_batch",   str(batch)]
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if result.returncode != 0:
        raise RuntimeError(f"Worker process failed:\n{result.stderr}")
    return json.loads(result.stdout.strip())["elapsed"]


# ---------------------------------------------------------------------------
# BLAS / thread-pool detection
# ---------------------------------------------------------------------------

def _blas_info() -> str:
    """Return a short human-readable BLAS back-end description."""
    try:
        from threadpoolctl import threadpool_info
        pools = threadpool_info()
        if pools:
            parts = []
            for p in pools:
                parts.append(
                    f"{p.get('internal_api','?')} "
                    f"(prefix={p.get('prefix','?')}, "
                    f"max_threads={p.get('num_threads','?')})"
                )
            return "; ".join(parts)
    except ImportError:
        pass
    return "unknown (threadpoolctl not installed)"


# ---------------------------------------------------------------------------
# Main benchmark / check logic
# ---------------------------------------------------------------------------

def _run_check(sizes: list[int], ncores: int) -> int:
    """
    Benchmark eigh for each size in `sizes` using 1 vs `ncores` threads.
    Returns 0 (OK / warn-only) or 1 (actionable warning).
    """
    print()
    print("=" * 62)
    print("  LinReTraCe — OpenMP / BLAS threading check")
    print("=" * 62)

    import numpy as np
    print(f"  NumPy        : {np.__version__}")
    print(f"  Logical CPUs : {os.cpu_count() or 1}")
    print(f"  Testing      : 1 thread vs {ncores} threads")
    print(f"  BLAS back-end: {_blas_info()}")
    print(f"  Matrix sizes : {sizes}  (batch={BATCH}, best of {REPEATS})")
    print()

    results: list[dict] = []
    header = f"  {'n':>6}  {'t(1)':>8}  {'t(N)':>8}  {'speedup':>8}"
    print(header)
    print("  " + "-" * (len(header) - 2))

    for n in sizes:
        t1 = min(_run_subprocess(n, 1,      BATCH) for _ in range(REPEATS))
        tN = min(_run_subprocess(n, ncores, BATCH) for _ in range(REPEATS))
        speedup = t1 / tN if tN > 0 else float("nan")
        results.append({"n": n, "t1": t1, "tN": tN, "speedup": speedup})
        print(f"  {n:>6}  {t1:>8.3f}s  {tN:>8.3f}s  {speedup:>7.2f}x")

    print()

    max_speedup = max(r["speedup"] for r in results)
    best_n      = results[[r["speedup"] for r in results].index(max_speedup)]["n"]

    print("  Result")
    print("  " + "-" * 40)

    if ncores == 1:
        print("  Only 1 logical CPU detected — nothing to parallelise.")
        print("  OpenMP check skipped (not a problem).")
        print("=" * 62)
        print()
        return 0

    if max_speedup >= SPEEDUP_WARN:
        print(f"  PASS  Threading is effective.")
        print(f"        Best speedup {max_speedup:.2f}x at n={best_n}.")
        print(f"        --openmp will benefit ltb run / ltb refine.")
        print("=" * 62)
        print()
        return 0
    else:
        print(f"  WARNING  No meaningful threading speedup detected")
        print(f"           (best: {max_speedup:.2f}x at n={best_n}).")
        print()
        print("  Possible causes:")
        print("    • numpy is linked against a single-threaded BLAS")
        print("      (e.g. the Netlib reference BLAS, or a no-thread OpenBLAS build).")
        print("    • OMP_NUM_THREADS / OPENBLAS_NUM_THREADS are overridden")
        print("      by a site-wide environment setting.")
        print()
        print("  Diagnosis:")
        print("    python3 -c \"import numpy; numpy.show_config()\"")
        print("    pip install threadpoolctl")
        print("    python3 -c \"from threadpoolctl import threadpool_info; "
              "print(threadpool_info())\"")
        print()
        print("  Note: --openmp will have no effect until this is resolved.")
        print("  The rest of the LinReTraCe test suite is unaffected.")
        print("=" * 62)
        print()
        return 1


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("--quick", action="store_true",
                        help="Use a smaller set of matrix sizes (faster, for CI).")
    # Internal worker arguments — not for direct use.
    parser.add_argument("--_worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--_worker_n",       type=int, help=argparse.SUPPRESS)
    parser.add_argument("--_worker_threads", type=int, help=argparse.SUPPRESS)
    parser.add_argument("--_worker_batch",   type=int, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    if args._worker:
        _worker_main(args._worker_n, args._worker_threads, args._worker_batch)
        return 0

    ncores = os.cpu_count() or 1
    sizes  = SIZES_QUICK if args.quick else SIZES_FULL
    return _run_check(sizes, ncores)


if __name__ == "__main__":
    sys.exit(main())
