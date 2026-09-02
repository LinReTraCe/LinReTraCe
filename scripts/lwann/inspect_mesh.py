#!/usr/bin/env python3
"""lwann inspect-mesh: thin wrapper around the generator-agnostic implementation.

The inspection needs nothing generator-specific -- energies, band-diagonal
moments, weights and the reciprocal vectors all come from the HDF5 -- so both
ltb and lwann share one implementation.
"""

import sys
from pathlib import Path

_root = Path(__file__).resolve().parents[2]
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))

from scripts.inspect_mesh import main  # noqa: E402,F401

if __name__ == "__main__":
    sys.exit(main())
