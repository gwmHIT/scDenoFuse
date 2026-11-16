#!/usr/bin/env python3
"""Internal bridge; original algorithm and defaults stay unchanged."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'bin'))
import scFusion
if __name__ == '__main__':
    scFusion.main(sys.argv[1:], diagnostics=False)
