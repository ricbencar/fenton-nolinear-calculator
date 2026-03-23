"""
==============================================================================
FENTON WAVELENGTH LIST GENERATOR (driver)
==============================================================================

DESCRIPTION:
  Generates a tab-separated dataset of nonlinear wavelengths (L) by calling the
  standalone solver in `function.py` (Fenton stream-function / Fourier method).

VALIDITY RULES:
  A case is written to list.txt if, and only if:

    1) `function.py` returns a finite wavelength L in the range:
           1.0 <= L <= 1000.0   [meters]

    2) The wave is NON-BREAKING by Miche breaking criterion:
           H/L <= 0.142 * tanh( 2π d / L )

    3) For adverse current (Uc < 0), Doppler blocking does NOT occur
       (i.e., the intrinsic linear phase speed exceeds |Uc|):
           c0 + Uc > 0
       where c0 is the linear intrinsic phase speed computed from the
       (k, d) pair:
           c0 = sqrt(g/k * tanh(k d))

RANGES (grid):
  - Height (H):   0.0 to 15.0 m   (step 2.5)
  - Period (T):   0.0 to 20.0 s   (step 1.0)
  - Depth (d):    1.0 to 100.0 m  (step 2.5)
  - Current (Uc): -3.5 to 3.5 m/s (step 0.5)

OUTPUT:
  - list.txt: Tab-separated columns (H, T, d, Uc, L)

USAGE:
  python list.py
==============================================================================
"""

from __future__ import annotations

import math
import multiprocessing as mp
import os
import sys
import time
from dataclasses import dataclass
from typing import List, Tuple

# --- ENVIRONMENT CONFIGURATION ------------------------------------------------
# Avoid oversubscription: each process should be single-threaded in BLAS/LAPACK.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import numpy as np

# Import the solver module. This is the ONLY place wavelength is computed.
# Ensure the local directory (where list.py and function.py live) is on sys.path.
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR and (_THIS_DIR not in sys.path):
    sys.path.insert(0, _THIS_DIR)

import function as fenton  # provides fenton.L(H, T, d, U)

# --- CONFIGURATION ------------------------------------------------------------

OUTPUT_FILE = "list.txt"

L_MIN = 1.0
L_MAX = 1000.0

# Miche (breaking) criterion constant:
MICHE_C = 0.142

# Gravity (m/s^2) for Doppler / linear phase speed screening:
G = 9.81

H_RANGE = np.round(np.arange(0.0, 15.0 + 0.5, 2.5), 10).tolist()
T_RANGE = np.round(np.arange(0.0, 20.0 + 0.5, 1.0), 10).tolist()
D_RANGE = np.round(np.arange(1.0, 100.0 + 0.5, 2.5), 10)
UC_VALUES = np.round(np.arange(-3.5, 3.5 + 0.5, 0.5), 10)

# --- HELPERS -----------------------------------------------------------------


def _is_valid_L(value: object) -> Tuple[bool, float]:
    """Return (ok, L). ok is True iff value is finite and L_MIN <= L <= L_MAX."""
    try:
        L = float(value)  # accepts float, numpy scalars, numeric strings
    except Exception:
        return False, 0.0

    if not math.isfinite(L):
        return False, 0.0

    if (L < L_MIN) or (L > L_MAX):
        return False, L

    return True, L


def _passes_miche_nonbreaking(H: float, L: float, d: float) -> bool:
    """
    Miche (1944) breaking limit in finite depth:
        H/L <= 0.142 * tanh( 2π d / L )

    Returns True if non-breaking (passes), False if breaking.
    """
    if H <= 0.0:
        return True
    if L <= 0.0 or d <= 0.0:
        return False

    steep = H / L
    limit = MICHE_C * math.tanh(2.0 * math.pi * d / L)

    # Small numerical tolerance for edge cases:
    return steep <= (limit * (1.0 + 1e-12))


def _passes_doppler_blocking(U: float, L: float, d: float) -> bool:
    """
    Adverse-current Doppler blocking screen (only for U < 0):

      Compute intrinsic linear phase speed from (k,d):
          c0 = sqrt(g/k * tanh(k d))

      Discard if the opposing current cancels/blocks propagation:
          c0 + U <= 0

    Returns True if allowed, False if blocked/cancelled.
    """
    if U >= 0.0:
        return True
    if L <= 0.0 or d <= 0.0:
        return False

    k = 2.0 * math.pi / L
    if not math.isfinite(k) or k <= 0.0:
        return False

    c0 = math.sqrt((G / k) * math.tanh(k * d))
    # If opposing current magnitude reaches/exceeds intrinsic phase speed -> block
    return (c0 + U) > 0.0


# --- PARALLEL WORKER ----------------------------------------------------------


@dataclass
class BlockStats:
    tried: int = 0
    written: int = 0
    solver_fail: int = 0
    invalid_L: int = 0
    breaking: int = 0
    doppler_blocked: int = 0


def _process_block(args: Tuple[float, float]) -> Tuple[float, float, str, int, BlockStats]:
    """Worker: compute all (H,T) for a fixed (U,d) block."""
    U, d = args
    U = float(U)
    d = float(d)

    lines: List[str] = []
    st = BlockStats()

    for T in T_RANGE:
        T = float(T)
        for H in H_RANGE:
            H = float(H)
            st.tried += 1

            try:
                L_raw = fenton.L(H, T, d, U)
            except Exception:
                st.solver_fail += 1
                continue

            ok, L = _is_valid_L(L_raw)
            if not ok:
                st.invalid_L += 1
                continue

            # Miche breaking criterion (discard breaking waves):
            if not _passes_miche_nonbreaking(H, L, d):
                st.breaking += 1
                continue

            # Adverse-current Doppler cancellation/blocking (discard if blocked):
            if not _passes_doppler_blocking(U, L, d):
                st.doppler_blocked += 1
                continue

            lines.append(f"{H:.2f}\t{T:.2f}\t{d:.2f}\t{U:.2f}\t{L:.10f}\n")
            st.written += 1

    text = "".join(lines)
    return U, d, text, st.written, st


# --- MAIN ---------------------------------------------------------------------


def main() -> int:
    print("--- FENTON LIST GENERATOR (L-range + Miche non-breaking + Doppler screen) ---")
    print(f"Output file  : {OUTPUT_FILE}")

    if os.path.exists(OUTPUT_FILE):
        print(f"ERROR: '{OUTPUT_FILE}' already exists. Rename/delete it before running.")
        return 2

    print(
        f"Grid         : U({len(UC_VALUES)}) x d({len(D_RANGE)}) x T({len(T_RANGE)}) x H({len(H_RANGE)})"
    )
    print(f"Total cases  : {len(UC_VALUES) * len(D_RANGE) * len(T_RANGE) * len(H_RANGE)}")
    print(f"Blocks       : {len(UC_VALUES) * len(D_RANGE)} (one block per (U,d))")

    tasks = [(float(U), float(d)) for U in UC_VALUES for d in D_RANGE]

    cpu_count = mp.cpu_count()
    print(f"CPU cores    : {cpu_count}")
    print(f"Strategy     : {cpu_count} processes x 1 BLAS thread/process")
    print("-" * 72)

    t0 = time.time()

    total_written = 0
    agg = BlockStats()

    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        f.write("H\tT\td\tUc\tL\n")

        with mp.Pool(processes=cpu_count) as pool:
            for i, (U, d, text, nrows, st) in enumerate(
                pool.imap_unordered(_process_block, tasks), start=1
            ):
                agg.tried += st.tried
                agg.written += st.written
                agg.solver_fail += st.solver_fail
                agg.invalid_L += st.invalid_L
                agg.breaking += st.breaking
                agg.doppler_blocked += st.doppler_blocked

                if text:
                    f.write(text)
                total_written += nrows

                elapsed = time.time() - t0
                pct = 100.0 * i / len(tasks)
                rate = total_written / elapsed if elapsed > 0 else 0.0
                print(
                    f"[{pct:6.2f}%] block {i:4d}/{len(tasks)}  U={U:>5.2f}  d={d:>7.2f}  "
                    f"+{nrows:4d} rows  total={total_written:8d}  ({rate:6.1f} rows/s)"
                )

    dt = time.time() - t0
    print("-" * 72)
    print(f"Finished in  : {dt:.2f} s")
    print(f"Attempts     : {agg.tried}")
    print(f"Written      : {agg.written}")
    print(
        "Rejected     : "
        f"invalid_L={agg.invalid_L}, breaking={agg.breaking}, doppler_blocked={agg.doppler_blocked}, "
        f"solver_fail={agg.solver_fail}"
    )
    print(f"Done         : {OUTPUT_FILE}")
    return 0


if __name__ == "__main__":
    mp.freeze_support()
    raise SystemExit(main())
