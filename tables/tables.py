"""
==============================================================================
TABLE GENERATOR USING EXTERNAL FENTON SOLVER
==============================================================================

DESCRIPTION:
  This script generates wavelength tables for combinations of:

    - Wave Height (H)
    - Wave Period (T)
    - Water Depth (d)
    - Current (Uc)
	
  Wavelengths are NOT computed internally in this script.
  Instead, each case is solved by the external module `fenton.py`, using:

      fenton.L(H, T, d, U)

  This guarantees that the table values come directly from the standalone
  Fenton wavelength solver used elsewhere in the project.

USAGE:
  Place `tables.py` and `fenton.py` in the same folder, then run:

      python tables.py

OUTPUT:
  A text file named `output.txt` containing the wavelength tables.
==============================================================================
"""

from __future__ import annotations

import math
import multiprocessing
import os
import shutil
import sys
import time
from typing import List, Optional, Tuple

# ------------------------------------------------------------------------------
# RUNTIME STABILITY
# ------------------------------------------------------------------------------
# Lock BLAS/OpenMP threads to 1 so multiprocessing can use the CPU efficiently.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

# Ensure this script can import the sibling module fenton.py.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR and SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

try:
    import fenton
except Exception as exc:  # pragma: no cover
    raise RuntimeError(
        "Could not import 'fenton.py'. Make sure 'tables.py' and 'fenton.py' "
        "are in the same folder."
    ) from exc


# ============================================================================== 
# USER CONFIGURATION
# ============================================================================== 

# 1. Currents [m/s]
#    Positive values follow wave direction, negative oppose it.
UC_VALUES = [0.0, 0.5, 1.0, -0.5, -1.0]

# 2. Wave Heights H [m]
H_RANGE = [1.0, 2.5, 5.0, 7.5, 10.0, 12.5]

# 3. Wave Periods T [s]
T_RANGE = list(range(1, 22, 2))

# 4. Water Depths d [m]
D_RANGE = list(range(5, 51, 5))

# Output filename
OUTPUT_FILE = "output.txt"


# ============================================================================== 
# CORE LOGIC
# ============================================================================== 

def _miche_breaking_limit(L: float, d: float) -> float:
    """Return Miche breaking limit H_break [m] for a given wavelength and depth."""
    if L <= 0.0 or d <= 0.0:
        return 0.0
    k = 2.0 * math.pi / L
    return 0.142 * L * math.tanh(k * d)


def solve_case(
    H: float,
    T: float,
    d: float,
    Uc: float,
    guess_vector: Optional[object] = None,
) -> Tuple[float, None, str]:
    """
    Solve a single case through external fenton.py.

    Parameters
    ----------
    H, T, d, Uc : float
        Wave height, period, depth, and current.
    guess_vector : ignored
        Kept only for backward compatibility with the previous tables.py API.

    Returns
    -------
    L : float
        Computed wavelength [m], or 0.0 if the solve failed.
    None : placeholder
        Kept for compatibility with existing table-generation logic.
    status : str
        Formatted wavelength, or BREAK / FAIL.
    """
    del guess_vector  # external fenton.py owns its own convergence strategy

    try:
        L = float(fenton.L(H, T, d, Uc))
    except Exception:
        return 0.0, None, "FAIL"

    if not math.isfinite(L) or L <= 0.0:
        return 0.0, None, "FAIL"

    breaking_limit = _miche_breaking_limit(L, d)
    if breaking_limit > 0.0 and H > breaking_limit:
        return L, None, "BREAK"

    return L, None, f"{L:.3f}"


# ============================================================================== 
# WORKER LOGIC
# ============================================================================== 

def process_depth_block(args: Tuple[float, float, List[float], List[float], int, int]):
    """
    Worker executed by the multiprocessing pool.
    Calculates a wavelength table block for a specific current and depth.
    """
    uc, d, h_range, t_range, block_id, total_blocks = args
    _ = (block_id, total_blocks)  # kept for API clarity / future use

    output_lines = []
    output_lines.append(f"\n  [ DEPTH d = {d} m ]")

    # Dynamic table formatting based on number of H columns.
    col_width_per_h = 15
    base_width = 12
    table_line_len = base_width + (len(h_range) * col_width_per_h)
    separator_line = "  " + "-" * (table_line_len - 2)

    output_lines.append(separator_line)

    row_label = "T \\ H"
    header = f"  {row_label:<8} |"
    for h in h_range:
        header += f" {f'H={h}m':^12} |"
    output_lines.append(header)
    output_lines.append(separator_line)

    for t in t_range:
        row_str = f"  {f'T={t}s':<8} |"
        for h in h_range:
            L, _, status = solve_case(h, t, d, uc)
            _ = L  # explicit: table prints status/value string only
            row_str += f" {status:^12} |"
        output_lines.append(row_str)

    output_lines.append(separator_line)
    return uc, d, "\n".join(output_lines)


# ============================================================================== 
# MAIN EXECUTION
# ============================================================================== 

def main() -> None:
    print("--- FENTON TABLE GENERATOR (EXTERNAL fenton.py) ---")
    print(f"Target Output : {OUTPUT_FILE}")

    tasks = []
    idx = 0
    total_tasks = len(UC_VALUES) * len(D_RANGE)

    for uc in UC_VALUES:
        for d in D_RANGE:
            idx += 1
            tasks.append((uc, d, H_RANGE, T_RANGE, idx, total_tasks))

    cpu_count = multiprocessing.cpu_count()
    print(f"System Cores  : {cpu_count}")
    print(f"Total Blocks  : {total_tasks}")
    print("Starting calculation pool...")
    print("-" * 60)

    t0 = time.time()
    results = []

    with multiprocessing.Pool(processes=cpu_count) as pool:
        for i, result in enumerate(pool.imap(process_depth_block, tasks)):
            uc_res, d_res, _ = result
            percent = ((i + 1) / total_tasks) * 100.0
            elapsed = time.time() - t0

            print(
                f"[{percent:5.1f}%] Computed Block {i+1}/{total_tasks}: "
                f"Uc={uc_res:<4.1f} d={d_res:<2} "
                f"(Elapsed: {elapsed:.1f}s)"
            )
            results.append(result)

    t1 = time.time()
    print("-" * 60)
    print(f"Calculation finished in {t1 - t0:.2f} seconds.")

    print("Writing to disk...")
    with open(os.path.join(SCRIPT_DIR, OUTPUT_FILE), "w", encoding="utf-8") as f:
        current_uc = None

        for uc, d, block_text in results:
            if current_uc != uc:
                current_uc = uc
                table_idx = UC_VALUES.index(uc) + 1
                f.write("=" * 90 + "\n")
                f.write(f"  TABLE {table_idx}: WAVELENGTH (L) [m]\n")
                f.write(f"  CURRENT (Uc) = {uc:.1f} m/s\n")
                f.write("  Legend: 'BREAK' = Wave Breaks, 'FAIL' = No Conv\n")
                f.write("=" * 90 + "\n")

            f.write(block_text)
            f.write("\n\n")

    print(f"Success. File '{OUTPUT_FILE}' is ready.")

    cache_dir = os.path.join(SCRIPT_DIR, "__pycache__")
    if os.path.exists(cache_dir):
        shutil.rmtree(cache_dir)
        print(f"Cleaned up '{cache_dir}' directory.")


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
