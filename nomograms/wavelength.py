"""
==============================================================================
  ENGINEERING TECHNICAL REFERENCE: LINEAR WAVE DISPERSION NOMOGRAM GENERATOR
==============================================================================
  PROGRAM:      wavelength.py
  DESCRIPTION:  Generates five high-precision nomograms for solving the Linear
                Wave Dispersion relation across shallow, intermediate and
                deeper-water depth ranges, then merges them into a single
                professional engineering PDF.
  THEORY:       Linear (Airy) Wave Theory
  METHOD:       Newton-Raphson Iterative Solver with PyNomo Geometry Generation
==============================================================================

  1. THEORETICAL FORMULATION (LINEAR AIRY THEORY)
  ----------------------------------------------------------------------------
  This software solves the implicit transcendental dispersion equation for
  water waves under the assumption of small-amplitude (linear) theory.

  A. Fundamental Equation (Dispersion Relation):
     The relationship between Wavelength (L), Period (T), and Water Depth (h)
     is given by:

     L = (g * T^2 / 2*pi) * tanh( 2*pi*h / L )

  B. Deep Water Asymptote (L0):
     In deep water (h/L > 0.5), the tanh term approaches 1.0, yielding the
     explicit deep-water wavelength:

     L0 = g * T^2 / 2*pi  (approx 1.56 * T^2)

  C. Numerical Solution Strategy:
     Since L appears on both sides of the dispersion equation, an iterative
     root-finding method is required.

     1. Initial Seed (Approximation):
        The script uses a weighted power-law approximation to estimate the
        initial tanh argument, ensuring rapid convergence:
        k0h = 2*pi*h / L0
        L/L0 = (6/5)^k0h * (k0h)^0.5

     2. Newton-Raphson Iteration:
        The solver refines the wavelength L by minimizing the residual:
        f(L) = L - L0 * tanh(2*pi*h / L)

        The iteration continues until the relative error (dL/L) < 1e-7.

==============================================================================

  2. OUTPUT CONFIGURATION
  ----------------------------------------------------------------------------
  This script generates five nomograms that together cover several water-depth
  regimes with overlap between charts to keep the graphical resolution useful:

     1. h =  1 to   5 m,  T =  5 to 20 s
     2. h =  5 to  13 m,  T = 10 to 20 s
     3. h = 13 to  25 m,  T =  5 to 20 s
     4. h = 25 to  50 m,  T =  7 to 20 s
     5. h = 50 to 100 m,  T = 10 to 20 s

  All five pages are merged into a single final file:
     wavelength.pdf

==============================================================================
"""

import math
import os
import sys
from pathlib import Path

import numpy as np
import scipy
from pypdf import PdfWriter

# --- FIX FOR SCIPY ATTRIBUTE ERROR ---
scipy.arange = np.arange
# ------------------------------------

# Ensure nomogen and pynomo are in the path
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(SCRIPT_DIR.parent))

from nomogen import Nomogen
from pynomo.nomographer import Nomographer
from pyx import text

# Ensure PyX uses the LaTeX engine
text.set(text.LatexEngine)


# =============================================================================
# SHARED WAVELENGTH CALCULATION FUNCTION
# =============================================================================
def LWave(h, T):
    """
    Calculates wavelength (L) from water depth (h) and period (T).
    """
    if T <= 0:
        return 0.0

    L0 = 9.80665 * T ** 2 / (2 * math.pi)
    k0h = 2 * math.pi * h / L0

    tanh_arg = (6 / 5) ** k0h * k0h ** 0.5
    L = L0 * math.tanh(tanh_arg)

    if L == 0:
        return 0.0

    dL = 1.0
    delta = 1e-7
    max_iter = 200
    count = 0

    while not abs(dL / L) < delta and count < max_iter:
        f1 = L - L0 * math.tanh(2 * math.pi * h / L)
        f2 = (L + delta) - L0 * math.tanh(2 * math.pi * h / (L + delta))

        if (f2 - f1) == 0:
            break

        dL = delta * f1 / (f2 - f1)
        L = L - dL
        count += 1

    return L


# =============================================================================
# CHART GENERATION HELPERS
# =============================================================================
def wavelength_limits(h_min, h_max, T_min, T_max):
    """
    Computes monotonic wavelength limits for the chart domain.
    """
    L_min = LWave(h_min, T_min)
    L_max = LWave(h_max, T_max)
    return L_min, L_max



def build_chart_params(spec, output_dir):
    """
    Builds a PyNomo/Nomogen parameter dictionary for one wavelength chart.
    """
    h_min, h_max = spec["h_range"]
    T_min, T_max = spec["T_range"]
    L_min, L_max = wavelength_limits(h_min, h_max, T_min, T_max)

    h_axis_params = {
        "u_min": h_min,
        "u_max": h_max,
        "title": r"Water Depth, h (m)",
        "tick_levels": 3,
        "tick_text_levels": spec.get("h_tick_text_levels", 2),
        "text_format": r"\Large{%g}",
    }
    T_axis_params = {
        "u_min": T_min,
        "u_max": T_max,
        "title": r"Wave Period, T (s)",
        "tick_levels": 3,
        "tick_text_levels": spec.get("T_tick_text_levels", 2),
        "text_format": r"\Large{%g}",
    }
    L_axis_params = {
        "u_min": L_min,
        "u_max": L_max,
        "title": r"Wavelength, L (m)",
        "tick_levels": 3,
        "tick_text_levels": spec.get("L_tick_text_levels", 2),
        "text_format": r"\Large\textbf{%g}",
    }

    block_params = {
        "block_type": "type_9",
        "f1_params": h_axis_params,
        "f2_params": L_axis_params,
        "f3_params": T_axis_params,
    }

    filename = output_dir / f"wavelength_{spec['index']}.pdf"
    title_str = (
        r"\Large Linear Wave Dispersion Nomogram: "
        rf"h = {h_min:g} to {h_max:g} m, "
        rf"T = {T_min:g} to {T_max:g} s"
    )

    return {
        "filename": str(filename),
        "paper_height": 29.7,
        "paper_width": 21,
        "block_params": [block_params],
        "transformations": [("scale paper",)],
        "title_str": title_str,
        "isopleth_params": [
            {"color": "Blue", "linewidth": "thick", "circle_size": 0.1}
        ],
        "npoints": spec.get("npoints", 12),
    }



def generate_chart(spec, output_dir):
    """
    Generates one chart PDF.
    """
    params = build_chart_params(spec, output_dir)
    filename = Path(params["filename"]).name
    print(f"\n--- Generating {filename} ---")
    print(
        "Calculating geometry for "
        f"{filename} (NN={params['npoints']}, "
        f"h={spec['h_range'][0]}-{spec['h_range'][1]} m, "
        f"T={spec['T_range'][0]}-{spec['T_range'][1]} s)..."
    )
    Nomogen(LWave, params)
    print(f"Rendering {filename}...")
    Nomographer(params)
    return Path(params["filename"])



def merge_pdfs(pdf_paths, output_filename):
    """
    Merges all generated pages into a single PDF file.
    """
    merger = PdfWriter()
    try:
        for pdf in pdf_paths:
            merger.append(str(pdf))
        merger.write(str(output_filename))
    finally:
        merger.close()


# =============================================================================
# MAIN EXECUTION
# =============================================================================
def main():
    chart_specs = [
        {
            "index": 1,
            "h_range": (1.0, 5.0),
            "T_range": (5.0, 20.0),
            "h_tick_text_levels": 3,
            "T_tick_text_levels": 2,
            "L_tick_text_levels": 2,
            "npoints": 12,
        },
        {
            "index": 2,
            "h_range": (5.0, 13.0),
            "T_range": (10.0, 20.0),
            "h_tick_text_levels": 2,
            "T_tick_text_levels": 2,
            "L_tick_text_levels": 2,
            "npoints": 12,
        },
        {
            "index": 3,
            "h_range": (13.0, 25.0),
            "T_range": (5.0, 20.0),
            "h_tick_text_levels": 2,
            "T_tick_text_levels": 2,
            "L_tick_text_levels": 2,
            "npoints": 12,
        },
        {
            "index": 4,
            "h_range": (25.0, 50.0),
            "T_range": (7.0, 20.0),
            "h_tick_text_levels": 2,
            "T_tick_text_levels": 2,
            "L_tick_text_levels": 2,
            "npoints": 12,
        },
        {
            "index": 5,
            "h_range": (50.0, 100.0),
            "T_range": (10.0, 20.0),
            "h_tick_text_levels": 2,
            "T_tick_text_levels": 2,
            "L_tick_text_levels": 2,
            "npoints": 12,
        },
    ]

    output_dir = SCRIPT_DIR
    generated_pdfs = []

    for spec in chart_specs:
        generated_pdfs.append(generate_chart(spec, output_dir))

    final_pdf = output_dir / "wavelength.pdf"
    print("\n--- Merging PDFs ---")
    merge_pdfs(generated_pdfs, final_pdf)
    print(f"Successfully created {final_pdf.name}")

    print("Cleaning up intermediate files...")
    for pdf in generated_pdfs:
        if pdf.exists():
            pdf.unlink()
            print(f"Deleted {pdf.name}")
        else:
            print(f"Warning: {pdf.name} not found for deletion.")

    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
