"""
==============================================================================
          WAVE-CURRENT INTERACTION NOMOGRAM GENERATOR (NOMOGEN EDITION)
==============================================================================
  PROGRAM:      Wave-Current Multiplier Calculator (Nomogen_Plots)
  DEPENDENCIES: nomogen.py (Must be in same directory), pynomo, pypdf
  
  DESCRIPTION:
  ----------------------------------------------------------------------------
  This script generates a professional 2-page engineering nomogram set on 
  A3 PAPER (29.7cm x 42.0cm) to solve the exact Wave-Current Interaction 
  dispersion relation.
  
  It uses the 'Nomogen' numerical relaxation engine to create curved scales 
  that perfectly match the non-linear equation provided.

  THEORETICAL BASIS (Exact Integer Fractions):
  ----------------------------------------------------------------------------
  The core dispersion relation is solved via a Multiplier (M) applied to the 
  Linear (Airy) wavelength.
  
    L_final = L_linear * Multiplier

  The Multiplier is calculated via a split-step formula:
  
    1. Pivot Calculation:
       Pivot = exp((345/509) * x5) + x6 / (x6 + (435/526))
       
    2. Correction Term (Z):
       Z = -((5/67) / x2) * tanh(x5 * x8 + (2117/961))
       
    3. Final Summation:
       Multiplier = Pivot + Z

  Where dimensionless parameters are:
    x2 = ln(H/d)             [Relative Height Log]
    x5 = U * T / L_linear    [Doppler Parameter]
    x6 = U / C0              [Velocity Ratio]
    x8 = ln(T * sqrt(g/d))   [Dimensionless Period Log]

  OUTPUTS:
  ----------------------------------------------------------------------------
  1. nomogen.txt:        A detailed ASCII engineering report of the calculation.
  2. nomogen_plots.pdf:  A merged 2-page PDF containing:
                         - Page 1: Pivot Calculation (Doppler & Velocity)
                         - Page 2: Multiplier Calculation (Pivot & Correction)

==============================================================================
"""

import sys
import math
import os

# ==============================================================================
#  1. SYSTEM CONFIGURATION & DEPENDENCY CHECK
# ==============================================================================
print(">>> Initializing Dependencies and Environment...")

# --- 1.1 Force PyX to use LaTeX Engine ---
# This is CRITICAL for rendering mathematical formulas on the nomogram.
# Without this, text generation may fail silently or produce artifacts.
try:
    from pyx import text
    text.set(text.LatexEngine)
    print("[SYSTEM] PyX LaTeX Engine initialized successfully.")
except ImportError:
    print("[WARNING] PyX not found or LaTeX engine configuration failed.")
    print("          Formulas may not render correctly in the output PDF.")

# --- 1.2 Import Numerical & Graphic Libraries ---
try:
    import numpy as np
    import scipy
    import numba
    from pynomo.nomographer import Nomographer
    # Import the local Nomogen engine (Assumes nomogen.py is fixed for A3)
    from nomogen import Nomogen
except ImportError as e:
    print("\n[CRITICAL ERROR] Missing Python Dependencies.")
    print(f"Details: {e}")
    print("Please install required libraries by running:")
    print("   pip install pynomo numpy scipy numba")
    sys.exit(1)

# --- 1.3 Check for PDF Merging Capability ---
TRY_MERGE = True
try:
    from pypdf import PdfWriter, PdfReader
except ImportError:
    TRY_MERGE = False
    print("[WARNING] 'pypdf' library not found.")
    print("   The script will generate separate PDF files (nomogram1.pdf, etc).")
    print("   To enable automatic merging, run: pip install pypdf")


# ==============================================================================
#  2. PHYSICS ENGINE: LINEAR THEORY & PARAMETERS
# ==============================================================================

def solve_linear_wavelength(T, d):
    """
    Iteratively solves the linear dispersion relation: L = (gT^2/2pi) * tanh(2pi*d/L)
    
    Args:
        T (float): Wave Period [s]
        d (float): Water Depth [m]
        
    Returns:
        float: Linear Wavelength L [m]
    """
    g = 9.80665
    pi = 3.14159265359
    
    # Deep water initial guess
    L0 = (g * T**2) / (2 * pi) 
    L = L0
    
    # Fixed-point iteration (stable for all depths)
    # Convergence usually within 5-10 iterations
    for _ in range(100):
        if L == 0: break
        L_new = L0 * math.tanh((2 * pi * d) / L)
        if abs(L_new - L) < 0.0001: 
            return L_new
        L = L_new
    return L

def check_constraints(H, d, L):
    """
    Checks physical realism using the Miche Breaking Criterion.
    
    Args:
        H (float): Wave Height [m]
        d (float): Water Depth [m]
        L (float): Wavelength [m]
        
    Returns:
        str: Status message ("PHYSICS OK" or Warning)
    """
    pi = 3.14159265359
    k = 2 * pi / L
    
    # Miche Breaking Limit (Theoretical maximum wave height)
    H_max = (0.142 * L) * math.tanh(k * d)
    
    if H > H_max:
        return f"WARNING: WAVE BREAKING (H={H:.2f} > Hmax={H_max:.2f})"
    
    # Steepness Check
    steepness = H / L
    if steepness > 0.142:
        return f"WARNING: TOO STEEP (H/L={steepness:.3f} > 0.142)"
        
    return "PHYSICS OK"

def calculate_physics_case(H, T, d, U):
    """
    Performs the full physics calculation for a specific test case.
    Computes all dimensionless parameters (x2, x5, etc.) and formula terms.
    
    Args:
        H, T, d, U: Input physical parameters.
        
    Returns:
        dict: A structured dictionary containing inputs, intermediate terms,
              and final results for reporting and plotting.
    """
    g = 9.80665
    pi = 3.14159265359
    
    # 1. Linear Theory Baseline
    L_lin = solve_linear_wavelength(T, d)
    C0 = (g * T) / (2 * pi)
    
    # 2. Dimensionless Parameters
    # Use max(val, 1e-9) to prevent log(0) errors during theoretical boundary checks
    x2 = math.log(max(H/d, 1e-9))
    x5 = (U * T) / L_lin
    x6 = U / C0
    dim_T = T * math.sqrt(g / d)
    x8 = math.log(max(dim_T, 1e-9))
    
    # 3. Kernel Constants (Exact Fractions as floats)
    # These constants are derived from the GEP Symbolic Regression model
    C1 = 345.0 / 509.0     # approx 0.6778
    C2 = 435.0 / 526.0     # approx 0.8270
    C3 = 5.0 / 67.0        # approx 0.0746
    C4 = 2117.0 / 961.0    # approx 2.2029
    
    # 4. Formula Terms Calculation
    
    # Term 1: Exponential Doppler
    t1 = math.exp(x5 * C1)
    
    # Term 2: Velocity Ratio (Singularity guard for x6 + C2 = 0)
    denom = x6 + C2
    if denom == 0: denom = 1e-9
    t2 = x6 / denom
    
    # Pivot (Intermediate Result for Page 1)
    pivot = t1 + t2
    
    # Term 3: Correction (Z)
    tanh_arg = (x5 * x8) + C4
    t3 = (-C3 / x2) * math.tanh(tanh_arg)
    
    # 5. Final Results
    multiplier = pivot + t3
    L_final = L_lin * multiplier
    status = check_constraints(H, d, L_lin)
    
    # 6. Return Data Structure
    return {
        'inputs': {'H':H, 'T':T, 'd':d, 'U':U},
        'physics': {'L_lin': L_lin, 'C0': C0, 'status':status},
        'dim': {'x2':x2, 'x5':x5, 'x6':x6, 'x8':x8},
        'terms': {'z':t3, 't1':t1, 't2':t2},
        'results': {'pivot':pivot, 'z':t3, 'multiplier':multiplier, 'L_final': L_final}
    }


# ==============================================================================
#  3. NOMOGEN EQUATIONS (GEOMETRY DEFINITIONS)
# ==============================================================================

def eq_pivot(x5, x6):
    """
    Page 1 Equation: Pivot Calculation.
    Pivot = exp(C1 * x5) + x6 / (x6 + C2)
    """
    C1 = 345.0 / 509.0
    C2 = 435.0 / 526.0
    return math.exp(C1 * x5) + x6 / (x6 + C2)

def eq_multiplier(pivot, z):
    """
    Page 2 Equation: Multiplier Calculation.
    Multiplier = Pivot + Z
    """
    return pivot + z

# ==============================================================================
#  4. GRAPHICS ENGINE (SMART BOUNDS & LAYOUT)
# ==============================================================================

def get_smart_bounds(val, typical_min, typical_max):
    """
    Calculates axis bounds that include the specific case value while 
    maintaining a reasonable minimum engineering range.
    Adds 15% padding for visual comfort.
    """
    vmin = min(val, typical_min)
    vmax = max(val, typical_max)
    
    span = vmax - vmin
    if span == 0: span = 0.1
    pad = span * 0.15 
    
    return vmin - pad, vmax + pad

def create_nomograms(d):
    """
    Main function to configure and generate the 3 nomogram pages.
    
    Args:
        d (dict): The physics case data returned by calculate_physics_case()
    """
    print("\n[GRAPHICS] Configuring Nomogen parameters for A3 Output...")

    # --- A3 Paper Dimensions (cm) ---
    PAPER_H = 42.0 
    PAPER_W = 29.7 
    
    # --- Visual Styles ---
    # Standard LaTeX formatting for axis labels
    FMT_INPUT  = r"{\bf\LARGE $%3.2f$}"
    FMT_RESULT = r"{\bf\huge $%3.2f$}" 
    
    # Titles
    TITLE_DOPPLER = r"\bf\Huge Doppler ($x_5$)"
    TITLE_VEL     = r"\bf\Huge Vel Ratio ($x_6$)"
    TITLE_PIVOT   = r"\bf\Huge PIVOT"
    TITLE_REF     = r"\bf\Huge PIVOT"
    TITLE_CORR    = r"\bf\Huge Correction ($Z$)"
    TITLE_MULT    = r"\bf\Huge MULTIPLIER (M)"

    # Common Configuration for all pages
    COMMON_OPTS = {
        'paper_height': PAPER_H,
        'paper_width': PAPER_W,
        'transformations': [('rotate', 0.01), ('scale paper',)],
        'npoints': 15, # High resolution Chebyshev nodes
        'muShape': 0,  # 0 = Standard optimization (no shape regularization)
    }
    
    # --- Dynamic Range Calculation ---
    # We ensure the charts cover the specific case provided in 'd'
    x5_min, x5_max = get_smart_bounds(d['dim']['x5'], -0.8, 0.8)
    x6_min, x6_max = get_smart_bounds(d['dim']['x6'], -0.5, 0.5)
    p_min, p_max   = get_smart_bounds(d['results']['pivot'], 0.3, 2.2)
    z_min, z_max   = get_smart_bounds(d['results']['z'], -0.4, 0.4)
    m_min, m_max   = get_smart_bounds(d['results']['multiplier'], 0.6, 2.6)

    # --------------------------------------------------------------------------
    # PAGE 1: PIVOT CALCULATION
    # --------------------------------------------------------------------------
    print("[GRAPHICS] Generating Page 1: Pivot Calculation...")
    
    # Define Text Blocks (LaTeX syntax)
    FORMULA_X = 2.0
    FORMULA_Y = 11.0
	
    text_formula = (
        r"\vbox{\hsize=20cm \parindent=0pt \huge" + "\n" 
        r"{\bf EXACT FORMULA SPECIFICATION}" + "\n"
        r"\vskip 10pt" + "\n"
        r"$L = L_{linear} \times M$" + "\n"
        r"\vskip 5pt" + "\n"
        # Line 1: Pivot terms (Fractional)
        r"$M = \exp({345 \over 509} \cdot x_5) + {x_6 \over x_6 + {435 \over 526}}$" + "\n"
        r"\vskip 8pt" + "\n"
        # Line 2: Correction term
        r"$\qquad - {{5 \over 67} \over x_2} \tanh(x_5 x_8 + {2117 \over 961})$" + "\n"
        r"\vskip 8pt" + "\n"
        r"{\bf Where:}" + "\n"
        r"$x_5 = U T / L_{lin}, \quad x_6 = U / C_0$" + "\n"
        r"\vskip 3pt" + "\n"
        r"$x_2 = \ln(H/d), \quad x_8 = \ln(T \sqrt{g/d})$" + "\n"
        r"}"
    )

    P1_TEXT_X = 15.0
    P1_TEXT_Y = 4.0
    P1_TEXT_W = 29.5 - P1_TEXT_X
	
    text_p1 = (
        r"\vbox{\hsize=" + str(P1_TEXT_W) + "cm \parindent=0pt \large" + "\n"
        r"{\bf\LARGE STEP 1: CALCULATE PIVOT TERM}" + "\n"
        r"\vskip 5pt" + "\n"
        r"SUB-FORMULA: $Pivot = \exp({345 \over 509} \cdot x_5) + x_6 / (x_6 + {435 \over 526})$" + "\n"
        r"\vskip 15pt" + "\n"
        r"\halign{#\hfil&\qquad #\hfil\cr" + "\n"
        r"\noalign{{\bf INPUT PARAMETERS:}\vskip 5pt}" + "\n"
        r"1. Doppler ($x_5$) $= U \cdot T / L(lin)$ \cr" + "\n"
        r"2. Vel Ratio ($x_6$) $= U / C_0$ \cr" + "\n"
        r"\noalign{\vskip 15pt{\bf INSTRUCTIONS:}\vskip 5pt}" + "\n"
        r"A. Draw line from $x_5$ (Left) to $x_6$ (Right). \cr" + "\n"
        r"B. Read PIVOT at Center intersection. \cr" + "\n"
        r"C. Write value down for Page 2. \cr}" + "\n"
        r"}"
    )

    # Define Axes
    axis_x5 = {
        'u_min': x5_min, 'u_max': x5_max,
        'title': TITLE_DOPPLER, 
        'scale_type': 'linear smart',
        'tick_levels': 3, 'tick_text_levels': 1, 'tick_side': 'left',
        'text_format': FMT_INPUT
    }
    
    axis_x6 = {
        'u_min': x6_min, 'u_max': x6_max,
        'title': TITLE_VEL,
        'scale_type': 'linear smart',
        'tick_levels': 3, 'tick_text_levels': 1, 'tick_side': 'right',
        'text_format': FMT_INPUT
    }
    
    axis_pivot_1 = {
        'u_min': p_min, 'u_max': p_max,
        'title': TITLE_PIVOT,
        'scale_type': 'linear smart',
        'tick_levels': 4, 'tick_text_levels': 1,
        'text_format': FMT_RESULT 
    }

    # Define Block (Type 9 for Curved Scales)
    block_1 = {
        'block_type': 'type_9',
        'f1_params': axis_x5,
        'f2_params': axis_pivot_1,
        'f3_params': axis_x6,
        'isopleth_values': [[d['dim']['x5'], 'x', d['dim']['x6']]], # Plot the solution line
    }

    params_1 = {
        **COMMON_OPTS,
        'filename': 'nomogram1.pdf',
        'block_params': [block_1],
        'extra_texts': [
            {'x':P1_TEXT_X, 'y':P1_TEXT_Y, 'width':P1_TEXT_W, 'color':'black', 'text':text_p1, 'fontsize':14},
            {'x':FORMULA_X, 'y':FORMULA_Y, 'width':20.0, 'color':'black', 'text':text_formula, 'fontsize':14}
        ]
    }
    
    # Run Generation for Page 1
    print("[GRAPHICS] Optimizing Page 1 Geometry...")
    Nomogen(eq_pivot, params_1)
    
    print("[GRAPHICS] Rendering Page 1 PDF...")
    Nomographer(params_1)


    # --------------------------------------------------------------------------
    # PAGE 2: FINAL MULTIPLIER
    # --------------------------------------------------------------------------
    print("[GRAPHICS] Generating Page 2: Multiplier Calculation...")

    P2_TEXT_X = 15.0
    P2_TEXT_Y = 2.5
    P2_TEXT_W = 29.5 - P2_TEXT_X

    text_p2 = (
        r"\vbox{\hsize=" + str(P2_TEXT_W) + "cm \parindent=0pt \large" + "\n"
        r"{\bf\LARGE STEP 2: CALCULATE FINAL MULTIPLIER}" + "\n"
        r"\vskip 5pt" + "\n"
        r"FORMULA: $Multiplier = Pivot + Correction(Z)$" + "\n\n"
        r"FORMULA: $Z = - ({5 \over 67} / x_2) \cdot \tanh(x_5 \cdot x_8 + {2117 \over 961})$" + "\n"
        r"\vskip 15pt" + "\n"
        r"\halign{#\hfil&\qquad #\hfil\cr" + "\n"
        r"\noalign{{\bf INPUT PARAMETERS:}\vskip 5pt}" + "\n"
        r"1. Rel Height ($x_2$) $= \ln(H/d)$ \cr" + "\n"
        r"2. Dim Period ($x_8$) $= \ln(T \cdot \sqrt{g/d})$ \cr" + "\n"
        r"\noalign{\vskip 15pt{\bf INSTRUCTIONS:}\vskip 5pt}" + "\n"
        r"A. Mark PIVOT (from Page 1) on Left Axis. \cr" + "\n"
        r"B. Calculate Z and mark on Right Axis. \cr" + "\n"
        r"C. Connect Pivot and Z to find MULTIPLIER. \cr}" + "\n"
        r"\vskip 10pt" + "\n"
        r"{\bf FINAL CALCULATION:} $L = L(lin) \cdot M$" + "\n"
        r"}"
    )

    axis_pivot_2 = {
        'u_min': p_min, 'u_max': p_max,
        'title': TITLE_REF,
        'scale_type': 'linear smart',
        'tick_levels': 3, 'tick_text_levels': 1, 'tick_side': 'left',
        'text_format': FMT_INPUT
    }

    axis_z = {
        'u_min': z_min, 'u_max': z_max,
        'title': TITLE_CORR,
        'scale_type': 'linear smart',
        'tick_levels': 3, 'tick_text_levels': 1, 'tick_side': 'right',
        'text_format': FMT_INPUT
    }

    axis_mult = {
        'u_min': m_min, 'u_max': m_max,
        'title': TITLE_MULT,
        'scale_type': 'linear smart',
        'tick_levels': 4, 'tick_text_levels': 1,
        'text_format': FMT_RESULT 
    }

    block_2 = {
        'block_type': 'type_9',
        'f1_params': axis_pivot_2,
        'f2_params': axis_mult,
        'f3_params': axis_z,
        'isopleth_values': [[d['results']['pivot'], 'x', d['results']['z']]],
    }

    params_2 = {
        **COMMON_OPTS,
        'filename': 'nomogram2.pdf',
        'block_params': [block_2],
        'extra_texts': [{'x':P2_TEXT_X, 'y':P2_TEXT_Y, 'width':P2_TEXT_W, 'color':'black', 'text':text_p2, 'fontsize':14}]
    }

    # Run Generation for Page 2
    print("[GRAPHICS] Optimizing Page 2 Geometry...")
    Nomogen(eq_multiplier, params_2)
    
    print("[GRAPHICS] Rendering Page 2 PDF...")
    Nomographer(params_2)

# ==============================================================================
#  5. REPORT GENERATION & PDF MERGING
# ==============================================================================

def generate_report_text(data):
    """
    Generates a formatted text report of the calculations.
    """
    d = data['dim']; r = data['results']; p = data['physics']; t = data['terms']
    i = data['inputs']
    
    lines = [
        "="*80,
        "             WAVE-CURRENT INTERACTION: CALCULATION MANUAL",
        "="*80,
        f"CASE INPUTS: H={i['H']}m, T={i['T']}s, d={i['d']}m, U={i['U']}m/s",
        "-" * 80,
        "PART A: PHYSICS PRE-CALCULATIONS",
        f"   Status Check:        {p['status']}",
        f"1. L(linear) [Depth d]: {p['L_lin']:.4f} m",
        f"   (Solved iteratively from dispersion relation)",
        f"2. C0 [Deep Water]:     {p['C0']:.4f} m/s",
        "",
        "PART B: DIMENSIONLESS PARAMETERS",
        f"1. Doppler Param (x5):  {d['x5']:.5f}",
        f"2. Velocity Ratio (x6): {d['x6']:.5f}",
        f"3. Rel. Height Log (x2):{d['x2']:.5f}",
        f"4. Dim. Period Log (x8):{d['x8']:.5f}",
        "",
        "PART C: FORMULA TERMS",
        f"   Term 1 (Doppler Exp):  {t['t1']:.5f}",
        f"   Term 2 (Vel Ratio):    {t['t2']:.5f}",
        f"   Term 3 (Z Correction): {t['z']:.5f}",
        "",
        "PART D: FINAL RESULTS",
        f"   Step 1 Pivot:          {r['pivot']:.5f}",
        f"   Final Multiplier:      {r['multiplier']:.5f}",
        f"   Final Wavelength:      {r['L_final']:.3f} m",
        "-" * 80,
        "GLOSSARY:",
        "   L(lin): Linear Airy wavelength at depth d.",
        "   C0:     Deep water phase speed.",
        "   Pivot:  Intermediate variable (Doppler + Velocity effects).",
        "   Z:      Correction term for wave steepness and depth.",
        "   M:      Final Multiplier (L_final / L_lin).",
        "="*80
    ]
    return "\n".join(lines)

def merge_pdfs_and_clean():
    """
    Merges the 3 generated PDF pages into a single file and deletes intermediates.
    """
    if not TRY_MERGE: 
        print("[INFO] Merge skipped (pypdf not installed).")
        return

    output_filename = "nomogen_plots.pdf"
    inputs = ["nomogram1.pdf", "nomogram2.pdf"]

    print(f"\n[MERGE] Attempting to combine: {inputs}")

    # Check which files actually exist
    found_files = [f for f in inputs if os.path.exists(f)]

    if found_files:
        try:
            writer = PdfWriter()
            for pdf_path in found_files:
                reader = PdfReader(pdf_path)
                for page in reader.pages:
                    writer.add_page(page)
            
            with open(output_filename, "wb") as f:
                writer.write(f)
            
            print(f"[SUCCESS] Merged {len(found_files)} files into '{output_filename}'")
            
            # Clean up intermediate files
            for pdf in found_files:
                try:
                    os.remove(pdf)
                except OSError as e:
                    print(f"[WARNING] Could not delete {pdf}: {e}")
            print("[CLEANUP] Deleted temporary single-page files.")
            
        except Exception as e:
            print(f"[ERROR] Failed to merge PDFs: {e}")
    else:
        print("[ERROR] No source PDF files found to merge.")


# ==============================================================================
#  6. MAIN EXECUTION BLOCK
# ==============================================================================

if __name__ == "__main__":
    print(">>> Generating Wave-Current Nomograms (Nomogen Edition)...")
    try:
        # 1. Calculation (Example Case)
        #    H=3.0m, T=9.0s, d=5.0m, U=1.0m/s
        case_data = calculate_physics_case(H=3.0, T=9.0, d=5.0, U=1.0)
        
        # 2. Generate and Save Report
        full_report = generate_report_text(case_data)
        print("\n" + full_report)
        
        with open("nomogen.txt", "w") as f:
            f.write(full_report)
        print("\n[INFO] Detailed report saved to 'nomogen.txt'")
        
        # 3. Generate Graphics (2 Pages)
        create_nomograms(case_data)
        
        # 4. Merge Output PDFs
        # We perform a check to see if files were actually generated before merging
        merge_pdfs_and_clean()
        
    except Exception as e:
        print(f"\n[CRITICAL EXECUTION ERROR] {e}")
        import traceback
        traceback.print_exc()