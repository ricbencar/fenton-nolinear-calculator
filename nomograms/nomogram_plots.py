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

# ------------------------------------------------------------------------------
# 1. DEPENDENCY VALIDATION
# ------------------------------------------------------------------------------
try:
    from pynomo.nomographer import Nomographer
except ImportError:
    print("\n[CRITICAL ERROR] 'pynomo' library is not installed.")
    print("Please install it by running: pip install pynomo\n")
    sys.exit(1)

TRY_MERGE = True
try:
    from PyPDF2 import PdfMerger
except ImportError:
    TRY_MERGE = False
    print("[INFO] PyPDF2 not found. Output will be generated as separate PDF files.")

# ------------------------------------------------------------------------------
# 2. PHYSICS ENGINE
# ------------------------------------------------------------------------------

def solve_linear_wavelength(T, d):
    """Iteratively solves the dispersion relation for linear waves."""
    g = 9.80665
    pi = 3.14159265359
    L0 = (g * T**2) / (2 * pi) 
    L = L0
    
    for _ in range(100):
        if L == 0: break
        L_new = L0 * math.tanh((2 * pi * d) / L)
        if abs(L_new - L) < 0.001: 
            return L_new
        L = L_new
    return L

def check_constraints(H, d, L):
    """Validates physical realism."""
    pi = 3.14159265359
    k = 2 * pi / L
    H_max = (0.142 * L) * math.tanh(k * d)
    
    if H > H_max:
        return f"WARNING: WAVE BREAKING (H={H:.2f} > Hmax={H_max:.2f})"
    
    steepness = H / L
    if steepness > 0.142:
        return f"WARNING: TOO STEEP (H/L={steepness:.3f} > 0.142)"
        
    return "PHYSICS OK"

def calculate_physics_case(H, T, d, U):
    """Calculates parameters using high-precision fractions for the kernel."""
    g = 9.80665
    pi = 3.14159265359
    
    L_lin = solve_linear_wavelength(T, d)
    C0 = (g * T) / (2 * pi)
    
    x2 = math.log(max(H/d, 1e-9))
    x5 = (U * T) / L_lin
    x6 = U / C0
    dim_T = T * math.sqrt(g / d)
    x8 = math.log(max(dim_T, 1e-9))
    
    # Kernel constants (keeping fractions for internal precision)
    C1 = 345.0 / 509.0
    C2 = 435.0 / 526.0
    C3 = 5.0 / 67.0
    C4 = 2117.0 / 961.0
    
    t1 = math.exp(x5 * C1)
    t2 = x6 / (x6 + C2)
    pivot = t1 + t2
    
    tanh_arg = (x5 * x8) + C4
    t3 = (-C3 / x2) * math.tanh(tanh_arg)
    
    multiplier = pivot + t3
    L_final = L_lin * multiplier
    status = check_constraints(H, d, L_lin)
    
    return {
        'inputs': {'H':H, 'T':T, 'd':d, 'U':U},
        'physics': {'L_lin': L_lin, 'C0': C0, 'status':status},
        'dim': {'x2':x2, 'x5':x5, 'x6':x6, 'x8':x8},
        'terms': {'z':t3, 't1':t1, 't2':t2},
        'results': {'pivot':pivot, 'multiplier':multiplier, 'L_final': L_final}
    }

# ------------------------------------------------------------------------------
# 3. REPORT GENERATOR
# ------------------------------------------------------------------------------

def generate_report_text(data):
    """Generates the engineering report."""
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

def get_smart_bounds(val, typical_min, typical_max):
    """Auto-calculates axis bounds."""
    vmin = min(val, typical_min)
    vmax = max(val, typical_max)
    span = vmax - vmin
    if span == 0: span = 0.1
    pad = span * 0.1
    return {'min': vmin - pad, 'max': vmax + pad}

# ------------------------------------------------------------------------------
# 4. NOMOGRAM GRAPHICS ENGINE
# ------------------------------------------------------------------------------

def create_nomograms(d):
    """Configures and renders the A3 nomograms."""
    
    # --- A. Transforms ---
    K1 = 345.0 / 509.0
    K2 = 435.0 / 526.0

    def func_t1(x): return math.exp(x * K1)
    def func_t2(x): return x / (x + K2)
    def func_neg(x): return -1.0 * x 
    def func_id(x): return x

    # --- B. Axis Bounds ---
    
    # 1. Input Axes (Fixed based on user requirement)
    # x5: [-0.8, 0.8] covers the operational Doppler range
    b_x5 = {'min': -0.8, 'max': 0.8}
    
    # x6: [-0.5, 0.5] clips the singularity at -0.827 for safety
    b_x6 = {'min': -0.5, 'max': 0.5}
    
    # Z: [-0.4, 0.4] covers typical correction magnitudes
    b_z  = {'min': -0.4, 'max': 0.4}

    # 2. Calculated Intermediate Axis (Pivot)
    # Max P (at 0.8, 0.5) approx 2.1
    # Min P (at -0.8, -0.5) approx -0.95
    # Bounds set to [-1.0, 2.2] to fully capture this.
    b_p  = {'min': -1.0, 'max': 2.2}

    # --- C. Visual Constants ---
    PAPER_H = 42.0 # A3 Height [cm]
    PAPER_W = 29.7 # A3 Width [cm]
    
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

    # --------------------------------------------------------------------------
    # PAGE 1: FORMULA DISPLAY & PIVOT
    # --------------------------------------------------------------------------
    
    # 1. Exact Formula Block (Updated to Decimal Coefficients & Large Font)
    FORMULA_X = 1.0
    FORMULA_Y = 36.0
    
    text_formula = (
        r"\vbox{\hsize=20cm \parindent=0pt \huge" + "\n" 
        r"{\bf EXACT FORMULA SPECIFICATION}" + "\n"
        r"\vskip 10pt" + "\n"
        r"$L = L_{linear} \times M$" + "\n"
        r"\vskip 5pt" + "\n"
        # Line 1: Pivot terms
        r"$M = \exp({345 \over 509} \cdot x_5) + {x_6 \over x_6 + {435 \over 526}}$" + "\n"
        r"\vskip 8pt" + "\n"
        # Line 2: Correction term (indented with \qquad)
        r"$\qquad - {{5 \over 67} \over x_2} \tanh(x_5 x_8 + {2117 \over 961})$" + "\n"
        r"\vskip 8pt" + "\n"
        r"{\bf Where:}" + "\n"
        r"$x_5 = U T / L_{lin}, \quad x_6 = U / C_0$" + "\n"
        r"\vskip 3pt" + "\n"
        r"$x_2 = \ln(H/d), \quad x_8 = \ln(T \sqrt{g/d})$" + "\n"
        r"}"
    )

    # 2. Instructions Block
    P1_TEXT_X = 18.0
    P1_TEXT_Y = 4.0
    P1_TEXT_W = 29.5 - P1_TEXT_X
    
    text_p1 = (
        r"\vbox{\hsize=" + str(P1_TEXT_W) + "cm \parindent=0pt \large" + "\n"
        r"{\bf\LARGE STEP 1: CALCULATE PIVOT TERM}" + "\n"
        r"\vskip 5pt" + "\n"
        # Replaced decimals 0.6778 and 0.8270
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

    block_1 = {
        'block_type': 'type_1',
        'width': 22.0, 
        'height': 16.0, 
        
        # --- LEFT AXIS: DOPPLER (x5) ---
        'f1_params': {
            'u_min': -0.8,      # Relaxed bound
            'u_max': 0.8,       # Relaxed bound
            'function': func_t1, 
            'title': TITLE_DOPPLER, 
            'tick_levels': 3, 
            'tick_text_levels': 1,
            'text_format': FMT_INPUT 
        },
        
        # --- MIDDLE AXIS: PIVOT (Result of Step 1) ---
        'f2_params': {
            'u_min': -1.0,      # Calculated to match input range
            'u_max': 2.2,       # Calculated to match input range
            'function': func_neg, 
            'title': TITLE_PIVOT,
            'tick_levels': 5,       # Adds many more subdivision marks
            'tick_text_levels': 3,  # Forces numbers to appear on minor ticks
            
            'text_format': FMT_RESULT 
        },
        
        # --- RIGHT AXIS: VELOCITY RATIO (x6) ---
        'f3_params': {
            'u_min': -0.5,      # Clipped to avoid singularity at -0.827
            'u_max': 0.5,       # Clipped for readability
            'function': func_t2, 
            'title': TITLE_VEL, 
            'tick_levels': 3, 
            'tick_text_levels': 1, 
            'text_format': FMT_INPUT
        },
        
        # --- ISOPLETH (The Red Line) ---
        'isopleth_values': [[d['dim']['x5'], 'x', d['dim']['x6']]],
        
        # --- PLACEMENT ON PAPER ---
        'transformations': [('translate', 3.0, 26.0)] 
    }

    params_1 = {
        'filename': 'nomogram1.pdf',
        'paper_height': PAPER_H, 'paper_width': PAPER_W,
        'block_params': [block_1],
        'transformations': [('rotate', 0.01), ('scale paper',)],
        # 'title_str': r'\bf\Huge Wave-Current Nomogram: Page 1',
        'extra_texts': [
            {'x':P1_TEXT_X, 'y':P1_TEXT_Y, 'width':P1_TEXT_W, 'color':'black', 'text':text_p1, 'fontsize':14},
            {'x':FORMULA_X, 'y':FORMULA_Y, 'width':20.0, 'color':'black', 'text':text_formula, 'fontsize':14}
        ]
    }

    # --------------------------------------------------------------------------
    # PAGE 2: MULTIPLIER CALCULATION
    # --------------------------------------------------------------------------

    P2_TEXT_X = 19.0
    P2_TEXT_Y = 1.5
    P2_TEXT_W = 32.75 - P2_TEXT_X
	
    text_p2 = (
        r"\vbox{\hsize=" + str(P2_TEXT_W) + "cm \parindent=0pt \large" + "\n"
        r"{\bf\LARGE STEP 2: CALCULATE FINAL MULTIPLIER}" + "\n"
        r"\vskip 5pt" + "\n"
        r"FORMULA: $Multiplier = Pivot + Correction(Z)$" + "\n\n"
        # Replaced decimals 0.0746 and 2.2029
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

    block_2 = {
        'block_type': 'type_1',
        'width': 22.0, 
        'height': 16.0,
        
        # --- LEFT AXIS: PIVOT (From Page 1) ---
        'f1_params': {
            'u_min': -1.0,      # Matches Page 1 output
            'u_max': 2.2,       # Matches Page 1 output
            'function': func_id, 
            'title': TITLE_REF, 
            
            # HIGH DENSITY TICKS (Kept high to match Page 1 precision)
            'tick_levels': 5,       
            'tick_text_levels': 3,
            
            'text_format': FMT_INPUT
        },
        
        # --- MIDDLE AXIS: MULTIPLIER (Final Result) ---
        'f2_params': {
            'u_min': -1.4,      # Covers blocking regime
            'u_max': 2.6,       # Covers following currents
            'function': func_neg, 
            'title': TITLE_MULT,
            'tick_levels': 5,       # Standard ticks
            'tick_text_levels': 3,  # Only major numbers
            'tick_side': 'right', 
            
            'text_format': FMT_RESULT
        },
        
        # --- RIGHT AXIS: CORRECTION (Z) ---
        'f3_params': {
            'u_min': -0.4,      # Zoomed in for typical steepness
            'u_max': 0.4,       # Zoomed in for typical steepness
            'function': func_id, 
            'title': TITLE_CORR, 
            
            # STANDARD TICKS
            'tick_levels': 3, 
            'tick_text_levels': 1, 
            'scale_type': 'linear', 
            
            'text_format': FMT_INPUT
        },
        
        # --- ISOPLETH (The Red Line) ---
        'isopleth_values': [[d['results']['pivot'], 'x', d['terms']['z']]],
        
        # --- PLACEMENT ON PAPER ---
        'transformations': [('translate', 3.0, 26.0)] 
    }

    params_2 = {
        'filename': 'nomogram2.pdf',
        'paper_height': PAPER_H, 'paper_width': PAPER_W,
        'block_params': [block_2],
        'transformations': [('rotate', 0.01), ('scale paper',)],
        # 'title_str': r'\bf\Huge Wave-Current Nomogram: Page 2',
        'extra_texts': [{'x':P2_TEXT_X, 'y':P2_TEXT_Y, 'width':P2_TEXT_W, 'color':'black', 'text':text_p2, 'fontsize':14}]
    }

    # --- Execution ---
    Nomographer(params_1)
    print("[GRAPHICS] Generated: nomogram1.pdf (A3 Format)")
    Nomographer(params_2)
    print("[GRAPHICS] Generated: nomogram2.pdf (A3 Format)")

# ------------------------------------------------------------------------------
# 5. MERGE UTILITIES
# ------------------------------------------------------------------------------

def merge_pdfs_and_clean():
    """Merges page 1 and 2 into nomogram_plots.pdf."""
    output_filename = "nomogram_plots.pdf"
    inputs = ["nomogram1.pdf", "nomogram2.pdf"]
    
    if TRY_MERGE:
        try:
            merger = PdfMerger()
            files_merged = 0
            for pdf in inputs:
                if os.path.exists(pdf):
                    merger.append(pdf)
                    files_merged += 1
            
            if files_merged > 0:
                merger.write(output_filename)
                merger.close()
                print(f"[SUCCESS] Merged nomograms into '{output_filename}'")
                for pdf in inputs:
                    if os.path.exists(pdf): os.remove(pdf)
                print("[CLEANUP] Deleted temporary single-page files.")
            else:
                print("[ERROR] No files were generated to merge.")
            
        except Exception as e:
            print(f"[MERGE ERROR] {e}")
    else:
        print("\n[INFO] PyPDF2 not installed. Merge skipped.")
        print("You have two separate files: 'nomogram1.pdf' and 'nomogram2.pdf'")

# ------------------------------------------------------------------------------
# 6. MAIN EXECUTION
# ------------------------------------------------------------------------------

if __name__ == "__main__":
    print(">>> Generating Wave-Current Interaction Nomograms (Decimal Formula)...")
    try:
        # 1. Calculation
        case_data = calculate_physics_case(H=3.0, T=9.0, d=5.0, U=1.0)
        
        # 2. Report
        full_report = generate_report_text(case_data)
        print(full_report)
        with open("nomogram.txt", "w") as f:
            f.write(full_report)
        print("\n[INFO] Detailed report saved to 'nomogram.txt'")
        
        # 3. Graphics
        create_nomograms(case_data)
        
        # 4. Merge
        merge_pdfs_and_clean()
        
    except Exception as e:
        print(f"\nCRITICAL EXECUTION ERROR: {e}")