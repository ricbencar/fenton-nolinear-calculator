import pandas as pd
import numpy as np
from scipy.optimize import least_squares
import random
import time
import matplotlib.pyplot as plt
from numba import njit
import sympy as sp

# ==============================================================================
#  USER CONFIGURATION SECTION
# ==============================================================================
# Global final-prediction filter: after fitting a regime model, samples whose
# final predicted absolute percentage error exceeds this threshold are removed
# and the regime is refit once using only the retained samples.
ERROR_LIMIT_PCT = 10.0

# Adjust model complexity and physics parameters here.
# Supported Degrees: 2 to 7.
USER_CONFIG = {
    # Configuration for Shallow Water (d/L < 0.05)
    'shallow': {
        'degree': 3,      
        'lam': 1e-15      
    },
    
    # Configuration for Intermediate Water (0.05 <= d/L < 0.5)
    'inter': {
        'degree': 3,
        'lam': 1e-15
    },
    
    # Configuration for Deep Water (d/L >= 0.5)
    'deep': {
        'degree': 5,
        'lam': 1e-15
    }
}

# ==============================================================================
#  SECTION 1: PHYSICS KERNEL (JIT COMPILED)
# ==============================================================================
@njit(fastmath=True, cache=True)
def solve_linear_doppler(T_arr, d_arr, U_arr):
    """
    Computes the Linear Wavelength (L) using a Newton-Raphson solver.
    This serves as the 'Physics Baseline' (L_linear).
    """
    n = len(T_arr)
    L_res = np.zeros(n, dtype=np.float64)
    g = 9.80665
    pi = 3.14159265359
    
    for i in range(n):
        T = T_arr[i]
        d = d_arr[i]
        U = U_arr[i]
        
        # Initial Guess
        L0 = (g * T**2) / (2 * pi)
        omega = 2 * pi / T
        C_shallow = np.sqrt(g * d)
        
        # Wave Blocking Check
        if U < -0.5 * C_shallow:
            k = omega / (C_shallow * 0.6) 
        else:
            k = 2 * pi / L0
            
        # Newton-Raphson
        for _ in range(40):
            # Clamp k to prevent numerical overflow
            if k < 1e-4: k = 1e-4
            if k > 200.0: k = 200.0
            
            kd = k * d
            th = np.tanh(kd)
            sigma = np.sqrt(g * k * th)
            f = sigma + k * U - omega
            sech2 = 1.0 - th**2
            
            d_sigma = (g * th + g * kd * sech2) / (2 * sigma) if sigma > 1e-9 else 0.0
            df = d_sigma + U
            
            if abs(df) < 1e-9: break
            
            k_new = k - f / df
            k = 0.8 * k + 0.2 * k_new # Dampening
            
            if abs(k_new - k) < 1e-7: break
        
        if k <= 1e-4: k = 1e-4
        L_res[i] = 2 * pi / k
        
    return L_res

# ==============================================================================
#  SECTION 2: FEATURE ENGINEERING (JIT COMPILED)
# ==============================================================================
@njit(fastmath=True, cache=True)
def build_features(H, T, d, U, L_base):
    """
    Constructs the feature matrix X.
    Applies the CRITICAL L < 0.1 clamp to ensure stability.
    """
    n = len(H)
    X = np.zeros((n, 4), dtype=np.float64)
    
    for i in range(n):
        # [CRITICAL] Clamp L to 0.1 matching generation logic
        L = L_base[i]
        if L < 0.1: L = 0.1
        
        # 1. Wave Steepness (ln(H/L))
        steep = max(H[i] / L, 1e-7)
        X[i, 0] = np.log(steep)
        
        # 2. Relative Depth (ln(d/L))
        rel_d = max(d[i] / L, 1e-7)
        X[i, 1] = np.log(rel_d)
        
        # 3. Doppler Factor
        X[i, 2] = (U[i] * T[i]) / L
        
        # 4. Ursell Number (ln(Ur))
        d_safe = max(d[i], 0.1)
        ur = max((H[i] * L**2) / (d_safe**3), 1e-7)
        X[i, 3] = np.log(ur)
        
    return X

@njit(fastmath=True, cache=True)
def expand_poly_jit(X, n_samples, degree):
    """
    Polynomial Expansion supporting Degrees 2 through 7.
    Explicitly unrolled to ensure correct term order.
    """
    if degree == 7:   n_terms = 330
    elif degree == 6: n_terms = 210
    elif degree == 5: n_terms = 126
    elif degree == 4: n_terms = 70
    elif degree == 3: n_terms = 35
    else:             n_terms = 15
        
    Poly = np.zeros((n_samples, n_terms), dtype=np.float64)
    col = 0
    
    # 1. Constant
    Poly[:, col] = 1.0; col += 1
    
    # 2. Linear
    for i in range(4):
        Poly[:, col] = X[:, i]; col += 1
        
    # 3. Quadratic
    for i in range(4):
        for j in range(i, 4):
            Poly[:, col] = X[:, i] * X[:, j]; col += 1
            
    # 4. Cubic
    if degree >= 3:
        for i in range(4):
            for j in range(i, 4):
                for k in range(j, 4):
                    Poly[:, col] = X[:, i] * X[:, j] * X[:, k]; col += 1
                    
    # 5. Quartic
    if degree >= 4:
        for i in range(4):
            for j in range(i, 4):
                for k in range(j, 4):
                    for l in range(k, 4):
                        Poly[:, col] = X[:, i] * X[:, j] * X[:, k] * X[:, l]; col += 1
                        
    # 6. Quintic
    if degree >= 5:
        for i in range(4):
            for j in range(i, 4):
                for k in range(j, 4):
                    for l in range(k, 4):
                        for m in range(l, 4):
                            Poly[:, col] = X[:, i] * X[:, j] * X[:, k] * X[:, l] * X[:, m]; col += 1

    # 7. Sextic
    if degree >= 6:
        for i in range(4):
            for j in range(i, 4):
                for k in range(j, 4):
                    for l in range(k, 4):
                        for m in range(l, 4):
                            for n in range(m, 4):
                                Poly[:, col] = X[:, i] * X[:, j] * X[:, k] * X[:, l] * X[:, m] * X[:, n]; col += 1

    # 8. Septic
    if degree >= 7:
        for i in range(4):
            for j in range(i, 4):
                for k in range(j, 4):
                    for l in range(k, 4):
                        for m in range(l, 4):
                            for n in range(m, 4):
                                for o in range(n, 4):
                                    Poly[:, col] = X[:, i] * X[:, j] * X[:, k] * X[:, l] * X[:, m] * X[:, n] * X[:, o]; col += 1
    return Poly

# ==============================================================================
#  SECTION 3: MODEL & LOSS
# ==============================================================================
@njit(fastmath=True, cache=True)
def eval_rational(params, X_num, X_den):
    n_terms = X_num.shape[1]
    p = params[:n_terms]
    q = params[n_terms:]
    num = X_num @ p
    den = 1.0 + (X_den @ q) 
    return num, den

@njit(fastmath=True, cache=True)
def residuals_jit(params, X_num, X_den, target_log, weights):
    num, den = eval_rational(params, X_num, X_den)
    
    # Limitless Guard
    for i in range(len(den)):
        if abs(den[i]) < 0.001:
            den[i] = 0.001 * np.sign(den[i])
            if den[i] == 0: den[i] = 0.001
    
    log_pred = num / den
    
    # Physics Clamp matching function.py logic
    for i in range(len(log_pred)):
        if log_pred[i] > 2.5: log_pred[i] = 2.5
        if log_pred[i] < -2.5: log_pred[i] = -2.5
        
    return (log_pred - target_log) * weights

@njit(fastmath=True, cache=True)
def jacobian_jit(params, X_num, X_den, target, weights):
    num, den = eval_rational(params, X_num, X_den)
    for i in range(len(den)):
        if abs(den[i]) < 0.001:
            den[i] = 0.001 * np.sign(den[i])
            if den[i] == 0: den[i] = 0.001
            
    inv_den = 1.0 / den
    pred = num * inv_den
    
    n_samples = X_num.shape[0]
    n_p = X_num.shape[1]
    n_q = X_den.shape[1]
    
    Jac = np.zeros((n_samples, n_p + n_q), dtype=np.float64)
    w_col = np.ascontiguousarray(weights).reshape(-1, 1)
    
    # Derivative P
    term_p = np.ascontiguousarray(inv_den).reshape(-1, 1) * w_col
    Jac[:, :n_p] = X_num * term_p
    
    # Derivative Q
    term_q = -pred * inv_den
    term_q = np.ascontiguousarray(term_q).reshape(-1, 1) * w_col
    Jac[:, n_p:] = X_den * term_q
    
    return Jac

# ==============================================================================
#  SECTION 4: CODE GENERATION (DIRECT EXPANSION)
# ==============================================================================
def get_term_count(degree):
    if degree == 7: return 330
    if degree == 6: return 210
    if degree == 5: return 126
    if degree == 4: return 70
    if degree == 3: return 35
    return 15

def get_term_names(degree):
    # Generates variable strings matching expand_poly_jit order
    v = ['x0', 'x1', 'x2', 'x3']
    names = ["1.0"] 
    
    # Linear
    for i in range(4): names.append(v[i])
    # Quadratic
    for i in range(4):
        for j in range(i, 4): names.append(f"{v[i]}*{v[j]}")
    # Cubic
    if degree >= 3:
        for i in range(4):
            for j in range(i, 4):
                for k in range(j, 4): names.append(f"{v[i]}*{v[j]}*{v[k]}")
    # Quartic
    if degree >= 4:
        for i in range(4):
            for j in range(i, 4):
                for k in range(j, 4):
                    for l in range(k, 4): names.append(f"{v[i]}*{v[j]}*{v[k]}*{v[l]}")
    # Quintic
    if degree >= 5:
        for i in range(4):
            for j in range(i, 4):
                for k in range(j, 4):
                    for l in range(k, 4):
                        for m in range(l, 4): names.append(f"{v[i]}*{v[j]}*{v[k]}*{v[l]}*{v[m]}")
    # Sextic
    if degree >= 6:
        for i in range(4):
            for j in range(i, 4):
                for k in range(j, 4):
                    for l in range(k, 4):
                        for m in range(l, 4):
                            for n in range(m, 4): names.append(f"{v[i]}*{v[j]}*{v[k]}*{v[l]}*{v[m]}*{v[n]}")
    # Septic
    if degree >= 7:
        for i in range(4):
            for j in range(i, 4):
                for k in range(j, 4):
                    for l in range(k, 4):
                        for m in range(l, 4):
                            for n in range(m, 4):
                                for o in range(n, 4): names.append(f"{v[i]}*{v[j]}*{v[k]}*{v[l]}*{v[m]}*{v[n]}*{v[o]}")
    return names

def get_exact_python_string(coeffs, names):
    # Generates Sum-of-Products string with 16-digit precision
    # This guarantees 1:1 match with matrix multiplication logic
    parts = []
    for c, n in zip(coeffs, names):
        if abs(c) < 1e-20: continue # Skip near-zero terms
        if n == "1.0":
            parts.append(f"{c:.16e}")
        else:
            parts.append(f"{c:.16e}*{n}")
    
    if not parts: return "0.0"
    return " + ".join(parts)

def generate_multi_model_py(results_dict, filename="function.py"):
    
    # Extract configs
    deg_s = results_dict['shallow']['degree']
    deg_i = results_dict['inter']['degree']
    deg_d = results_dict['deep']['degree']

    # Helper to get P/Q strings using exact expansion
    def get_pq_strs(res, deg):
        n_p = get_term_count(deg)
        p_coeffs = res['params'][:n_p]
        q_coeffs = res['params'][n_p:]
        
        all_names = get_term_names(deg)
        # P uses all terms
        str_p = get_exact_python_string(p_coeffs, all_names)
        
        # Q uses terms starting from index 1, constant is 1.0 + ...
        # q_coeffs align with names[1:] (denominator doesn't fit bias)
        str_q = get_exact_python_string(q_coeffs, all_names[1:])
        str_q = f"1.0 + {str_q}" 
        return str_p, str_q

    p_s, q_s = get_pq_strs(results_dict['shallow'], deg_s)
    p_i, q_i = get_pq_strs(results_dict['inter'], deg_i)
    p_d, q_d = get_pq_strs(results_dict['deep'], deg_d)
    
    m_s, s_s = results_dict['shallow']['mean'], results_dict['shallow']['std']
    m_i, s_i = results_dict['inter']['mean'], results_dict['inter']['std']
    m_d, s_d = results_dict['deep']['mean'], results_dict['deep']['std']
    
    code = f'''# -*- coding: utf-8 -*-
"""
WAVELENGTH CALCULATOR (Non-Linear + Doppler)
============================================
Generated by pade.py
Precision: Float64 (16 digits)
Logic: Direct Expanded Polynomial
"""
import numpy as np

def L(H, T, d, Uc):
    # 1. Physics Baseline
    g = 9.80665
    pi = 3.14159265359
    omega = 2 * pi / T
    C_shallow = np.sqrt(g * d)
    
    if Uc < -0.5 * C_shallow:
        k = omega / (C_shallow * 0.6)
    else:
        L0 = (g * T**2) / (2 * pi)
        k = 2 * pi / L0
        
    for _ in range(40):
        k = max(1e-4, min(k, 200.0))
        kd = k * d
        th = np.tanh(kd)
        sigma = np.sqrt(g * k * th)
        f = sigma + k * Uc - omega
        sech2 = 1.0 - th**2
        d_sigma = (g * th + g * kd * sech2) / (2 * sigma) if sigma > 1e-9 else 0.0
        df = d_sigma + Uc
        if abs(df) < 1e-9: break
        k_new = k - f / df
        k = 0.8 * k + 0.2 * k_new 
        if abs(k_new - k) < 1e-7: break
        
    L_lin = 2 * pi / k
    
    # 2. Feature Extraction (Safe Baseline)
    L_feat = L_lin
    if L_feat < 0.1: L_feat = 0.1

    x0_raw = np.log(max(H/L_feat, 1e-7))
    x1_raw = np.log(max(d/L_feat, 1e-7))
    x2_raw = (Uc * T) / L_feat
    d_safe = max(d, 0.1)
    x3_raw = np.log(max((H * L_feat**2) / (d_safe**3), 1e-7))
    
    rel_depth = d / L_lin
    
    # 3. Model Selection
    if rel_depth < 0.05:
        m0, m1, m2, m3 = {m_s[0]:.16e}, {m_s[1]:.16e}, {m_s[2]:.16e}, {m_s[3]:.16e}
        s0, s1, s2, s3 = {s_s[0]:.16e}, {s_s[1]:.16e}, {s_s[2]:.16e}, {s_s[3]:.16e}
        x0, x1, x2, x3 = (x0_raw-m0)/s0, (x1_raw-m1)/s1, (x2_raw-m2)/s2, (x3_raw-m3)/s3
        
        P = {p_s}
        Q = {q_s}

    elif rel_depth < 0.5:
        m0, m1, m2, m3 = {m_i[0]:.16e}, {m_i[1]:.16e}, {m_i[2]:.16e}, {m_i[3]:.16e}
        s0, s1, s2, s3 = {s_i[0]:.16e}, {s_i[1]:.16e}, {s_i[2]:.16e}, {s_i[3]:.16e}
        x0, x1, x2, x3 = (x0_raw-m0)/s0, (x1_raw-m1)/s1, (x2_raw-m2)/s2, (x3_raw-m3)/s3
        
        P = {p_i}
        Q = {q_i}

    else:
        m0, m1, m2, m3 = {m_d[0]:.16e}, {m_d[1]:.16e}, {m_d[2]:.16e}, {m_d[3]:.16e}
        s0, s1, s2, s3 = {s_d[0]:.16e}, {s_d[1]:.16e}, {s_d[2]:.16e}, {s_d[3]:.16e}
        x0, x1, x2, x3 = (x0_raw-m0)/s0, (x1_raw-m1)/s1, (x2_raw-m2)/s2, (x3_raw-m3)/s3
        
        P = {p_d}
        Q = {q_d}
    
    # 4. Final Calculation
    if abs(Q) < 0.001: Q = 0.001 * np.sign(Q)
    if Q == 0: Q = 0.001
    
    log_corr = P / Q
    log_corr = max(-2.5, min(log_corr, 2.5))
    
    # Return L_feat (Safe Baseline) * Correction
    return L_feat * np.exp(log_corr)

if __name__ == "__main__":
    print(f"Test Case [H=3, T=9, d=5, U=1]: L = {{L(3, 9, 5, 1):.4f}} m")
'''
    with open(filename, "w", encoding="utf-8") as f: f.write(code)
    print(f"\n[SUCCESS] Independent function saved to '{filename}'")



def _term_name_to_vba(name):
    if name == "1.0":
        return "1#"
    counts = [0, 0, 0, 0]
    for part in name.split("*"):
        if not part.startswith("x"):
            continue
        idx = int(part[1])
        counts[idx] += 1
    pieces = [f"p{i}({power})" for i, power in enumerate(counts) if power > 0]
    return " * ".join(pieces) if pieces else "1#"


def _vba_poly_lines(coeffs, names, indent="        ", max_len=100, start_at_one=False):
    lines = [f"{indent}poly = {'1#' if start_at_one else '0#'}"]
    for coeff, name in zip(coeffs, names):
        if abs(coeff) < 1e-20:
            continue
        coeff_str = f"{coeff:.16E}"
        term_str = _term_name_to_vba(name)
        if term_str == "1#":
            lines.append(f"{indent}poly = poly + {coeff_str}")
            continue
        expr = f"{coeff_str} * {term_str}"
        line = f"{indent}poly = poly + {expr}"
        if len(line) <= max_len:
            lines.append(line)
        else:
            lines.append(f"{indent}poly = poly + {coeff_str} * _")
            lines.append(f"{indent}    {term_str}")
    return lines


def generate_multi_model_vba(results_dict, filename="pade.bas"):
    max_degree = max(
        results_dict['shallow']['degree'],
        results_dict['inter']['degree'],
        results_dict['deep']['degree'],
    )

    def build_case_block(case_name, reg_label):
        res = results_dict[case_name]
        degree = res['degree']
        n_p = get_term_count(degree)
        names = get_term_names(degree)
        p_coeffs = res['params'][:n_p]
        q_coeffs = res['params'][n_p:]
        means = res['mean']
        stds = res['std']

        block = []
        block.append(f"    Case {reg_label}")
        for i, value in enumerate(means):
            block.append(f"        m{i} = {value:.16E}")
        for i, value in enumerate(stds):
            block.append(f"        s{i} = {value:.16E}")
        block.append("")
        block.append("        x0 = (x0_raw - m0) / s0")
        block.append("        x1 = (x1_raw - m1) / s1")
        block.append("        x2 = (x2_raw - m2) / s2")
        block.append("        x3 = (x3_raw - m3) / s3")
        block.append("")
        block.append(
            f"        Call ComputePowers({degree}, x0, x1, x2, x3, p0, p1, p2, p3)"
        )
        block.append("")
        block.extend(_vba_poly_lines(p_coeffs, names))
        block.append("        P = poly")
        block.append("")
        block.extend(_vba_poly_lines(q_coeffs, names[1:], start_at_one=True))
        block.append("        Q = poly")
        block.append("")
        return "\n".join(block)

    case_blocks = [
        build_case_block('shallow', 'REG_SHALLOW'),
        build_case_block('inter', 'REG_INTER'),
        build_case_block('deep', 'REG_DEEP'),
    ]
    case_text = "\n".join(case_blocks)

    code = f'''Attribute VB_Name = "pade"
Option Explicit

Private Const G As Double = 9.80665#
Private Const PI As Double = 3.1415926535897931#
Private Const REL_SHALLOW As Double = 0.05#
Private Const REL_DEEP As Double = 0.5#
Private Const EPS_K As Double = 0.0001#
Private Const EPS_LOG As Double = 0.0000001#
Private Const EPS_DEN As Double = 0.001#
Private Const LOG_CORR_MIN As Double = -2.5#
Private Const LOG_CORR_MAX As Double = 2.5#

Private Const REG_SHALLOW As Long = 0
Private Const REG_INTER As Long = 1
Private Const REG_DEEP As Long = 2

Public Function L_pade(ByVal H As Double, ByVal T As Double, _
                       ByVal d As Double, ByVal Uc As Double) As Double
    Dim L_lin As Double
    Dim L_base As Double
    Dim rel_depth As Double
    Dim reg As Long

    If T <= 0# Or d <= 0# Then
        L_pade = 0#
        Exit Function
    End If

    L_lin = SolveLinearDoppler(T, d, Uc)
    L_base = MaxVal(L_lin, 0.1#)
    rel_depth = d / MaxVal(L_lin, EPS_K)

    If rel_depth < REL_SHALLOW Then
        reg = REG_SHALLOW
    ElseIf rel_depth < REL_DEEP Then
        reg = REG_INTER
    Else
        reg = REG_DEEP
    End If

    L_pade = PredictRegime(H, T, d, Uc, L_base, reg)
End Function

Public Function L(ByVal H As Double, ByVal T As Double, _
                  ByVal d As Double, ByVal Uc As Double) As Double
    L = L_pade(H, T, d, Uc)
End Function

Private Function SolveLinearDoppler(ByVal T As Double, ByVal d As Double, _
                                    ByVal Uc As Double) As Double
    Dim L0 As Double
    Dim omega As Double
    Dim C_shallow As Double
    Dim k As Double
    Dim kd As Double
    Dim th As Double
    Dim sigma As Double
    Dim f As Double
    Dim sech2 As Double
    Dim d_sigma As Double
    Dim df As Double
    Dim k_new As Double
    Dim diff_val As Double
    Dim it As Long

    omega = 2# * PI / T
    C_shallow = Sqr(G * d)

    If Uc < -0.5# * C_shallow Then
        k = omega / (C_shallow * 0.6#)
    Else
        L0 = (G * T * T) / (2# * PI)
        k = 2# * PI / L0
    End If

    For it = 1 To 40
        If k < EPS_K Then k = EPS_K
        If k > 200# Then k = 200#

        kd = k * d
        th = FnTanh(kd)
        sigma = Sqr(G * k * th)
        f = sigma + k * Uc - omega
        sech2 = 1# - th * th

        If sigma > 0.000000001# Then
            d_sigma = (G * th + G * kd * sech2) / (2# * sigma)
        Else
            d_sigma = 0#
        End If

        df = d_sigma + Uc
        If Abs(df) < 0.000000001# Then Exit For

        k_new = k - f / df
        diff_val = Abs(k_new - k)
        k = 0.8# * k + 0.2# * k_new

        If diff_val < 0.0000001# Then Exit For
    Next it

    If k <= EPS_K Then k = EPS_K
    SolveLinearDoppler = 2# * PI / k
End Function

Private Sub Features(ByVal H As Double, ByVal T As Double, ByVal d As Double, _
                     ByVal Uc As Double, ByVal L_base As Double, _
                     ByRef x0_raw As Double, ByRef x1_raw As Double, _
                     ByRef x2_raw As Double, ByRef x3_raw As Double)
    Dim L_feat As Double
    Dim d_safe As Double

    L_feat = MaxVal(L_base, 0.1#)
    x0_raw = Log(MaxVal(H / L_feat, EPS_LOG))
    x1_raw = Log(MaxVal(d / L_feat, EPS_LOG))
    x2_raw = (Uc * T) / L_feat
    d_safe = MaxVal(d, 0.1#)
    x3_raw = Log(MaxVal((H * L_feat * L_feat) / _
                        (d_safe * d_safe * d_safe), EPS_LOG))
End Sub

Private Sub ComputePowers(ByVal max_degree As Long, ByVal x0 As Double, _
                          ByVal x1 As Double, ByVal x2 As Double, _
                          ByVal x3 As Double, ByRef p0() As Double, _
                          ByRef p1() As Double, ByRef p2() As Double, _
                          ByRef p3() As Double)
    Dim ip As Long

    p0(0) = 1#: p1(0) = 1#: p2(0) = 1#: p3(0) = 1#
    For ip = 1 To max_degree
        p0(ip) = p0(ip - 1) * x0
        p1(ip) = p1(ip - 1) * x1
        p2(ip) = p2(ip - 1) * x2
        p3(ip) = p3(ip - 1) * x3
    Next ip
End Sub

Private Function PredictRegime(ByVal H As Double, ByVal T As Double, _
                               ByVal d As Double, ByVal Uc As Double, _
                               ByVal L_base As Double, ByVal reg As Long) As Double
    Dim x0_raw As Double, x1_raw As Double
    Dim x2_raw As Double, x3_raw As Double
    Dim x0 As Double, x1 As Double, x2 As Double, x3 As Double
    Dim m0 As Double, m1 As Double, m2 As Double, m3 As Double
    Dim s0 As Double, s1 As Double, s2 As Double, s3 As Double
    Dim p0(0 To {max_degree}) As Double
    Dim p1(0 To {max_degree}) As Double
    Dim p2(0 To {max_degree}) As Double
    Dim p3(0 To {max_degree}) As Double
    Dim poly As Double
    Dim P As Double, Q As Double
    Dim log_corr As Double

    Call Features(H, T, d, Uc, L_base, x0_raw, x1_raw, x2_raw, x3_raw)

    Select Case reg
{case_text}    End Select

    If Abs(Q) < EPS_DEN Then
        If Q < 0# Then
            Q = -EPS_DEN
        Else
            Q = EPS_DEN
        End If
    End If

    log_corr = P / Q
    If log_corr < LOG_CORR_MIN Then log_corr = LOG_CORR_MIN
    If log_corr > LOG_CORR_MAX Then log_corr = LOG_CORR_MAX

    PredictRegime = MaxVal(L_base, 0.1#) * Exp(log_corr)
End Function

Private Function MaxVal(ByVal a As Double, ByVal b As Double) As Double
    If a >= b Then
        MaxVal = a
    Else
        MaxVal = b
    End If
End Function

Private Function FnTanh(ByVal x As Double) As Double
    If x > 20# Then
        FnTanh = 1#
    ElseIf x < -20# Then
        FnTanh = -1#
    Else
        FnTanh = (Exp(x) - Exp(-x)) / (Exp(x) + Exp(-x))
    End If
End Function
'''

    with open(filename, 'w', encoding='utf-8') as f:
        f.write(code)
    print(f"[SUCCESS] VBA module saved to '{filename}'")

# ==============================================================================
#  SECTION 5: TRAINING LOGIC
# ==============================================================================
def train_regime(name, X_phys, target_log, config):
    n_samples = len(target_log)
    degree = config['degree']
    lam = config['lam']
    min_fit_samples = max(10, get_term_count(degree) + 5)
    max_filter_iterations = 10
    print(f"   > Training {name} Model (N={n_samples}) [Degree={degree}, Lam={lam}]")

    if n_samples < min_fit_samples:
        n_p = get_term_count(degree)
        params = np.zeros(n_p + n_p - 1)
        params[0] = 0.0
        keep_mask = np.ones(n_samples, dtype=bool)
        return params, np.zeros(4), np.ones(4), keep_mask

    def fit_subset(X_sub, y_sub, init_params=None, max_nfev_second=2000):
        n_sub = len(y_sub)
        X_mean_sub = np.mean(X_sub, axis=0)
        X_std_sub = np.std(X_sub, axis=0)
        X_std_sub[X_std_sub == 0] = 1.0
        X_norm_sub = (X_sub - X_mean_sub) / X_std_sub

        X_poly_sub = expand_poly_jit(X_norm_sub, n_sub, degree)
        X_num_sub = np.ascontiguousarray(X_poly_sub)
        X_den_sub = np.ascontiguousarray(X_poly_sub[:, 1:])

        n_p = X_num_sub.shape[1]
        n_q = X_den_sub.shape[1]
        if init_params is None:
            ATA = X_num_sub.T @ X_num_sub
            ATA[np.diag_indices_from(ATA)] += lam
            ATb = X_num_sub.T @ y_sub
            p_init = np.linalg.solve(ATA, ATb)
            q_init = np.zeros(n_q)
            init_params = np.concatenate([p_init, q_init])

        weights = np.ones(n_sub)
        res = least_squares(
            residuals_jit,
            init_params,
            jac=jacobian_jit,
            args=(X_num_sub, X_den_sub, y_sub, weights),
            method='trf',
            loss='linear',
            f_scale=0.1,
            max_nfev=500,
        )

        num, den = eval_rational(res.x, X_num_sub, X_den_sub)
        den = np.where(np.abs(den) < 0.001, 0.001 * np.sign(den), den)
        den = np.where(den == 0.0, 0.001, den)
        log_pred = np.clip(num / den, -2.5, 2.5)
        errs = np.abs(log_pred - y_sub)
        weights = 1.0 + (errs * 5.0)

        res = least_squares(
            residuals_jit,
            res.x,
            jac=jacobian_jit,
            args=(X_num_sub, X_den_sub, y_sub, weights),
            method='trf',
            loss='linear',
            ftol=1e-15,
            xtol=1e-15,
            gtol=1e-15,
            max_nfev=max_nfev_second,
        )
        return res.x, X_mean_sub, X_std_sub

    def predict_pct_error(params_sub, X_sub, y_sub, mean_sub, std_sub):
        X_norm_sub = (X_sub - mean_sub) / std_sub
        X_poly_sub = expand_poly_jit(X_norm_sub, len(X_sub), degree)
        X_num_sub = np.ascontiguousarray(X_poly_sub)
        X_den_sub = np.ascontiguousarray(X_poly_sub[:, 1:])

        num, den = eval_rational(params_sub, X_num_sub, X_den_sub)
        den = np.where(np.abs(den) < 0.001, 0.001 * np.sign(den), den)
        den = np.where(den == 0.0, 0.001, den)
        log_pred = np.clip(num / den, -2.5, 2.5)
        return np.abs(np.exp(log_pred - y_sub) - 1.0) * 100.0

    active_mask = np.ones(n_samples, dtype=bool)
    params = None
    X_mean = np.zeros(4)
    X_std = np.ones(4)

    for iteration in range(1, max_filter_iterations + 1):
        X_fit = X_phys[active_mask]
        y_fit = target_log[active_mask]
        params, X_mean, X_std = fit_subset(X_fit, y_fit)

        pct_err_fit = predict_pct_error(params, X_fit, y_fit, X_mean, X_std)
        local_keep_mask = pct_err_fit <= ERROR_LIMIT_PCT
        removed_now = int(np.sum(~local_keep_mask))
        kept_now = int(np.sum(local_keep_mask))

        if removed_now == 0:
            if iteration == 1:
                print(f"     -> No fitted samples exceeded {ERROR_LIMIT_PCT:.0f}% final prediction error.")
            else:
                print(
                    f"     -> Filtering converged after {iteration - 1} refit(s); "
                    f"kept {int(np.sum(active_mask))}/{n_samples} samples within {ERROR_LIMIT_PCT:.0f}%."
                )
            return params, X_mean, X_std, active_mask

        if kept_now < min_fit_samples:
            final_keep_mask = active_mask.copy()
            active_indices = np.where(active_mask)[0]
            final_keep_mask[active_indices[~local_keep_mask]] = False
            print(
                f"     -> Iteration {iteration}: removing {removed_now} samples would leave only "
                f"{kept_now}, below the minimum {min_fit_samples}."
            )
            print(
                f"     -> Keeping current model and marking {int(np.sum(final_keep_mask))}/{n_samples} "
                f"samples as fitted for output."
            )
            return params, X_mean, X_std, final_keep_mask

        active_indices = np.where(active_mask)[0]
        active_mask[active_indices[~local_keep_mask]] = False
        print(
            f"     -> Iteration {iteration}: removed {removed_now} samples with final prediction error "
            f"> {ERROR_LIMIT_PCT:.0f}% (kept {int(np.sum(active_mask))}/{n_samples})."
        )

    print(
        f"     -> Reached maximum of {max_filter_iterations} filter iterations; "
        f"keeping {int(np.sum(active_mask))}/{n_samples} samples."
    )
    return params, X_mean, X_std, active_mask
# ==============================================================================
#  SECTION 6: PLOTTING
# ==============================================================================
def plot_diagnostics(df_subset, title, filename):
    if len(df_subset) == 0: return
    plt.figure(figsize=(14, 10))
    plt.suptitle(f"{title} Diagnostics (N={len(df_subset)})", fontsize=16)
    mape = np.mean(df_subset['Error_Pct'])
    
    plt.subplot(2, 2, 1)
    plt.scatter(df_subset['L_fenton'], df_subset['L_pred'], alpha=0.5, s=15, c='blue')
    max_val = max(df_subset['L_fenton'].max(), df_subset['L_pred'].max())
    min_val = min(df_subset['L_fenton'].min(), df_subset['L_pred'].min())
    plt.plot([min_val, max_val], [min_val, max_val], 'k--', lw=2)
    plt.xlabel('True L')
    plt.ylabel('Predicted L')
    plt.title(f'Parity (MAPE={mape:.3f}%)')
    plt.grid(True, alpha=0.3)

    plt.subplot(2, 2, 2)
    plt.scatter(df_subset['L_pred'], df_subset['Error_Pct'], alpha=0.5, s=15, c='purple')
    plt.axhline(0, color='k', linestyle='--')
    plt.xlabel('Predicted L')
    plt.ylabel('Error (%)')
    plt.title('Residuals vs Prediction')
    plt.grid(True, alpha=0.3)

    steepness = df_subset['H'] / df_subset['L_pred']
    plt.subplot(2, 2, 3)
    plt.scatter(steepness, df_subset['Error_Pct'], alpha=0.5, s=15, c='green')
    plt.axhline(0, color='k', linestyle='--')
    plt.xlabel('Steepness (H/L)')
    plt.ylabel('Error (%)')
    plt.title('Residuals vs Steepness')
    plt.grid(True, alpha=0.3)
    
    plt.subplot(2, 2, 4)
    plt.hist(df_subset['Error_Pct'], bins=30, color='orange', edgecolor='black', alpha=0.7)
    plt.xlabel('Error (%)')
    plt.title('Error Distribution')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(filename)
    plt.close()
    print(f"Saved: {filename}")

# ==============================================================================
#  SECTION 7: MAIN
# ==============================================================================
def load_data(filename):
    try:
        df = pd.read_csv(filename, sep='\t')
        if df.shape[1] < 5: df = pd.read_csv(filename, sep=r'\s+')
    except:
        df = pd.read_csv(filename, sep=r'\s+')
    if 's' in df.columns and 'd' not in df.columns:
        df.rename(columns={'s': 'd'}, inplace=True)
    return df

def main():
    print("--- Loading Data ---")
    df = load_data('list.txt')
    H = df['H'].values.astype(np.float64)
    T = df['T'].values.astype(np.float64)
    d = df['d'].values.astype(np.float64)
    Uc = df['Uc'].values.astype(np.float64)
    L_actual = df['L'].values.astype(np.float64)
    
    # 1. Physics Baseline
    L_linear = solve_linear_doppler(T, d, Uc)
    
    # [CRITICAL] Safe Baseline for Training Targets
    L_safe_base = np.maximum(L_linear, 0.1)
    
    # 2. Target
    target_log = np.log(L_actual / L_safe_base)
    
    print("--- Generating Features ---")
    X_phys_all = build_features(H, T, d, Uc, L_safe_base)
    
    # 3. Regimes (Raw Linear)
    rel_depth_lin = d / L_linear  
    mask_shallow = rel_depth_lin < 0.05
    mask_inter   = (rel_depth_lin >= 0.05) & (rel_depth_lin < 0.5)
    mask_deep    = rel_depth_lin >= 0.5
    
    regimes = [('shallow', mask_shallow), ('inter', mask_inter), ('deep', mask_deep)]
    results_dict = {}
    
    print("--- Training Models ---")
    for name, mask in regimes:
        if np.sum(mask) == 0:
            results_dict[name] = {
                'params': np.zeros(1),
                'mean': np.zeros(4),
                'std': np.ones(4),
                'degree': 3,
                'keep_mask': np.zeros(0, dtype=bool),
            }
            continue
        params, mean, std, keep_mask = train_regime(
            name,
            X_phys_all[mask],
            target_log[mask],
            USER_CONFIG[name],
        )
        results_dict[name] = {
            'params': params,
            'mean': mean,
            'std': std,
            'degree': USER_CONFIG[name]['degree'],
            'keep_mask': keep_mask,
        }

    # --- Reconstruction ---
    L_pred = np.zeros_like(L_actual)
    for name, mask in regimes:
        if np.sum(mask) == 0: continue
        deg = results_dict[name]['degree']
        params = results_dict[name]['params']
        mean = results_dict[name]['mean']
        std = results_dict[name]['std']
        
        X_sub = X_phys_all[mask]
        X_norm = (X_sub - mean) / std
        X_poly = expand_poly_jit(X_norm, len(X_norm), deg)
        X_num = np.ascontiguousarray(X_poly)
        X_den = np.ascontiguousarray(X_poly[:, 1:])
        
        num, den = eval_rational(params, X_num, X_den)
        den = np.where(np.abs(den)<0.001, 0.001*np.sign(den), den)
        log_corr = np.clip(num/den, -2.5, 2.5)
        L_pred[mask] = L_safe_base[mask] * np.exp(log_corr)

    df['L_base'] = L_linear
    df['L_fenton'] = L_actual
    df['L_pred'] = L_pred
    df['Error_Pct'] = np.abs((df['L_fenton'] - df['L_pred']) / df['L_fenton']) * 100

    fitted_mask = np.zeros(len(df), dtype=bool)
    for name, mask in regimes:
        if np.sum(mask) == 0:
            continue
        regime_indices = np.where(mask)[0]
        regime_keep_mask = results_dict[name]['keep_mask']
        fitted_mask[regime_indices[regime_keep_mask]] = True

    # Final safety gate: only rows retained by the regime filter AND within the
    # maximum allowed final prediction error are considered fitted outputs.
    fitted_mask &= (df['Error_Pct'].values <= ERROR_LIMIT_PCT)

    df_fited = df.loc[fitted_mask].copy()
    df_not_fited = df.loc[~fitted_mask].copy()

    df_fited[['H', 'T', 'd', 'Uc', 'L']].to_csv(
        'fited.txt',
        sep='	',
        index=False,
        float_format='%.10g',
    )
    print(
        f"Saved fitted-only parameter list to 'fited.txt' "
        f"({len(df_fited)}/{len(df)} rows; columns: H, T, d, Uc, L)."
    )
    print(
        f"Excluded {len(df_not_fited)} rows from fitted outputs because they exceeded "
        f"the {ERROR_LIMIT_PCT:.0f}% final prediction error limit or were removed during refit."
    )

    generate_multi_model_py(results_dict, 'function.py')
    generate_multi_model_vba(results_dict, 'pade.bas')

    if len(df_fited) == 0:
        print("\nNo fitted rows remain after filtering; skipping statistics and plots.")
        return

    mse = np.mean((df_fited['L_fenton'] - df_fited['L_pred'])**2)
    mape = np.mean(df_fited['Error_Pct'])

    print("\n==========================================")
    print(" FINAL STATISTICS (FITTED ONLY)")
    print("==========================================")
    print(f"Global MSE:  {mse:.6f}")
    print(f"Global MAPE: {mape:.2f}%")
    print(f"Global MaxErr: {df_fited['Error_Pct'].max():.2f}%")
    print("-" * 42)

    rel_depth_actual = df_fited['d'] / df_fited['L_fenton']
    stat_mask_shallow = rel_depth_actual < 0.05
    stat_mask_inter = (rel_depth_actual >= 0.05) & (rel_depth_actual < 0.5)
    stat_mask_deep = rel_depth_actual >= 0.5

    def print_regime_stats(name, mask):
        subset = df_fited.loc[mask]
        if len(subset) > 0:
            regime_mse = np.mean((subset['L_fenton'] - subset['L_pred'])**2)
            regime_mape = np.mean(subset['Error_Pct'])
            regime_max = np.max(subset['Error_Pct'])
            count = len(subset)
            print(
                f"{name:<15} (N={count:<4}): MSE={regime_mse:.6f}, "
                f"MAPE={regime_mape:.2f}%, MaxErr={regime_max:.2f}%"
            )
        else:
            print(f"{name:<15} (N=0   ): No samples")

    print_regime_stats('Shallow', stat_mask_shallow)
    print_regime_stats('Intermediate', stat_mask_inter)
    print_regime_stats('Deep', stat_mask_deep)
    print("==========================================")

    print("\n[Generating Visualization...]")
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.scatter(df_fited.loc[stat_mask_deep, 'L_fenton'], df_fited.loc[stat_mask_deep, 'L_pred'], alpha=0.3, s=5, c='blue', label='Deep')
    plt.scatter(df_fited.loc[stat_mask_inter, 'L_fenton'], df_fited.loc[stat_mask_inter, 'L_pred'], alpha=0.3, s=5, c='green', label='Intermediate')
    plt.scatter(df_fited.loc[stat_mask_shallow, 'L_fenton'], df_fited.loc[stat_mask_shallow, 'L_pred'], alpha=0.3, s=5, c='red', label='Shallow')
    max_val = max(df_fited['L_fenton'].max(), df_fited['L_pred'].max())
    plt.plot([0, max_val], [0, max_val], 'k--', linewidth=2, label='Ideal')
    plt.xlabel('True Wavelength (m)')
    plt.ylabel('Predicted Wavelength (m)')
    plt.title(f'Parity Plot (Fitted MAPE={mape:.2f}%)')
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.subplot(1, 2, 2)
    plt.hist(df_fited['Error_Pct'], bins=50, color='purple', alpha=0.7, edgecolor='black')
    plt.xlabel('Absolute Percentage Error (%)')
    plt.ylabel('Count')
    plt.title('Fitted Error Distribution')
    plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('plot_all.png')
    plt.close()
    print("Global Visualization saved as 'plot_all.png'")

    plot_diagnostics(df_fited.loc[stat_mask_shallow], 'Shallow Water (d/L < 0.05)', 'plot_shallow.png')
    plot_diagnostics(df_fited.loc[stat_mask_inter], 'Intermediate Water (0.05 <= d/L < 0.5)', 'plot_intermediate.png')
    plot_diagnostics(df_fited.loc[stat_mask_deep], 'Deep Water (d/L >= 0.5)', 'plot_deep.png')

    print("\n==========================================")
    print(" TOP 20 WORST FITTED PREDICTIONS")
    print("==========================================")
    df_worst = df_fited.sort_values(by='Error_Pct', ascending=False).head(20)
    print(f"{'Idx':<6} | {'H':<5} {'T':<5} {'d':<5} {'Uc':<5} | {'L_base':<9} {'L_true':<9} {'L_pred':<9} | {'Err%':<6}")
    print("-" * 84)
    for idx, row in df_worst.iterrows():
        print(f"{idx:<6} | {row['H']:<5.2f} {row['T']:<5.2f} {row['d']:<5.2f} {row['Uc']:<5.2f} | {row['L_base']:<9.4f} {row['L_fenton']:<9.4f} {row['L_pred']:<9.4f} | {row['Error_Pct']:<6.2f}")

if __name__ == "__main__":
    main()
