"""
==============================================================================
                    WAVE REGRESSION SYSTEM
==============================================================================

PROGRAM DESCRIPTION:
------------------------------------------------------------------------------
This program performs Symbolic Regression using Gene Expression Programming (GEP)
to discover an explicit mathematical formula for the "Correction Factor" (Multiplier).
This factor corrects standard Linear Wave Theory predictions to account for 
non-linear effects and Doppler shifts caused by currents.

THE HYBRID PHYSICS-AI APPROACH:
1.  **Physics Baseline:** The system first calculates the linear wavelength ($L_{linear}$)
    using the dispersion relation: $\omega^2 = g k \tanh(kd)$. This baseline 
    deliberately ignores the current velocity ($U$) to establish a stable starting point.
    It acts as a "scaffold" for the AI.
    
2.  **AI Correction:** It then evolves a function $f(X)$ using GEP such that:
    $$ L_{final} = L_{linear} \times f(X) $$
    where $f(X)$ is composed of dimensionless parameters representing physical regimes.
    By learning a multiplier rather than the raw value, the AI only needs to learn
    the *deviation* from linear theory, which is much easier and more robust.

FEATURES & DIMENSIONLESS PARAMETERS:
------------------------------------------------------------------------------
[Group 1: Regime & Nonlinearity]
1.  x0 = ln(d/L): Relative Depth (Shallow vs Deep water proxy)
2.  x1 = ln(H/L): Wave Steepness (Non-linearity indicator)
3.  x2 = ln(H/d): Relative Height (Breaking wave indicator)
4.  x3 = ln(Ur):  Ursell Number (Regime classifier)

[Group 2: Wave-Current Interaction]
5.  x4 = Fr:      Froude Number (Inertial vs Gravitational forces)
6.  x5 = Doppler: (U*T)/L (Direct Doppler shift proxy)
7.  x6 = U/C0:    Velocity Ratio (Current relative to Phase Speed)

[Group 3: Input Proxies]
8.  x7 = ln(H/L0): Deep Water Steepness
9.  x8 = ln(T*sqrt(g/d)): Dimensionless Period

OPTIMIZATION OPTIONS:
------------------------------------------------------------------------------
The user can choose to minimize one of the following metrics:
1. MAPE (Mean Absolute Percentage Error) - Good for overall fit.
2. MAX ERROR % - Minimizes the worst-case scenario outlier.
3. RMSE (Root Mean Square Error) - Penalizes large errors heavily.
4. P99 ERROR - Minimizes the error threshold that covers 99% of the data.
5. RMSE * MAX ERROR - A hybrid metric balancing general accuracy and extreme outliers.

INSTALLATION:
    pip install numpy sympy deap geppy numba

==============================================================================
"""

import warnings
# Suppress all warnings for a clean output
warnings.filterwarnings("ignore")
import sys
import operator
import math
import random
import numpy as np
# Configure numpy to ignore divide-by-zero errors (handled manually in code)
np.seterr(all='ignore') 
import sympy as sp
import geppy as gep
from deap import creator, base, tools
from functools import reduce
from numba import vectorize, float64, njit

# --- 1. Global Configuration ---
FILENAME = 'list.txt'        # Input data file
OUTPUT_FILE = 'output.txt'   # Output log file
POPULATION_SIZE = 300        # Size of the GEP population
N_GENERATIONS = 100000       # Max generations (stop with Ctrl+C)
HEAD_LENGTH = 7              # Length of the head of the gene
N_GENES = 4                  # Number of genes per chromosome

# --- Logger Class for Dual Output (Console + File) ---
class DualLogger(object):
    """
    Redirects stdout to both the console and a text file.
    Useful for keeping a permanent record of the evolution process.
    """
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "w", encoding='utf-8')
    
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush() 
        
    def flush(self):
        self.terminal.flush()
        self.log.flush()

# ==============================================================================
#  SECTION 1: PHYSICS KERNEL (LINEAR THEORY - NO CURRENT)
# ==============================================================================

@njit(fastmath=True, cache=True)
def solve_linear_no_current(T_arr, d_arr):
    """
    Computes the Linear Wavelength (L_base) using the Newton-Raphson method
    for the dispersion relation.
    
    CRITICAL NOTE: This function deliberately DISREGARDS the current velocity ($U$).
    It provides the "Physics Baseline" that the AI will correct later.
    """
    n = len(T_arr)
    L_res = np.zeros(n, dtype=np.float64)
    g = 9.80665
    pi = 3.14159265359
    
    for i in range(n):
        T = T_arr[i]
        d = d_arr[i]
        omega = 2 * pi / T
        L0 = (g * T**2) / (2 * pi)
        k = 2 * pi / L0 # Initial guess (Deep water)
            
        # Newton-Raphson Iteration for k (Wave Number)
        for _ in range(40):
            if k < 1e-4: k = 1e-4
            if k > 200.0: k = 200.0
            kd = k * d
            th = np.tanh(kd)
            sigma = np.sqrt(g * k * th)
            f = sigma - omega
            sech2 = 1.0 - th**2
            if sigma > 1e-9:
                d_sigma = (g * th + g * kd * sech2) / (2 * sigma)
            else:
                d_sigma = 0.0
            df = d_sigma
            if abs(df) < 1e-9: break
            k_new = k - f / df
            k = 0.8 * k + 0.2 * k_new # Damped update
            if abs(k_new - k) < 1e-7: break
        
        if k <= 1e-4: k = 1e-4
        L_res[i] = 2 * pi / k
        
    return L_res

@njit(fastmath=True, cache=True)
def build_features_jit(H, T, d, U, L_base):
    """
    Extracts 9 Dimensionless Features from the raw input data.
    These features serve as the inputs (terminals) for the GEP algorithm.
    """
    n = len(H)
    X = np.zeros((9, n), dtype=np.float64) 
    g = 9.80665
    pi = 3.14159265359
    
    for i in range(n):
        L = L_base[i]
        if L < 0.1: L = 0.1
        
        # [Group 1: Regime Parameters]
        X[0, i] = np.log(max(d[i] / L, 1e-7))      # x0: Relative Depth
        X[1, i] = np.log(max(H[i] / L, 1e-7))      # x1: Wave Steepness
        X[2, i] = np.log(max(H[i] / d[i], 1e-7))   # x2: Relative Height
        
        # Ursell Number Calculation
        d_safe = max(d[i], 0.1)
        ur = max((H[i] * L**2) / (d_safe**3), 1e-7)
        X[3, i] = np.log(ur)                       # x3: Log Ursell
        
        # [Group 2: Wave-Current Parameters]
        c_shallow = np.sqrt(g * d_safe)
        X[4, i] = U[i] / c_shallow                 # x4: Froude Number
        X[5, i] = (U[i] * T[i]) / L                # x5: Doppler Shift Proxy
        L0 = (g * T[i]**2) / (2 * pi)
        C0 = L0 / T[i]
        X[6, i] = U[i] / C0                        # x6: Velocity Ratio

        # [Group 3: Input Data Proxies]
        X[7, i] = np.log(max(H[i] / L0, 1e-7))     # x7: Deep Water Steepness
        dim_T = T[i] * np.sqrt(g / d_safe)
        X[8, i] = np.log(max(dim_T, 1e-7))         # x8: Dimensionless Period
        
    return X

# ==============================================================================
#  SECTION 2: NUMBA PRIMITIVES & PROTECTED FUNCTIONS
# ==============================================================================

# Fast vectorized operations for GEP expression tree evaluation
@vectorize([float64(float64, float64)], nopython=True)
def fast_add(x, y): return x + y

@vectorize([float64(float64, float64)], nopython=True)
def fast_sub(x, y): return x - y

@vectorize([float64(float64, float64)], nopython=True)
def fast_mul(x, y): return x * y

@vectorize([float64(float64, float64)], nopython=True)
def protected_div(x, y):
    """Protected division to avoid singularities."""
    if abs(y) < 1e-6: return 1.0
    return x / y

@vectorize([float64(float64)], nopython=True)
def protected_sqrt(x):
    """Protected Sqrt: Returns sqrt of ABSOLUTE value to prevent domain errors."""
    return math.sqrt(abs(x))

@vectorize([float64(float64)], nopython=True)
def sq(x):
    """Protected Square: Clamps extremely large values to prevent overflow."""
    if abs(x) > 1e5: return 1e10
    return x * x

# Linker functions (connect genes to form a chromosome)
def link_add(*args): return sum(args)
def link_mult(*args): return reduce(operator.mul, args, 1)

# ==============================================================================
#  SECTION 3: GP SETUP & FITNESS FUNCTIONS
# ==============================================================================

def setup_gep(linker_func):
    """
    Configures the GEP system:
    - Defines the primitive set (operators + arguments)
    - Sets up the individual structure (genes + head/tail)
    - Registers genetic operators
    """
    input_names = [f'x{i}' for i in range(9)]
    pset = gep.PrimitiveSet('Main', input_names=input_names)
    
    # Registering Operators
    pset.add_function(fast_add, 2, name='add')
    pset.add_function(fast_sub, 2, name='sub')
    pset.add_function(fast_mul, 2, name='mul')
    pset.add_function(protected_div, 2, name='div')
    pset.add_function(protected_sqrt, 1, name='sqrt')
    pset.add_function(sq, 1, name='sq')
    
    # Ephemeral Random Constants (ERCs)
    pset.add_ephemeral_terminal(name='c', gen=lambda: random.uniform(-5, 5))

    # Fitness Definition: We minimize the error (weight = -1.0)
    if not hasattr(creator, 'FitnessMin'):
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
        creator.create("Individual", gep.Chromosome, fitness=creator.FitnessMin)

    toolbox = gep.Toolbox()
    toolbox.register('gene_gen', gep.Gene, pset=pset, head_length=HEAD_LENGTH)
    # Note: Linker is passed here (Addition is hardcoded in Main)
    toolbox.register('individual', creator.Individual, gene_gen=toolbox.gene_gen, n_genes=N_GENES, linker=linker_func)
    toolbox.register('population', tools.initRepeat, list, toolbox.individual)
    toolbox.register('compile', gep.compile_, pset=pset)
    
    return toolbox, pset

# --- JIT Compiled Fitness Functions for Speed ---

@njit(float64(float64[:], float64[:], float64[:], float64[:]))
def calculate_mape_jit(y_true, y_base, y_multiplier, y_safe):
    """Calculates Mean Absolute Percentage Error (MAPE)"""
    n = len(y_true)
    total_err = 0.0
    valid_count = 0
    
    for i in range(n):
        mult = y_multiplier[i]
        # Sanity checks for explosion/NaN
        if math.isinf(mult) or math.isnan(mult) or abs(mult) > 100.0:
            return 1e6 
        
        pred = y_base[i] * mult
        diff = abs(y_true[i] - pred)
        total_err += (diff / abs(y_safe[i]))
        valid_count += 1
        
    if valid_count == 0: return 1e6
    return (total_err / valid_count) * 100.0

@njit(float64(float64[:], float64[:], float64[:], float64[:]))
def calculate_max_err_jit(y_true, y_base, y_multiplier, y_safe):
    """Calculates Maximum Percentage Error (Worst Case)"""
    n = len(y_true)
    max_err = 0.0
    valid_count = 0
    
    for i in range(n):
        mult = y_multiplier[i]
        if math.isinf(mult) or math.isnan(mult) or abs(mult) > 100.0:
            return 1e6 
        
        pred = y_base[i] * mult
        diff = abs(y_true[i] - pred)
        pct = (diff / abs(y_safe[i])) * 100.0
        
        if pct > max_err:
            max_err = pct
        valid_count += 1
        
    if valid_count == 0: return 1e6
    return max_err

@njit(float64(float64[:], float64[:], float64[:], float64[:]))
def calculate_rmse_jit(y_true, y_base, y_multiplier, y_safe):
    """Calculates Root Mean Square Error (RMSE)"""
    n = len(y_true)
    sum_sq_err = 0.0
    valid_count = 0
    
    for i in range(n):
        mult = y_multiplier[i]
        if math.isinf(mult) or math.isnan(mult) or abs(mult) > 100.0:
            return 1e6 
        
        pred = y_base[i] * mult
        diff = y_true[i] - pred
        sum_sq_err += (diff * diff)
        valid_count += 1
        
    if valid_count == 0: return 1e6
    return math.sqrt(sum_sq_err / valid_count)

@njit(float64(float64[:], float64[:], float64[:], float64[:]))
def calculate_p99_jit(y_true, y_base, y_multiplier, y_safe):
    """
    Calculates the 99th Percentile Error.
    This ensures that 99% of the data points have an error below this value.
    """
    n = len(y_true)
    errors = np.zeros(n, dtype=np.float64)
    valid_count = 0
    
    for i in range(n):
        mult = y_multiplier[i]
        if math.isinf(mult) or math.isnan(mult) or abs(mult) > 100.0:
            return 1e6 
            
        pred = y_base[i] * mult
        diff = abs(y_true[i] - pred)
        pct = (diff / abs(y_safe[i])) * 100.0
        errors[i] = pct
        valid_count += 1

    if valid_count == 0: return 1e6
    
    # Sort errors to identify the percentile (Numba optimized sort)
    sorted_errors = np.sort(errors)
    
    # Calculate index for 99th percentile
    # If N=100, we want index 98 (99th item), covering 0-98 (99 items)
    idx = int(math.ceil(0.99 * n)) - 1
    if idx >= n: idx = n - 1
    if idx < 0: idx = 0
    
    return sorted_errors[idx]

@njit(float64(float64[:], float64[:], float64[:], float64[:]))
def calculate_product_score_jit(y_true, y_base, y_multiplier, y_safe):
    """
    Calculates the product of RMSE and Maximum Error Percent.
    Score = RMSE * Max_Error(%)
    Minimizing this balances general accuracy (RMSE) with outlier suppression (Max Error).
    """
    n = len(y_true)
    sum_sq_err = 0.0
    max_err_pct = 0.0
    valid_count = 0
    
    for i in range(n):
        mult = y_multiplier[i]
        if math.isinf(mult) or math.isnan(mult) or abs(mult) > 100.0:
            return 1e9 # Heavy penalty
        
        pred = y_base[i] * mult
        diff = abs(y_true[i] - pred)
        
        # Accumulate RMSE part
        sum_sq_err += (diff * diff)
        
        # Track Max Error part
        pct = (diff / abs(y_safe[i])) * 100.0
        if pct > max_err_pct:
            max_err_pct = pct
            
        valid_count += 1
        
    if valid_count == 0: return 1e9
    
    rmse = math.sqrt(sum_sq_err / valid_count)
    return rmse * max_err_pct

# ==============================================================================
#  SECTION 4: EVALUATION & DETAILED STATS
# ==============================================================================

def calculate_detailed_stats(y_true, y_pred, y_safe, d_arr):
    """
    Calculates comprehensive statistics for the final report.
    Includes global metrics and regime-specific breakdown.
    """
    # 1. Global Errors
    residuals = y_true - y_pred
    abs_res = np.abs(residuals)
    pct_err = (abs_res / np.abs(y_safe)) * 100.0
    
    mape = np.mean(pct_err)
    max_err = np.max(pct_err)
    rmse = np.sqrt(np.mean(residuals**2))
    bias = np.mean(residuals)
    std_dev = np.std(residuals)
    
    # 2. Percentiles
    p50 = np.percentile(pct_err, 50)
    p90 = np.percentile(pct_err, 90)
    p95 = np.percentile(pct_err, 95)
    p99 = np.percentile(pct_err, 99)
    
    # 3. Regime Classification (using d/L_true for categorization)
    rel_depth = d_arr / y_safe
    mask_shallow = rel_depth < 0.05
    mask_inter   = (rel_depth >= 0.05) & (rel_depth < 0.5)
    mask_deep    = rel_depth >= 0.5
    
    def get_regime_metrics(mask):
        if np.sum(mask) == 0: return 0.0, 0.0, 0
        r_pct = pct_err[mask]
        return np.mean(r_pct), np.max(r_pct), np.sum(mask)

    shallow_mape, shallow_max, n_s = get_regime_metrics(mask_shallow)
    inter_mape, inter_max, n_i = get_regime_metrics(mask_inter)
    deep_mape, deep_max, n_d = get_regime_metrics(mask_deep)
    
    stats = {
        'mape': mape, 'max_err': max_err, 'rmse': rmse, 'bias': bias, 'std_dev': std_dev,
        'p50': p50, 'p90': p90, 'p95': p95, 'p99': p99,
        'shallow': (shallow_mape, shallow_max, n_s),
        'inter': (inter_mape, inter_max, n_i),
        'deep': (deep_mape, deep_max, n_d),
        'pct_err': pct_err,
        'residuals': residuals
    }
    return stats

def simplify_formula(best_ind):
    """
    Converts GEP individual to a readable symbolic formula.
    Includes fix for SQRT arguments to ensure mathematical validity.
    """
    sym_map = {
        'add': operator.add, 
        'sub': operator.sub, 
        'mul': operator.mul, 
        'div': lambda x,y: x/y,
        # THE FIX: Wrap argument in Abs so formula output is valid in standard math tools
        'sqrt': lambda x: sp.sqrt(sp.Abs(x)), 
        'sq': lambda x: x**2
    }
    
    try:
        simple_expr = gep.simplify(best_ind, symbolic_function_map=sym_map)
        # Limit float precision for readability
        simple_expr = simple_expr.replace(
            lambda x: x.is_Float,
            lambda x: sp.Rational(x).limit_denominator(100)
        )
        return simple_expr
    except:
        return str(best_ind)

def print_full_report(stats, formula_str, gen_num=None):
    if gen_num is not None:
        print(f"\n[{gen_num:6d}] >>> NEW BEST FOUND <<<")
    
    print("-" * 80)
    print(f"FORMULA: {formula_str}")
    print("-" * 80)
    
    print(f"1. GLOBAL ACCURACY")
    print(f"   MAPE (Mean % Error):   {stats['mape']:.4f}%")
    print(f"   Max Error Percent:     {stats['max_err']:.4f}%")
    print(f"   RMSE (Avg Unit Error): {stats['rmse']:.6f}")
    print(f"   Bias (Mean Residual):  {stats['bias']:.6f}")
    
    print(f"\n2. PERCENTILE BREAKDOWN")
    print(f"   50% of data has error < {stats['p50']:.4f}%")
    print(f"   90% of data has error < {stats['p90']:.4f}%")
    print(f"   99% of data has error < {stats['p99']:.4f}%")
    
    print(f"\n3. REGIME ANALYSIS")
    print(f"   {'REGIME':<15} | {'COUNT':<6} | {'MAPE':<10} | {'MAX ERROR':<10}")
    print(f"   {'-'*15}-+-{'-'*6}-+-{'-'*10}-+-{'-'*10}")
    print(f"   {'Shallow':<15} | {stats['shallow'][2]:<6} | {stats['shallow'][0]:<9.4f}% | {stats['shallow'][1]:<9.4f}%")
    print(f"   {'Intermediate':<15} | {stats['inter'][2]:<6} | {stats['inter'][0]:<9.4f}% | {stats['inter'][1]:<9.4f}%")
    print(f"   {'Deep':<15} | {stats['deep'][2]:<6} | {stats['deep'][0]:<9.4f}% | {stats['deep'][1]:<9.4f}%")
    print("="*80)

# ==============================================================================
#  SECTION 5: MAIN EXECUTION
# ==============================================================================

def main():
    sys.stdout = DualLogger(OUTPUT_FILE)
    
    # --- CONFIGURATION: Linker & Metric ---
    print("\n--- CONFIGURATION ---")
    
    # 1. LINKER: Hardcoded to Addition as requested
    global LINKER_FUNC
    LINKER_FUNC = link_add
    print(">> Linker Function: ADDITION (+) [Fixed]")

    # 2. METRIC: User Selection
    print("\nSelect Optimization Metric:")
    print("  1) MAPE (Mean Absolute Percentage Error)")
    print("  2) Max Error Percent (Worst Case)")
    print("  3) RMSE (Root Mean Square Error)")
    print("  4) P99 ERROR (99% of data error max error)")
    print("  5) RMSE * MAX ERROR (Combined Hybrid Metric)")
    
    metric_choice = '1'
    while True:
        metric_choice = input("Enter choice (1, 2, 3, 4, or 5): ").strip()
        if metric_choice in ['1', '2', '3', '4', '5']:
            break
        print("Invalid selection. Please enter 1, 2, 3, 4, or 5.")

    metric_name = "MAPE"
    if metric_choice == '2': metric_name = "MAX ERROR"
    elif metric_choice == '3': metric_name = "RMSE"
    elif metric_choice == '4': metric_name = "P99 ERROR (99th Percentile)"
    elif metric_choice == '5': metric_name = "RMSE * MAX ERROR"
    
    print(f">> Selected Metric: {metric_name}")

    print(f"\n--- Loading {FILENAME} ---")
    try:
        with open(FILENAME, 'r') as f: lines = f.readlines()
    except FileNotFoundError: sys.exit(1)

    # Parsing Headers
    idx_map = {'H':0, 'T':1, 'd':2, 'U':3, 'L':4}
    if len(lines) > 0:
        headers = lines[0].strip().replace('\t', ' ').split()
        for col in idx_map.keys():
            for h in headers:
                if h == col or (col == 'd' and h == 'D'): idx_map[col] = headers.index(h)

    # Load Data
    H, T, d, U, L = [], [], [], [], []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) < 5: continue
        try:
            H.append(float(parts[idx_map['H']]))
            T.append(float(parts[idx_map['T']]))
            d.append(float(parts[idx_map['d']]))
            U.append(float(parts[idx_map['U']]))
            L.append(float(parts[idx_map['L']]))
        except: continue

    if len(L) < 5: sys.exit("Not enough data.")
    
    # Convert to Numpy Arrays for Speed
    H_arr = np.array(H, dtype=np.float64)
    T_arr = np.array(T, dtype=np.float64)
    d_arr = np.array(d, dtype=np.float64)
    U_arr = np.array(U, dtype=np.float64)
    L_true = np.array(L, dtype=np.float64)
    L_safe = np.where(L_true == 0, 1e-6, L_true)

    print(f"Loaded {len(L_true)} rows.")
    
    # 1. Calculate Physics Baseline
    print("Calculating Linear Baseline (Disregarding Current)...")
    L_base = solve_linear_no_current(T_arr, d_arr)
    
    # 2. Extract Features
    print("Extracting 9 Dimensionless Features...")
    X_feats = build_features_jit(H_arr, T_arr, d_arr, U_arr, L_base)
    
    # Unpack for GEP
    x_args = tuple(X_feats[i] for i in range(9))

    print("\n" + "="*80)
    print("                        FORMULA DEFINITION")
    print("="*80)
    print("Features (x0 - x8):")
    print("  x0 = ln(d/L)        [Relative Depth]")
    print("  x1 = ln(H/L)        [Steepness]")
    print("  x2 = ln(H/d)        [Relative Height]")
    print("  x3 = ln(Ur)         [Ursell Number]")
    print("  x4 = Fr             [Froude Number]")
    print("  x5 = Doppler        [U*T/L]")
    print("  x6 = U/C0           [Velocity Ratio]")
    print("  x7 = ln(H/L0)       [Deep Water Steepness]")
    print("  x8 = ln(T*sqrt(g/d))[Dimensionless Period]")
    print("="*80 + "\n")
    
    toolbox, pset = setup_gep(LINKER_FUNC)

    # Dynamic Evaluation Function based on user choice
    def evaluate(individual):
        func = toolbox.compile(individual)
        try:
            multiplier = func(*x_args)
            
            # Select the appropriate JIT function
            if metric_choice == '1':
                # MAPE
                score = calculate_mape_jit(L_true, L_base, multiplier, L_safe)
            elif metric_choice == '2':
                # Max Error
                score = calculate_max_err_jit(L_true, L_base, multiplier, L_safe)
            elif metric_choice == '3':
                # RMSE
                score = calculate_rmse_jit(L_true, L_base, multiplier, L_safe)
            elif metric_choice == '4':
                # 99th Percentile Error (P99)
                score = calculate_p99_jit(L_true, L_base, multiplier, L_safe)
            else:
                # RMSE * MAX ERROR (Combined)
                score = calculate_product_score_jit(L_true, L_base, multiplier, L_safe)
                
            return (score,)
        except:
            return (1e6,)

    toolbox.register('evaluate', evaluate)
    toolbox.register('select', tools.selTournament, tournsize=5)
    toolbox.register('mut_uniform', gep.mutate_uniform, pset=pset, ind_pb=0.1, pb=1)
    toolbox.register('mut_invert', gep.invert, pb=0.1)
    toolbox.register('mut_is_transpose', gep.is_transpose, pb=0.1)
    toolbox.register('cx_1p', gep.crossover_one_point, pb=0.4)
    
    stats = tools.Statistics(key=lambda ind: ind.fitness.values[0])
    stats.register("min", np.min)

    pop = toolbox.population(n=POPULATION_SIZE)
    hof = tools.HallOfFame(1)
    
    # Evolution Loop with Improvement Detection
    best_score_so_far = float('inf')
    
    print(f"Starting Evolution ({N_GENERATIONS} gens)... Press Ctrl+C to stop.\n")

    try:
        for gen in range(1, N_GENERATIONS + 1):
            # Run 1 generation silently
            pop, log = gep.gep_simple(pop, toolbox, n_generations=1, n_elites=3, 
                                      stats=stats, hall_of_fame=hof, verbose=False)
            
            current_best = hof[0]
            current_score = current_best.fitness.values[0]
            
            # CHECK FOR IMPROVEMENT
            if current_score < best_score_so_far - 1e-5: 
                best_score_so_far = current_score
                
                # Calculate full stats for this new best (Always show full stats for report)
                func = toolbox.compile(current_best)
                mult = func(*x_args)
                L_pred = L_base * mult
                s = calculate_detailed_stats(L_true, L_pred, L_safe, d_arr)
                expr_str = simplify_formula(current_best)
                
                # Print FULL Report immediately
                print_full_report(s, expr_str, gen_num=gen)
                
    except KeyboardInterrupt:
        print("\n>>> INTERRUPTED BY USER <<<")

    if not hof: sys.exit(0)
    best = hof[0]
    
    # Final Report
    func = toolbox.compile(best)
    multiplier = func(*x_args)
    L_pred = L_base * multiplier
    s = calculate_detailed_stats(L_true, L_pred, L_safe, d_arr)
    final_expr = simplify_formula(best)

    print("\n" + "="*80)
    print("                        FINAL ANALYSIS REPORT")
    print("="*80)
    print_full_report(s, final_expr)

    print(f"\n4. TOP 20 WORST PREDICTIONS")
    print("-" * 80)
    print(f"   {'Index':<6} {'H':<8} {'T':<8} {'d':<8} {'U':<8} {'L_True':<10} {'L_Pred':<10} {'Diff':<10} {'Error%':<8}")
    print("-" * 80)
    
    sorted_idx = np.argsort(s['pct_err'])[::-1]
    for i in sorted_idx[:20]:
        print(f"   {i:<6} {H_arr[i]:<8.2f} {T_arr[i]:<8.2f} {d_arr[i]:<8.2f} {U_arr[i]:<8.2f} "
              f"{L_true[i]:<10.2f} {L_pred[i]:<10.2f} {s['residuals'][i]:<10.2f} {s['pct_err'][i]:<8.2f}%")

    print("\n" + "="*80)
    print("L_Final = L_linear_no_current * (Multiplier)")
    print("="*80)

if __name__ == "__main__":
    main()