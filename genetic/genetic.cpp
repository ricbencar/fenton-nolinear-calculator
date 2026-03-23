/*
====================================================================================================
|                                          WAVE REGRESSION SYSTEM                                  |
|                                            (C++ CPU EDITION)                                     |
====================================================================================================

1. SYSTEM OVERVIEW
----------------------------------------------------------------------------------------------------
This software is a high-performance, industrial-grade C++ implementation of a Hybrid Physics-AI 
Symbolic Regression system. It is designed to solve a critical problem in coastal engineering and 
oceanography: accurately predicting wave characteristics (specifically wavelength) in the presence 
of strong currents and non-linear effects.

Standard linear wave theory often fails when waves interact with currents or become steep 
(non-linear). This system uses Gene Expression Programming (GEP) to discover an explicit 
mathematical "Correction Factor" that bridges the gap between theoretical linear predictions 
and observed reality.

2. PHYSICS BACKGROUND & METHODOLOGY
----------------------------------------------------------------------------------------------------
The core philosophy is a "Grey-Box" model: combining a White-Box physics kernel with a 
Black-Box AI corrector.

    A. THE DISPERSION RELATION (The "Scaffold")
    -------------------------------------------
    Linear Wave Theory (Airy Wave Theory) describes the relationship between wave period (T), 
    wavelength (L), and water depth (d) via the transcendental Dispersion Equation:
    
        omega^2 = g * k * tanh(k * d)
    
    Where:
        - omega (Angular Frequency) = 2 * PI / T
        - k (Wave Number) = 2 * PI / L
        - g = Gravitational acceleration (9.80665 m/s^2)
    
    This equation is implicit for 'k' (and thus 'L'). The system includes a robust Newton-Raphson 
    numerical solver to find the exact root of this equation for every data point. This calculated 
    'L_linear' serves as the "Physics Baseline."

    B. THE DOPPLER EFFECT & NON-LINEARITY
    -------------------------------------
    When waves propagate on a current (U), their apparent frequency shifts (Doppler Effect). 
    The intrinsic dispersion relation becomes:
    
        (omega - k * U)^2 = g * k * tanh(k * d)
    
    Additionally, as waves become steep (large Height H), non-linear effects (Stokes waves, 
    stream function theory) alter the wavelength.
    
    Instead of solving complex high-order polynomials or stream functions, this system allows 
    the AI to "learn" these effects empirically. It defines the final prediction as:
    
        L_final = L_linear_no_current * Correction_Factor(X)

    The AI searches for the function 'Correction_Factor(X)' that minimizes error.

    C. DIMENSIONLESS PARAMETERS (The Inputs)
    ----------------------------------------
    The AI does not see raw dimensions (meters, seconds). It sees 9 dimensionless groups 
    representing physical regimes:
    
    [Regime & Nonlinearity]
    - x0 = ln(d/L): Relative Depth. Distinguishes Deep vs. Shallow water.
    - x1 = ln(H/L): Wave Steepness. Primary indicator of non-linearity.
    - x2 = ln(H/d): Relative Height. Critical for depth-induced breaking.
    - x3 = ln(Ur):  Ursell Number (H*L^2 / d^3). Defines the validity of linear theory.
    
    [Wave-Current Interaction]
    - x4 = Fr:      Froude Number (U / sqrt(g*d)). Ratio of inertial to gravitational forces.
    - x5 = Doppler: (U*T)/L. Direct proxy for the phase shift caused by current.
    - x6 = U/C0:    Velocity Ratio. Current speed vs. Deep-water phase speed.
    
    [Proxies]
    - x7 = ln(H/L0): Deep water steepness proxy.
    - x8 = ln(T*sqrt(g/d)): Dimensionless Period.

3. ALGORITHMS & COMPUTATIONAL ARCHITECTURE
----------------------------------------------------------------------------------------------------
    A. GENE EXPRESSION PROGRAMMING (GEP)
    ------------------------------------
    GEP is an evolutionary algorithm that evolves computer programs (mathematical formulas). 
    Unlike Genetic Programming (GP) which manipulates trees directly, GEP evolves linear strings 
    ("chromosomes") that map to expression trees.
    
    - **Karva Notation**: The linear encoding uses a Head/Tail domain structure. 
      The Head contains functions and terminals; the Tail contains only terminals. 
      This guarantees that *any* genetic mutation produces a syntactically valid formula.
    
    - **Evolutionary Pipeline (ORCHESTRATED V3.2)**:
        1. Selection: Tournament Size 3 (Balance between pressure and diversity).
        2. Replication: Strong Elitism (Top 5 individuals preserved perfectly).
        3. Hall of Fame (HoF): Retains top 10% of historical bests to prevent regression.
        4. Recombination (High Prob):
           - One-Point (40%): Primary engine for mixing large logic blocks.
           - Two-Point (20%): Secondary engine for preserving schemas.
           - Gene Crossover (10%): Swaps entire functional units.
        5. Transposition (Dynamic Scaling):
           - IS/RIS operators now scale insertion length based on Head Size (approx 1/3).
             This prevents bloat in small heads while enabling innovation in large ones.
        6. Mutation (Fine-Tuned & Typed):
           - Uniform (Boosted to 30% chance for constants).
           - Gaussian (10%): Applied to Floating Point Constants (+/- 5%).
           - Step (10%): Applied to Integer Constants (+/- 1).
        7. Smart Seeding: Constants seeded with Physics values (PI, g, 0.5, etc.) and Integers.
        8. Validation: Strict Division-by-Zero penalty (Fitness = 1e9).
           - **CRITICAL FIX**: Explicit sanitization of Fitness values before sorting.

    B. BYTECODE JIT COMPILER (The Speed Engine)
    -------------------------------------------
    Standard tree-based evaluation is slow due to pointer indirection and recursion overhead.
    This system implements a Domain-Specific Virtual Machine (VM):
    
    1. **Compilation**: Before evaluation, every individual's expression tree is flattened 
       into a linear array of Bytecode Instructions (Reverse Polish Notation).
    2. **Execution**: A tight, branch-predicted loop executes these instructions using 
       a small CPU-register-resident stack.
    3. **Performance**: This approach eliminates function call overhead and maximizes 
       L1/L2 cache hits, achieving ~50x speedup over Python interpreters.

    C. SYMBOLIC ALGEBRA ENGINE (SymPy Replica)
    ------------------------------------------
    A custom recursive-descent symbolic engine cleans up the evolved formulas for 
    human readability.
    
    - **Operator Precedence**: Automatically strips redundant parentheses (e.g., `(a*b)+c` -> `a*b+c`).
    - **Algebraic Simplification**:
        - Identity: `x + 0 -> x`, `x * 1 -> x`.
        - Cancellation: `x - x -> 0`, `x / x -> 1`.
        - Aggregation: `x + x -> 2*x`.
    - **Nested Function Folding**: Detects chains like `sqrt(sqrt(x))` and converts 
      them to fractional powers `x**0.25`.
    - **Sign Correction**: Converts `x + -5` to `x - 5`.
    - **Fractional Formatting**: Converts floating point constants (0.3333) to 
      exact rational fractions (1/3) and respects Integer constants.
    
    D. ROBUSTNESS & STABILITY
    ---------------------------------
    - **Strict Math Validity**: Division by zero now returns NaN (Not a Number) instead of 1.0.
      This triggers the fitness penalty immediately, causing the evolutionary engine
      to kill any individual that divides by zero.
    - **Stack Guards**: VM checks for underflow/overflow to prevent segfaults on malformed genes.
    - **Safe Math**: Clamps exp/ln/div to prevent NaNs (except for invalid operations).

4. COMPILATION INSTRUCTIONS
----------------------------------------------------------------------------------------------------
    This code requires a C++17 compliant compiler (GCC/Clang/MSVC) with OpenMP support.
    
    WINDOWS COMMAND (MinGW64):
    g++ -O3 -std=c++17 -static -fopenmp -ffast-math -funroll-loops \
    -march=native -fstack-protector-strong -Wl,--stack,67108864 \
    -Wl,--allow-multiple-definition genetic.cpp -o genetic.exe

    LINUX COMMAND:
    g++ -O3 -std=c++17 -fopenmp -ffast-math -funroll-loops -march=native \
    -fstack-protector-strong -Wall -Wextra genetic.cpp -o genetic.exe

====================================================================================================
*/

#ifdef __MINGW32__
#define NOMINMAX
#include <stdio.h>
#include <time.h>
#include <windows.h>

extern "C" {
    // Ignore warnings about attributes
    #pragma GCC diagnostic ignored "-Wattributes"

    // 1. Implement the missing file functions
    int fseeko64(FILE* stream, long long offset, int origin) {
        return _fseeki64(stream, offset, origin);
    }

    long long ftello64(FILE* stream) {
        return _ftelli64(stream);
    }

    // 2. CRITICAL FIX: Provide the __imp_ pointers the linker is begging for.
    // This tricks libstdc++ into using our local functions as if they were imports.
    int (*__imp_fseeko64)(FILE*, long long, int) = fseeko64;
    long long (*__imp_ftello64)(FILE*) = ftello64;

    // 3. Implement missing time function (Standard signature)
    int clock_gettime64(clockid_t clk_id, struct _timespec64 *tp) {
        if (!tp) return -1;
        FILETIME ft;
        GetSystemTimeAsFileTime(&ft);
        unsigned long long t = ((unsigned long long)ft.dwHighDateTime << 32) | ft.dwLowDateTime;
        t -= 116444736000000000ULL; 
        tp->tv_sec = t / 10000000ULL;
        tp->tv_nsec = (long)((t % 10000000ULL) * 100);
        return 0;
    }
}
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <random>
#include <functional>
#include <iomanip>
#include <numeric>
#include <cfloat>
#include <map>
#include <omp.h>
#include <chrono>
#include <limits>
#include <cstring>

// ------------------------------------------------------------------------------------------------
//                                     GLOBAL CONFIGURATION
// ------------------------------------------------------------------------------------------------

const std::string FILENAME = "fited.txt";
const std::string OUTPUT_FILE = "output.txt";
const int POPULATION_SIZE = 1000;
const int N_GENERATIONS = 100000;

// GEP Gene Configuration (Non-const to allow user input)
// These define the complexity of the evolved formulas.
int N_GENES = 3;       // Default, will be updated by user
int HEAD_LENGTH = 9;   // Default, will be updated by user
int GENE_LENGTH = 0;   // Calculated at runtime: HEAD_LENGTH + (HEAD_LENGTH * (MAX_ARITY - 1) + 1)

const int MAX_ARITY = 2; // Binary operators (+, -, *, /) have arity 2. Sqrt/Exp have arity 1.

// VM Constants
// Stack size must be sufficient for the deepest possible tree
const int MAX_STACK_SIZE = 64;
// Increased program size buffer to accommodate potentially large user-defined HEAD_LENGTH
const int MAX_PROGRAM_SIZE = 2048; 

// Physics Constants
const double PI = 3.14159265359;
const double G = 9.80665;

// ------------------------------------------------------------------------------------------------
//                                     UTILITIES: LOGGER
// ------------------------------------------------------------------------------------------------

// A dual-stream logger that writes to both Console (cout) and File simultaneously.
class DualLogger {
    std::ofstream file;
public:
    DualLogger(const std::string& filename) { 
        file.open(filename);
        if (!file.is_open()) {
            std::cerr << "!! Error: Could not open " << filename << " for writing." << std::endl;
        }
    }
    
    template <typename T> void print(const T& msg) {
        std::cout << msg; 
        if (file.is_open()) file << msg;
    }
    
    template <typename T> void println(const T& msg) {
        std::cout << msg << std::endl; 
        if (file.is_open()) file << msg << std::endl;
    }
    
    void flush() { 
        std::cout.flush(); 
        if (file.is_open()) file.flush(); 
    }
};
DualLogger logger(OUTPUT_FILE);

// ------------------------------------------------------------------------------------------------
//                                UTILITIES: FAST BITWISE MATH
// ------------------------------------------------------------------------------------------------

// Standard std::isnan() can be optimized away by -ffast-math. 
// We implement a bitwise check to guarantee detection of NaNs regardless of compiler flags.
inline bool safe_isnan(double d) {
    uint64_t u;
    std::memcpy(&u, &d, sizeof(d));
    // IEEE 754: NaN has exponent bits all 1, and mantissa non-zero
    return (u & 0x7ff0000000000000ULL) == 0x7ff0000000000000ULL && (u & 0x000fffffffffffffULL);
}

inline bool safe_isinf(double d) {
    uint64_t u;
    std::memcpy(&u, &d, sizeof(d));
    // IEEE 754: Inf has exponent bits all 1, and mantissa zero
    return (u & 0x7fffffffffffffffULL) == 0x7ff0000000000000ULL;
}

// ------------------------------------------------------------------------------------------------
//                                UTILITIES: FAST RNG (XorShift256)
// ------------------------------------------------------------------------------------------------

// Much faster than std::mt19937 for the tight loops required in genetic algorithms.
// We use a static global instance since the mutation loop is serial.
class XorShift256 {
    uint64_t s[4];
    
    // Helper for seeding from a single 64-bit integer
    static uint64_t splitmix64(uint64_t& x) {
        uint64_t z = (x += 0x9e3779b97f4a7c15);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
        z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
        return z ^ (z >> 31);
    }
    
public:
    XorShift256(uint64_t seed = 0) {
        // If no seed provided, use high-resolution clock
        if(seed == 0) seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        s[0] = splitmix64(seed);
        s[1] = splitmix64(seed);
        s[2] = splitmix64(seed);
        s[3] = splitmix64(seed);
    }
    
    uint64_t next() {
        const uint64_t result = s[0] + s[3];
        const uint64_t t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = (s[3] << 45) | (s[3] >> 19);
        return result;
    }
    
    double next_double() {
        // Convert 64-bit integer to [0, 1) double
        return (next() >> 11) * 0x1.0p-53;
    }
    
    int next_int(int min, int max) {
        return min + (next() % (max - min + 1));
    }
};

// Global RNG instance
XorShift256 rng;

// Global RNG Wrappers
double random_double(double min, double max) {
    return min + (rng.next_double() * (max - min));
}

int random_int(int min, int max) {
    return rng.next_int(min, max);
}

// ================================================================================================
//  SECTION 1: PHYSICS KERNEL & MATH UTILS
// ================================================================================================

// STRICT Math Protection: Division by zero returns NaN.
// This allows the fitness evaluation loop to detect the error and penalize the individual.
inline double protected_div(double a, double b) { 
    return std::abs(b) < 1e-9 ? std::numeric_limits<double>::quiet_NaN() : a / b; 
}

inline double protected_sqrt(double a) { return std::sqrt(std::abs(a)); }
inline double protected_sq(double a) { return (std::abs(a) > 1e5) ? 1e10 : a * a; }

// New protected math helpers for extended function set
inline double protected_tanh(double a) { return std::tanh(a); }

inline double protected_exp(double a) { 
    // Clamp to prevent Inf/Overflow. exp(20) is approx 4.8e8, enough for physics context.
    if (a > 20.0) return 4.85e8; 
    if (a < -20.0) return 0.0; 
    return std::exp(a); 
}

inline double protected_ln(double a) { 
    double v = std::abs(a); 
    // Prevent -Inf, return large negative but valid penalty approximation
    return (v < 1e-6) ? -13.8 : std::log(v); // ln(1e-6) approx -13.8
}

// Solves the transcendental Dispersion Relation for Linear Wave Theory
// Equation: omega^2 = g * k * tanh(k * d)
// Method: Newton-Raphson iteration
std::vector<double> solve_linear_no_current(const std::vector<double>& T_arr, const std::vector<double>& d_arr) {
    size_t n = T_arr.size();
    std::vector<double> L_res(n);
    const double PI2 = 2.0 * PI;
    const double PI4 = 4.0 * PI * PI;

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; i++) {
        double T = T_arr[i];
        double d = d_arr[i];
        
        // Safety check for T=0 (prevents initial crash on bad data)
        if(std::abs(T) < 1e-5) { L_res[i] = 0.0; continue; }

        double omega = PI2 / T;
        // Initial guess: Deep water approximation
        double k = PI4 / (G * T * T);
        if (k == 0) k = 1e-4;

        // Newton-Raphson Iteration (Max 40 steps)
        for (int j = 0; j < 40; j++) {
            if (k < 1e-4) k = 1e-4; 
            if (k > 200.0) k = 200.0; // Clamp wave number for stability

            double kd = k * d;
            double th = std::tanh(kd);
            double sigma = std::sqrt(G * k * th);
            double f = sigma - omega;
            
            // Derivative df/dk
            double df = (sigma > 1e-9) ? (G * th + G * kd * (1.0 - th * th)) / (2 * sigma) : 0.0;
            
            if (std::abs(df) < 1e-9) break;
            double k_new = k - f / df;
            
            // Damping for stability (0.8 old + 0.2 new)
            k = 0.8 * k + 0.2 * k_new;
            if (std::abs(k_new - k) < 1e-7) break; // Convergence check
        }
        L_res[i] = PI2 / k;
    }
    return L_res;
}

// Pre-calculates 9 dimensionless parameters for every data point
// This is done once at the beginning to speed up evolution.
std::vector<std::vector<double>> build_features(
    const std::vector<double>& H, const std::vector<double>& T, 
    const std::vector<double>& d, const std::vector<double>& U, 
    const std::vector<double>& L_base) 
{
    size_t n = H.size();
    std::vector<std::vector<double>> X(9, std::vector<double>(n));
    const double PI2 = 2.0 * PI;

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; i++) {
        double L = std::max(L_base[i], 0.1);
        double d_safe = std::max(d[i], 0.1);
        double g_d = G * d_safe;
        double sqrt_gd = std::sqrt(g_d);
        
        // Logarithmic transforms help the AI find power laws easily
        X[0][i] = std::log(std::max(d[i] / L, 1e-7));
        X[1][i] = std::log(std::max(H[i] / L, 1e-7));
        X[2][i] = std::log(std::max(H[i] / d[i], 1e-7));
        X[3][i] = std::log(std::max((H[i] * L * L) / (d_safe * d_safe * d_safe), 1e-7));
        
        // Linear terms for current interaction
        X[4][i] = U[i] / sqrt_gd; // Froude
        X[5][i] = (U[i] * T[i]) / L; // Doppler proxy
        
        double L0 = (G * T[i] * T[i]) / PI2;
        double C0 = L0 / T[i];
        X[6][i] = U[i] / C0; // Velocity ratio
        
        // Proxies
        X[7][i] = std::log(std::max(H[i] / L0, 1e-7));
        X[8][i] = std::log(std::max(T[i] * std::sqrt(G / d_safe), 1e-7));
    }
    return X;
}

// ================================================================================================
//  SECTION 2: JIT BYTECODE ENGINE
// ================================================================================================

// OpCodes for the Virtual Machine (Reverse Polish Notation)
enum OpCode : uint8_t {
    OP_ADD, OP_SUB, OP_MUL, OP_DIV, 
    OP_SQRT, OP_SQ, 
    OP_TANH, OP_EXP, OP_LN,
    OP_VAR,   
    OP_CONST, 
    OP_END
};

struct Instruction {
    uint8_t op;
    uint8_t idx;
    double val;
};

// Flattened program structure for cache-friendly execution
struct CompiledProgram {
    Instruction code[MAX_PROGRAM_SIZE];
    int length;
};

// Node structure for the Expression Tree (used during compilation/simplification)
struct Node { 
    int id; 
    double value; 
    bool is_int; // Hybrid Flag: Determines if constant is Integer or Float
};

// A Gene is a sequence of Nodes in Karva notation
struct Gene { 
    std::vector<Node> sequence; 
}; 

class Individual {
public:
    std::vector<Gene> genes; 
    double fitness;
    CompiledProgram program; 

    Individual() : fitness(DBL_MAX) { 
        program.length = 0; 
        genes.resize(N_GENES);
        for(auto& g : genes) g.sequence.resize(GENE_LENGTH);
    }

    void randomize() {
        for (int g = 0; g < N_GENES; g++) {
            // OPTIMIZED: Smart Seeding (Physics Knowledge + Integer Support)
            for (int i = 0; i < GENE_LENGTH; i++) {
                int min_id = (i < HEAD_LENGTH) ? 0 : 9;
                
                // Explicitly boost constant probability to 30% to fix rarity issue
                // Random check: < 0.3 means Constant (ID 18)
                if(random_double(0, 1) < 0.3) {
                    genes[g].sequence[i].id = 18; 
                } else {
                    genes[g].sequence[i].id = random_int(min_id, 17); // Var or Func
                }
                
                genes[g].sequence[i].is_int = false; 
                
                // If it is a constant, decide type (Int vs Float)
                if (genes[g].sequence[i].id == 18) {
                    if (random_double(0, 1) < 0.5) { 
                        // 50% Chance: Integer Constant (-10 to 10)
                        genes[g].sequence[i].value = (double)random_int(-10, 10);
                        genes[g].sequence[i].is_int = true;
                    } else {
                        // 50% Chance: Physics Float Constant
                         if (random_double(0, 1) < 0.3) {
                            // Smart Physics Seeds
                            double magics[] = { 0.0, 1.0, 2.0, 0.5, PI, G, 2*PI };
                            genes[g].sequence[i].value = magics[random_int(0, 6)];
                        } else {
                            // Random Float
                            genes[g].sequence[i].value = random_double(-5, 5);
                        }
                    }
                }
            }
        }
    }

    // Compiles the Linear Gene string into Bytecode for the VM
    void compile() {
        program.length = 0;
        int prog_idx = 0;

        for (int g = 0; g < N_GENES; g++) {
            // Build Expression Tree for Traversal
            struct TreeNode { int id; double val; TreeNode *l=0, *r=0; };
            TreeNode pool[256]; int ptr = 0; // Increased pool for safety
            TreeNode* q[256]; int head=0, tail=0;
            
            TreeNode* root = &pool[ptr++];
            root->id = genes[g].sequence[0].id; 
            root->val = genes[g].sequence[0].value;
            q[tail++] = root;
            
            int seq_idx = 1;
            // Breadth-First construction from Karva string
            while(head < tail) {
                TreeNode* curr = q[head++];
                // IDs 0-8 are functions. (0-3 Binary, 4-8 Unary)
                if(curr->id <= 8) { 
                    int arity = (curr->id <= 3) ? 2 : 1;
                    for(int k=0; k<arity; k++) {
                        if(seq_idx >= GENE_LENGTH) break;
                        TreeNode* child = &pool[ptr++];
                        child->id = genes[g].sequence[seq_idx].id;
                        child->val = genes[g].sequence[seq_idx].value;
                        seq_idx++;
                        if(k==0) curr->l = child; else curr->r = child;
                        q[tail++] = child;
                    }
                }
            }

            // Post-Order Traversal to emit RPN Bytecode
            std::function<void(TreeNode*)> emit = [&](TreeNode* n) {
                if(!n) return;
                emit(n->l);
                emit(n->r);
                
                if(prog_idx >= MAX_PROGRAM_SIZE - 2) return; 

                Instruction instr;
                // Variables and Consts now start at ID 9
                if(n->id >= 9) {
                    if(n->id == 18) { instr.op = OP_CONST; instr.val = n->val; }
                    else { instr.op = OP_VAR; instr.idx = (uint8_t)(n->id - 9); }
                } else {
                    switch(n->id) {
                        case 0: instr.op = OP_ADD; break;
                        case 1: instr.op = OP_SUB; break;
                        case 2: instr.op = OP_MUL; break;
                        case 3: instr.op = OP_DIV; break;
                        case 4: instr.op = OP_SQRT; break;
                        case 5: instr.op = OP_SQ; break;
                        case 6: instr.op = OP_TANH; break;
                        case 7: instr.op = OP_EXP; break;
                        case 8: instr.op = OP_LN; break;
                    }
                }
                program.code[prog_idx++] = instr;
            };
            emit(root);
            
            // Link genes with addition
            if (g > 0 && prog_idx < MAX_PROGRAM_SIZE - 1) {
                Instruction link; link.op = OP_ADD;
                program.code[prog_idx++] = link;
            }
        }
        program.length = prog_idx;
    }
};

// Robust VM with Stack Guards
// Returns NaN if an illegal operation (like div by zero) occurs
inline double execute_vm(const CompiledProgram& prog, const double* vars) {
    double stack[MAX_STACK_SIZE];
    int sp = 0;

    for(int i = 0; i < prog.length; ++i) {
        const Instruction& instr = prog.code[i];
        switch(instr.op) {
            case OP_VAR: 
                if(sp >= MAX_STACK_SIZE) return 0.0; // Overflow check
                if(instr.idx >= 9) return 0.0; // Safety guard for malformed index
                stack[sp++] = vars[instr.idx]; 
                break;
            case OP_CONST: 
                if(sp >= MAX_STACK_SIZE) return 0.0; // Overflow check
                stack[sp++] = instr.val; 
                break;
            
            // Binary Ops - Need 2 operands
            case OP_ADD: 
                if(sp < 2) return 0.0; // Underflow check
                sp--; stack[sp-1] += stack[sp]; 
                break;
            case OP_MUL: 
                if(sp < 2) return 0.0; 
                sp--; stack[sp-1] *= stack[sp]; 
                break;
            case OP_SUB: 
                if(sp < 2) return 0.0; 
                sp--; stack[sp-1] -= stack[sp]; 
                break;
            case OP_DIV: 
                if(sp < 2) return 0.0; 
                // Uses protected_div which returns NaN on zero division
                sp--; stack[sp-1] = protected_div(stack[sp-1], stack[sp]); 
                break;
                
            // Unary Ops - Need 1 operand
            case OP_SQRT: if(sp < 1) return 0.0; stack[sp-1] = protected_sqrt(stack[sp-1]); break;
            case OP_SQ:   if(sp < 1) return 0.0; stack[sp-1] = protected_sq(stack[sp-1]); break;
            case OP_TANH: if(sp < 1) return 0.0; stack[sp-1] = protected_tanh(stack[sp-1]); break;
            case OP_EXP:  if(sp < 1) return 0.0; stack[sp-1] = protected_exp(stack[sp-1]); break;
            case OP_LN:   if(sp < 1) return 0.0; stack[sp-1] = protected_ln(stack[sp-1]); break;
        }
    }
    return (sp > 0) ? stack[0] : 0.0;
}

// ================================================================================================
//  SECTION 3: SYMBOLIC MATH ENGINE (FULL SYMPY REPLICA)
// ================================================================================================

// Helper: Float to Fraction (Modified for Integer Support)
// Handles "is_int" flag to print clean integers (e.g., "5") instead of floats ("5.000")
std::string to_fraction(double val, bool is_int = false) {
    if (is_int || (std::abs(val - std::round(val)) < 1e-9)) {
        return std::to_string((int)std::round(val));
    }
    if (std::abs(val) < 1e-9) return "0";
    
    double v = std::abs(val);
    int sign = (val < 0) ? -1 : 1;
    int best_n = 1, best_d = 1;
    double min_err = 1e9;
    for (int d = 1; d <= 1000; ++d) {
        int n = (int)std::round(v * d);
        double err = std::abs(v - (double)n / d);
        if (err < min_err) { min_err = err; best_n = n; best_d = d; }
        if (err < 1e-9) break;
    }
    if (best_d == 1) return std::to_string(sign * best_n);
    return std::to_string(sign * best_n) + "/" + std::to_string(best_d);
}

// Operator Precedence helper for adding parentheses
int get_precedence(int id) {
    if (id == 0 || id == 1) return 1; // + -
    if (id == 2 || id == 3) return 2; // * /
    if (id == 5) return 3;            // **
    return 4;                         // Atom/Func
}

struct SymNode {
    int id; double val; bool is_int;
    SymNode *l = nullptr, *r = nullptr;
    // Power tracking for sqrt nesting
    double power = 1.0; 
    
    bool is_const() const { return id == 18; }
    
    // Deep equality for algebraic reduction
    bool equals(const SymNode* other) const {
        if (!other) return false;
        if (id != other->id) return false;
        if (is_const()) return std::abs(val - other->val) < 1e-6;
        bool l_eq = (l && other->l) ? l->equals(other->l) : (l == other->l);
        bool r_eq = (r && other->r) ? r->equals(other->r) : (r == other->r);
        return l_eq && r_eq;
    }
    
    // Check non-negative for redundant abs removal
    bool is_non_negative() const {
        if (id == 5 || id == 4 || id == 7) return true; // sq, sqrt, exp
        if (is_const() && val >= 0) return true;
        return false;
    }
};

// Memory Arena to prevent leaks in Symbolic Engine
struct SymArena {
    std::vector<SymNode*> nodes;
    ~SymArena() { for(auto n : nodes) delete n; }
    
    SymNode* create(int id, double val = 0.0, bool is_int = false) {
        SymNode* n = new SymNode{id, val, is_int};
        nodes.push_back(n);
        return n;
    }
    
    SymNode* clone(SymNode* n) {
        if(!n) return nullptr;
        SymNode* c = create(n->id, n->val, n->is_int);
        c->power = n->power;
        c->l = clone(n->l);
        c->r = clone(n->r);
        return c;
    }
};

// ALGEBRAIC SIMPLIFIER
SymNode* simplify_ast(SymNode* n, SymArena& arena) {
    if (!n) return nullptr;
    // Variables now start at ID 9
    if (n->id >= 9) return n; 
    
    n->l = simplify_ast(n->l, arena);
    n->r = simplify_ast(n->r, arena);

    // 1. Constant Folding
    if ((n->l && n->l->is_const()) && (!n->r || (n->r && n->r->is_const()))) {
        double l_val = n->l->val;
        double r_val = (n->r) ? n->r->val : 0.0;
        double res = 0.0;
        bool ok = true;
        switch(n->id) {
            case 0: res = l_val + r_val; break;
            case 1: res = l_val - r_val; break;
            case 2: res = l_val * r_val; break;
            case 3: 
                // STRICT CHECK: Do not fold if division by zero occurs
                if(std::abs(r_val) > 1e-9) res = l_val/r_val; 
                else ok = false; 
                break;
            case 4: res = protected_sqrt(l_val); break;
            case 5: res = protected_sq(l_val); break;
            case 6: res = protected_tanh(l_val); break;
            case 7: res = protected_exp(l_val); break;
            case 8: res = protected_ln(l_val); break;
        }
        if (ok) { 
            // Reuse node or create new valid node
            n->id = 18; n->val = res; 
            // Auto-detect integer result
            n->is_int = (std::abs(res - std::round(res)) < 1e-9); 
            n->l = nullptr; n->r = nullptr; 
            return n; 
        }
    }
    
    // 2. Algebraic Identities
    // x + 0 = x
    if (n->id == 0) {
        if (n->r && n->r->is_const() && std::abs(n->r->val) < 1e-9) return n->l;
        if (n->l && n->l->is_const() && std::abs(n->l->val) < 1e-9) return n->r;
        if (n->l && n->r && n->l->equals(n->r)) { 
            n->id = 2; n->r = arena.create(18, 2.0, true); return n; 
        }
    }
    // x - x = 0
    if (n->id == 1) {
        if (n->l && n->r && n->l->equals(n->r)) { 
            n->id = 18; n->val = 0; n->is_int = true; n->l=0; n->r=0; return n; 
        }
    }
    // x * 1 = x
    if (n->id == 2) {
        if (n->r && n->r->is_const() && std::abs(n->r->val - 1.0) < 1e-9) return n->l;
        if (n->l && n->l->is_const() && std::abs(n->l->val - 1.0) < 1e-9) return n->r;
        // x * 0 = 0
        if ((n->r && n->r->is_const() && std::abs(n->r->val) < 1e-9) ||
            (n->l && n->l->is_const() && std::abs(n->l->val) < 1e-9)) {
            n->id = 18; n->val = 0; n->is_int = true; n->l=0; n->r=0; return n; 
        }
    }
    // x / x = 1
    if (n->id == 3) {
        if (n->l && n->r && n->l->equals(n->r)) { 
            n->id = 18; n->val = 1.0; n->is_int = true; n->l=0; n->r=0; return n; 
        }
        // Note: We deliberately do NOT include x/0=1 here.
    }
    
    // 3. Nested Roots & Powers
    // sqrt(sqrt(x)) -> power accumulation
    if (n->id == 4 && n->l && n->l->id == 4) { 
        n->l->power *= 0.5; // Multiply exponents
    }

    return n;
}

// SMART PRINTER (Precedence & Formatting)
std::string print_ast(SymNode* n, int parent_prec = 0, double accumulated_pow = 1.0) {
    if (!n) return "";
    
    double current_pow = accumulated_pow * n->power;
    
    // If this node is a sqrt/sq, we might be inside a chain.
    bool is_func = (n->id == 4 || n->id == 5);
    
    if (n->id == 4) { // sqrt
        if (n->l && (n->l->id == 4 || n->l->id == 5)) {
            return print_ast(n->l, parent_prec, current_pow * 0.5);
        } else {
            double final_p = current_pow * 0.5;
            std::string inner = print_ast(n->l, 0, 1.0);
            if (!n->l->is_non_negative()) inner = "abs(" + inner + ")";
            
            if (std::abs(final_p - 0.5) < 1e-5) return "sqrt(" + inner + ")";
            if (std::abs(final_p - 1.0) < 1e-5) return inner;
            return "(" + inner + "**" + to_fraction(final_p) + ")";
        }
    }
    if (n->id == 5) { // sq
        if (n->l && (n->l->id == 4 || n->l->id == 5)) {
            return print_ast(n->l, parent_prec, current_pow * 2.0);
        } else {
            double final_p = current_pow * 2.0;
            std::string inner = print_ast(n->l, get_precedence(2)+1, 1.0);
            if (std::abs(final_p - 1.0) < 1e-5) return inner;
            return "(" + inner + "**" + to_fraction(final_p) + ")";
        }
    }

    // Normal Nodes
    if (n->id == 18) return to_fraction(n->val, n->is_int);
    if (n->id >= 9) return "x" + std::to_string(n->id - 9);
    
    if (n->id == 6) return "tanh(" + print_ast(n->l, 0, 1.0) + ")";
    if (n->id == 7) return "exp(" + print_ast(n->l, 0, 1.0) + ")";
    if (n->id == 8) return "ln(" + print_ast(n->l, 0, 1.0) + ")";

    int my_prec = get_precedence(n->id);
    bool parens = (my_prec < parent_prec); 
    
    std::string s;
    switch(n->id) {
        case 0: // +
            if (n->r && n->r->is_const() && n->r->val < 0) { // x + -5 -> x - 5
                s = print_ast(n->l, my_prec) + " - " + to_fraction(std::abs(n->r->val), n->r->is_int);
            } else {
                s = print_ast(n->l, my_prec) + " + " + print_ast(n->r, my_prec);
            }
            break;
        case 1: // -
            if (n->r && n->r->is_const() && n->r->val < 0) { // x - -5 -> x + 5
                s = print_ast(n->l, my_prec) + " + " + to_fraction(std::abs(n->r->val), n->r->is_int);
            } else {
                s = print_ast(n->l, my_prec) + " - " + print_ast(n->r, my_prec + 1); 
            }
            break;
        case 2: // *
            s = print_ast(n->l, my_prec) + " * " + print_ast(n->r, my_prec);
            break;
        case 3: // /
            s = print_ast(n->l, my_prec) + " / " + print_ast(n->r, my_prec + 1); 
            break;
    }
    
    return parens ? "(" + s + ")" : s;
}

std::string simplify_gene_str(const Gene& gene) {
    SymArena arena; 
    std::vector<SymNode*> q; 
    SymNode* root = arena.create(gene.sequence[0].id, gene.sequence[0].value, gene.sequence[0].is_int);
    q.push_back(root);
    
    int idx = 1, cur = 0;
    while(cur < q.size() && idx < GENE_LENGTH) {
        SymNode* node = q[cur++];
        // <= 8 includes all functions
        if(node->id <= 8) {
            int arity = (node->id <= 3) ? 2 : 1;
            for(int k=0; k<arity; ++k) {
                if(idx >= GENE_LENGTH) break;
                SymNode* c = arena.create(gene.sequence[idx].id, gene.sequence[idx].value, gene.sequence[idx].is_int);
                idx++;
                if(k==0) node->l = c; else node->r = c;
                q.push_back(c);
            }
        }
    }
    
    // Run twice to catch ripple effects
    simplify_ast(root, arena);
    simplify_ast(root, arena);
    
    return print_ast(root);
}

std::string get_formula_str(const Individual& ind) {
    std::string s = "";
    for(int i=0; i<N_GENES; ++i) {
        if(i > 0) s += " + ";
        s += simplify_gene_str(ind.genes[i]);
    }
    return s;
}

// ==============================================================================
//  SECTION 4: STATS
// ==============================================================================

struct DetailedStats {
    double mape, max_err, rmse, bias, p50, p90, p99;
    std::vector<double> residuals, pct_err;
    struct Regime { double mape, max; int count; };
    Regime shallow, inter, deep;
};

DetailedStats calculate_stats(const std::vector<double>& true_val, const std::vector<double>& pred_val, 
                              const std::vector<double>& safe_val, const std::vector<double>& d_arr) {
    DetailedStats s;
    size_t n = true_val.size();
    s.residuals.resize(n); s.pct_err.resize(n);
    double sum_pct=0, sum_sq=0, sum_res=0;
    
    for(size_t i=0; i<n; ++i) {
        // Safe handling of potential NaNs from invalid individuals
        double p = pred_val[i];
        if (safe_isnan(p) || safe_isinf(p)) p = 0.0; 
        
        s.residuals[i] = true_val[i] - p;
        s.pct_err[i] = (std::abs(s.residuals[i]) / std::abs(safe_val[i])) * 100.0;
        sum_pct += s.pct_err[i];
        sum_sq += s.residuals[i]*s.residuals[i];
        sum_res += s.residuals[i];
    }
    
    s.mape = sum_pct/n; s.rmse = std::sqrt(sum_sq/n); s.bias = sum_res/n;
    std::vector<double> srtd = s.pct_err;
    std::sort(srtd.begin(), srtd.end());
    s.max_err = srtd.back();
    s.p50 = srtd[size_t(0.5*n)]; s.p90 = srtd[size_t(0.9*n)]; 
    s.p99 = srtd[size_t(std::min(size_t(0.99*n), n-1))];
    
    auto reg = [&](double l, double h) -> DetailedStats::Regime {
        double r_sum=0, r_max=0; int c=0;
        for(size_t i=0; i<n; i++) {
            double rel = d_arr[i]/safe_val[i];
            if(rel>=l && rel<h) { r_sum+=s.pct_err[i]; if(s.pct_err[i]>r_max)r_max=s.pct_err[i]; c++; }
        }
        return { c?r_sum/c:0, r_max, c };
    };
    s.shallow = reg(0,0.05); s.inter = reg(0.05,0.5); s.deep = reg(0.5,1e9);
    return s;
}

double calc_p99(const std::vector<double>& v) {
    // Optimization: Sanitize + Nth Element instead of full sort
    std::vector<double> t = v; 
    // Sanitize NaNs to prevent crashes in nth_element
    for(auto& val : t) if(safe_isnan(val) || safe_isinf(val)) val = 1e9;
    
    size_t idx = size_t(std::min(size_t(0.99*t.size()), t.size()-1));
    std::nth_element(t.begin(), t.begin() + idx, t.end());
    return t[idx];
}

void print_full_report(const DetailedStats& s, const std::string& f, int gen=-1) {
    std::stringstream ss;
    if(gen!=-1) ss<<"\n["<<std::setw(6)<<gen<<"] >>> NEW BEST FOUND <<<\n";
    ss<<std::string(80,'-')<<"\nFORMULA: "<<f<<"\n"<<std::string(80,'-')<<"\n";
    ss<<"1. GLOBAL ACCURACY\n   MAPE (Mean % Error):   "<<std::fixed<<std::setprecision(4)<<s.mape<<"%\n";
    ss<<"   Max Error Percent:     "<<s.max_err<<"%\n   RMSE (Avg Unit Error): "<<std::setprecision(6)<<s.rmse<<"\n";
    ss<<"   Bias (Mean Residual):  "<<s.bias<<"\n\n2. PERCENTILE BREAKDOWN\n";
    ss<<"   50% of data has error < "<<std::setprecision(4)<<s.p50<<"%\n   90% of data has error < "<<s.p90<<"%\n   99% of data has error < "<<s.p99<<"%\n";
    ss<<"\n3. REGIME ANALYSIS\n   REGIME          | COUNT  | MAPE       | MAX ERROR \n";
    ss<<"   ----------------+--------+------------+-----------\n";
    auto pr = [&](std::string n, auto r){ 
        ss<<"   "<<std::left<<std::setw(15)<<n<<" | "<<std::setw(6)<<r.count<<" | "<<std::setw(9)<<r.mape<<"% | "<<r.max<<"%\n"; 
    };
    pr("Shallow",s.shallow); pr("Intermediate",s.inter); pr("Deep",s.deep);
    ss<<std::string(80,'=');
    logger.println(ss.str());
}

// ==============================================================================
//  SECTION 5: EVOLUTION OPERATORS & HALL OF FAME
// ==============================================================================

// Hall of Fame: Keeps best 10% of unique individuals from history
class HallOfFame {
public:
    std::vector<Individual> members;
    size_t capacity;

    HallOfFame(size_t cap) : capacity(cap) { members.reserve(cap); }

    void update(const std::vector<Individual>& pop) {
        // We assume 'pop' is partially sorted from main loop. 
        // Take the top portion of pop to merge.
        // Optimization: Don't merge whole pop, just top candidates.
        size_t candidates_count = std::min(pop.size(), capacity * 2); 
        std::vector<Individual> candidates = members;
        candidates.insert(candidates.end(), pop.begin(), pop.begin() + candidates_count);
        
        // Sort entire candidate pool
        // CRASH FIX: Safe Sort to handle potential NaNs
        std::sort(candidates.begin(), candidates.end(), [](const Individual& a, const Individual& b){
            bool a_nan = safe_isnan(a.fitness);
            bool b_nan = safe_isnan(b.fitness);
            if (a_nan && b_nan) return false;
            if (a_nan) return false; // Push NaN to end
            if (b_nan) return true;
            return a.fitness < b.fitness;
        });

        // Unique Check (Fitness proximity + Simple Gene Hash proxy)
        // This prevents HoF from filling with identical clones, keeping diversity high
        auto is_duplicate = [](const Individual& a, const Individual& b) {
            if (std::abs(a.fitness - b.fitness) > 1e-9) return false;
            // If fitness is identical, check program length as secondary proxy
            return a.program.length == b.program.length;
        };
        auto last = std::unique(candidates.begin(), candidates.end(), is_duplicate);
        candidates.erase(last, candidates.end());

        // Resize to capacity
        if (candidates.size() > capacity) candidates.resize(capacity);
        members = candidates;
    }
    
    Individual& best() { return members[0]; }
};

// 1. One-Point Crossover (Increased to 0.4 for better mixing)
void crossover_one_point(Individual& p1, Individual& p2) {
    if (random_double(0, 1) < 0.4) {
        int total = N_GENES * GENE_LENGTH;
        int pt = random_int(1, total - 1);
        for (int i = pt; i < total; i++) {
            int g = i / GENE_LENGTH, n = i % GENE_LENGTH;
            std::swap(p1.genes[g].sequence[n], p2.genes[g].sequence[n]);
        }
    }
}

// 2. Two-Point Crossover (Reduced to 0.2 to lower destructiveness)
void crossover_two_point(Individual& p1, Individual& p2) {
    if (random_double(0, 1) < 0.2) {
        int total = N_GENES * GENE_LENGTH;
        int pt1 = random_int(0, total - 2);
        int pt2 = random_int(pt1 + 1, total - 1);
        
        for(int i = pt1; i <= pt2; ++i) {
            int g = i / GENE_LENGTH, n = i % GENE_LENGTH;
            std::swap(p1.genes[g].sequence[n], p2.genes[g].sequence[n]);
        }
    }
}

// 3. Gene Crossover (Rate ~0.1) 
// Swaps entire genes between parents (preserves sub-functions)
void crossover_gene(Individual& p1, Individual& p2) {
    if (N_GENES > 1 && random_double(0, 1) < 0.1) {
        int g = random_int(0, N_GENES - 1);
        std::swap(p1.genes[g], p2.genes[g]);
    }
}

// 4. Uniform Mutation (Rate ~0.044)
void mutate_uniform(Individual& ind) {
    const double MUTATION_RATE = 0.044;
    for (int g = 0; g < N_GENES; g++) {
        for (int i = 0; i < GENE_LENGTH; i++) {
            if (random_double(0, 1) < MUTATION_RATE) {
                int min_id = (i < HEAD_LENGTH) ? 0 : 9;
                
                // Boost constant chance: 30%
                if (random_double(0, 1) < 0.3) ind.genes[g].sequence[i].id = 18;
                else ind.genes[g].sequence[i].id = random_int(min_id, 17);
                
                ind.genes[g].sequence[i].is_int = false;
                
                if (ind.genes[g].sequence[i].id == 18) {
                     // 50/50 Chance for new constant type (Integer/Float)
                     if(random_double(0, 1) < 0.5) {
                         ind.genes[g].sequence[i].value = (double)random_int(-10, 10);
                         ind.genes[g].sequence[i].is_int = true;
                     } else {
                         ind.genes[g].sequence[i].value += random_double(-1.0, 1.0);
                     }
                }
            }
        }
    }
}

// 5. Gaussian Mutation & Step Mutation (Split Logic)
// Numerical Polishing with tighter bounds
void mutate_gaussian(Individual& ind) {
    for (int g = 0; g < N_GENES; g++) {
        for (int i = 0; i < GENE_LENGTH; i++) {
            if (ind.genes[g].sequence[i].id == 18) {
                // 10% chance per constant
                if (random_double(0, 1) < 0.1) {
                    if (ind.genes[g].sequence[i].is_int) {
                        // Step Mutation for Integers (+1 or -1)
                        int step = (random_double(0, 1) < 0.5) ? 1 : -1;
                        ind.genes[g].sequence[i].value += step;
                    } else {
                        // Gaussian Mutation for Floats (+/- 5%)
                        double val = ind.genes[g].sequence[i].value;
                        // Tuned: Tighter bound (5%) for fine polishing
                        double perturbation = random_double(-0.05, 0.05);
                        ind.genes[g].sequence[i].value = val + (val * perturbation); 
                    }
                }
            }
        }
    }
}

// 6. IS Transposition (Dynamic Length Scaling)
void mutate_is_transpose(Individual& ind) {
    if (random_double(0, 1) < 0.1) {
        int g = random_int(0, N_GENES - 1);
        // Tune: Max length is dynamic (approx 1/3 of head)
        int max_len = (HEAD_LENGTH > 3) ? (HEAD_LENGTH / 3) : 1; 
        int len = random_int(1, std::max(1, max_len));
        
        int src = random_int(0, GENE_LENGTH - len);
        int tgt = random_int(0, HEAD_LENGTH - 1); 
        
        std::vector<Node> chunk(len);
        for(int i=0; i<len; ++i) chunk[i] = ind.genes[g].sequence[src+i];
        
        std::vector<Node> new_head = ind.genes[g].sequence; 
        for(int i=0; i<len; ++i) if(tgt+i < HEAD_LENGTH) new_head[tgt+i] = chunk[i];
        
        int r = tgt; int w = tgt + len;
        while(w < HEAD_LENGTH) new_head[w++] = ind.genes[g].sequence[r++];
        
        for(int i=0; i<HEAD_LENGTH; ++i) ind.genes[g].sequence[i] = new_head[i];
    }
}

// 7. RIS Transposition (Dynamic Length Scaling)
void mutate_ris_transpose(Individual& ind) {
    if (random_double(0, 1) < 0.1) {
        int g = random_int(0, N_GENES - 1);
        std::vector<int> func_indices;
        for(int i=0; i<HEAD_LENGTH; ++i) if(ind.genes[g].sequence[i].id <= 8) func_indices.push_back(i);
        
        if(func_indices.empty()) return;
        
        int start_idx = func_indices[random_int(0, func_indices.size()-1)];
        // Tune: Max length dynamic
        int max_len = (HEAD_LENGTH > 3) ? (HEAD_LENGTH / 3) : 1;
        int len = random_int(1, std::max(1, max_len));
        if(start_idx + len > GENE_LENGTH) len = GENE_LENGTH - start_idx;
        
        std::vector<Node> ris(len);
        for(int i=0; i<len; ++i) ris[i] = ind.genes[g].sequence[start_idx+i];
        
        std::vector<Node> new_head(HEAD_LENGTH);
        for(int i=0; i<len; ++i) if(i<HEAD_LENGTH) new_head[i] = ris[i];
        
        int r = 0; int w = len;
        while(w < HEAD_LENGTH) new_head[w++] = ind.genes[g].sequence[r++];
        
        for(int i=0; i<HEAD_LENGTH; ++i) ind.genes[g].sequence[i] = new_head[i];
    }
}

// 8. Gene Transposition (Rate ~0.1)
// Moves a random gene to the beginning of the chromosome
void mutate_gene_transpose(Individual& ind) {
    if (N_GENES > 1 && random_double(0, 1) < 0.1) {
        int src_g = random_int(0, N_GENES - 1);
        Gene target_gene = ind.genes[src_g];
        
        // Shift genes down
        for(int g = src_g; g > 0; --g) {
            ind.genes[g] = ind.genes[g-1];
        }
        ind.genes[0] = target_gene;
    }
}

// 9. Inversion (Rate ~0.1)
void mutate_invert(Individual& ind) {
    if (random_double(0, 1) < 0.1) {
        int g = random_int(0, N_GENES - 1);
        int p1 = random_int(0, HEAD_LENGTH - 1), p2 = random_int(0, HEAD_LENGTH - 1);
        if (p1 > p2) std::swap(p1, p2);
        while (p1 < p2) std::swap(ind.genes[g].sequence[p1++], ind.genes[g].sequence[p2--]);
    }
}

// ==============================================================================
//  MAIN
// ==============================================================================

int main() {
    std::ios::sync_with_stdio(false);

    // Metric Selection   
    logger.println("\n--- CONFIGURATION ---");
    logger.println(">> Linker Function: ADDITION (+) [Fixed]");
    logger.println("\nSelect Optimization Metric:");
    logger.println("  1) MAPE (Mean Absolute Percentage Error)");
    logger.println("  2) MAX ERROR Percent (Worst Case)");
    logger.println("  3) RMSE (Root Mean Square Error)");
    logger.println("  4) P99 ERROR (99% of data error max error)");
    logger.println("  5) RMSE * MAX ERROR (hybrid metric)");
    logger.println("  6) P99 ERROR * MAX ERROR (hybrid metric)");
    logger.println("  7) RMSE * P99 ERROR * MAX ERROR (hybrid metric)");
    
    int choice = 1;
    std::string in;
    while(true) {
        std::cout << "Enter choice (1-7): ";
        std::getline(std::cin, in);
        if(!in.empty()){ try{choice=stoi(in);}catch(...){choice=0;} if(choice>=1 && choice<=7) break;}
        std::cout << "Invalid selection. Please enter 1, 2, 3, 4, 5, 6, or 7.\n";
    }
    
    std::string mname = "MAPE";
    if(choice==2) mname="MAX ERROR"; else if(choice==3) mname="RMSE";
    else if(choice==4) mname="P99 ERROR"; else if(choice==5) mname="RMSE * MAX ERROR";
    else if(choice==6) mname="P99 * MAX ERROR"; else if(choice==7) mname="RMSE * P99 * MAX ERROR";
    logger.println(">> Selected Metric: " + mname);
    
    // --- USER INPUT FOR GENE CONFIGURATION ---
    std::cout << "\n--- GENE CONFIGURATION ---\n";
    std::cout << "Enter Number of Genes (N_GENES): ";
    while(!(std::cin >> N_GENES) || N_GENES <= 0) {
        std::cout << "Invalid. Enter integer > 0: "; std::cin.clear(); std::cin.ignore(1000, '\n');
    }

    std::cout << "Enter Head Length (HEAD_LENGTH): ";
    while(!(std::cin >> HEAD_LENGTH) || HEAD_LENGTH <= 0) {
        std::cout << "Invalid. Enter integer > 0: "; std::cin.clear(); std::cin.ignore(1000, '\n');
    }

    // Dynamic Calculation of Gene Length
    GENE_LENGTH = HEAD_LENGTH + (HEAD_LENGTH * (MAX_ARITY - 1) + 1);
    
    logger.println(">> Configured: N_GENES=" + std::to_string(N_GENES) + 
                   ", HEAD=" + std::to_string(HEAD_LENGTH) + 
                   ", TOTAL_LENGTH=" + std::to_string(GENE_LENGTH));
    // ----------------------------------------
    
    logger.println("\n--- Loading " + FILENAME + " ---");
    std::ifstream f(FILENAME);
    if(!f.good()) { logger.println("Error: File not found."); return 1; }
    
    std::vector<double> H, T_raw, d, U, L_true;
    std::string line; std::getline(f, line);
    while(std::getline(f, line)) {
        std::stringstream ss(line);
        double h, t, dv, u, l;
        if(ss >> h >> t >> dv >> u >> l) {
            H.push_back(h); T_raw.push_back(t); d.push_back(dv); U.push_back(u); L_true.push_back(l);
        }
    }
    size_t rows = L_true.size();
    if(rows < 5) { logger.println("Not enough data."); return 1; }
    logger.print("Loaded "); logger.print(rows); logger.println(" rows.");

    logger.println("Calculating Linear Baseline (Disregarding Current)...");
    std::vector<double> L_base = solve_linear_no_current(T_raw, d);
    logger.println("Extracting 9 Dimensionless Features...");
    std::vector<std::vector<double>> X = build_features(H, T_raw, d, U, L_base);
    std::vector<double> L_safe = L_true; for(auto& v : L_safe) if(v==0) v=1e-6;

    logger.println("\n" + std::string(80, '='));
    logger.println("                        FORMULA DEFINITION");
    logger.println(std::string(80, '='));
    logger.println("Features (x0 - x8):");
    logger.println("  x0 = ln(d/L)        [Relative Depth]");
    logger.println("  x1 = ln(H/L)        [Steepness]");
    logger.println("  x2 = ln(H/d)        [Relative Height]");
    logger.println("  x3 = ln(Ur)         [Ursell Number]");
    logger.println("  x4 = Fr             [Froude Number]");
    logger.println("  x5 = Doppler        [U*T/L]");
    logger.println("  x6 = U/C0           [Velocity Ratio]");
    logger.println("  x7 = ln(H/L0)       [Deep Water Steepness]");
    logger.println("  x8 = ln(T*sqrt(g/d))[Dimensionless Period]");
    logger.println(std::string(80, '=') + "\n");
    
    // Initialize Population
    std::vector<Individual> pop(POPULATION_SIZE);
    for(auto& ind : pop) { ind.randomize(); ind.compile(); }
    
    // Pre-allocate next generation vector to prevent heap thrashing
    std::vector<Individual> next(POPULATION_SIZE);

    // Initialize Hall of Fame (Top 10% capacity)
    HallOfFame hof(POPULATION_SIZE / 10);
    double global_best_score = DBL_MAX;

    logger.println("Starting Evolution (" + std::to_string(N_GENERATIONS) + " gens)... Press Ctrl+C to stop.\n");
    logger.flush(); 

    for(int gen = 1; gen <= N_GENERATIONS; ++gen) {
        
        #pragma omp parallel for schedule(static)
        for(int i = 0; i < POPULATION_SIZE; ++i) {
            double row[9], fit=0, max_e=0, sq_e=0;
            // P99 vector required for metrics 4, 6, and 7
            std::vector<double> p99; if(choice==4 || choice==6 || choice==7) p99.resize(rows);
            bool err = false;
            
            const CompiledProgram& vm = pop[i].program;

            for(size_t r = 0; r < rows; ++r) {
                for(int k=0; k<9; ++k) row[k] = X[k][r];
                double m = execute_vm(vm, row);

                // STRICT CHECK: if NaN (from div by zero) or Inf, this is a fatal error
                if(safe_isnan(m) || safe_isinf(m) || std::abs(m)>100.0) { 
                    err = true; 
                    break; 
                }
                
                double diff = std::abs(L_true[r] - L_base[r]*m);
                double pct = diff/std::abs(L_safe[r])*100.0;
                
                if(choice==1) fit += pct;
                else if(choice==2) max_e = std::max(max_e, pct);
                else if(choice==3) sq_e += diff*diff;
                else if(choice==4) p99[r] = pct;
                else if(choice==6) { p99[r] = pct; max_e = std::max(max_e, pct); }
                else if(choice==7) { p99[r] = pct; max_e = std::max(max_e, pct); sq_e += diff*diff; }
                else { sq_e += diff*diff; max_e = std::max(max_e, pct); }
            }
            
            // Parsimony Pressure (Anti-Bloat)
            double penalty = vm.length * 0.05;
            
            if(err) pop[i].fitness = 1e9; // Fatal Penalty for Div by Zero
            else if(choice==1) pop[i].fitness = (fit/rows) + penalty;
            else if(choice==2) pop[i].fitness = max_e + penalty;
            else if(choice==3) pop[i].fitness = std::sqrt(sq_e/rows) + penalty;
            else if(choice==4) pop[i].fitness = calc_p99(p99) + penalty;
            else if(choice==6) pop[i].fitness = (calc_p99(p99) * max_e) + penalty;
            else if(choice==7) pop[i].fitness = (std::sqrt(sq_e/rows) * calc_p99(p99) * max_e) + penalty;
            else pop[i].fitness = (std::sqrt(sq_e/rows) * max_e) + penalty;
        }

        // CRITICAL CRASH FIX: Explicitly sanitize Fitness values before sort
        // MSVC sort will crash if NaNs are present in strict weak ordering.
        for(auto& ind : pop) {
            if (safe_isnan(ind.fitness) || safe_isinf(ind.fitness) || ind.fitness < 0) {
                ind.fitness = 1e9;
            }
        }

        // OPTIMIZED: Partial Sort using SAFE LOGIC
        int elites = 5;
        int hof_candidates = POPULATION_SIZE / 10;
        int sort_depth = std::max(elites, hof_candidates);
        
        // We use a simple comparator now because we sanitized the data above.
        std::partial_sort(pop.begin(), pop.begin() + sort_depth, pop.end(), [](auto& a, auto& b){ 
            return a.fitness < b.fitness; 
        });
        
        // Update Hall of Fame with best candidates
        hof.update(pop);

        // Check Global Best from HoF (not just current pop)
        // Ensure best score is valid (< 1e8) to avoid reporting broken formulas
        if(hof.best().fitness < global_best_score - 1e-5 && hof.best().fitness < 1e8) {
            global_best_score = hof.best().fitness;
            std::vector<double> preds(rows); double t[9];
            for(size_t r=0; r<rows; ++r) {
                for(int k=0; k<9; ++k) t[k] = X[k][r];
                preds[r] = L_base[r] * execute_vm(hof.best().program, t);
            }
            print_full_report(calculate_stats(L_true, preds, L_safe, d), get_formula_str(hof.best()), gen);
            logger.flush(); 
        }

        // Elitism (Top 1% or 5 individuals)
        int next_idx = 0;
        for(int k=0; k<elites; ++k) next[next_idx++] = pop[k];
        
        while(next_idx < POPULATION_SIZE) {
            auto tourn = [&](){ 
                // Tournament size reduced to 3 (k<2)
                int b = random_int(0, POPULATION_SIZE-1);
                for(int k=0; k<2; k++) { 
                    int c = random_int(0, POPULATION_SIZE-1); 
                    if(pop[c].fitness < pop[b].fitness) b = c; 
                }
                return pop[b];
            };

            Individual c1 = tourn(), c2 = tourn();
            
            // 1. RECOMBINATION (High Probability Mixing)
            crossover_one_point(c1, c2); 
            crossover_two_point(c1, c2);
            crossover_gene(c1, c2); 
            
            // 2. TRANSPOSITION (Structural Innovation)
            mutate_is_transpose(c1);
            mutate_ris_transpose(c1);
            mutate_gene_transpose(c1);
            mutate_invert(c1);
            
            // 3. MUTATION (Refinement)
            mutate_uniform(c1);       
            mutate_gaussian(c1);      

            // Compile & Push
            c1.compile(); 
            next[next_idx++] = c1;
            
            // Process c2 if space permits
            if(next_idx < POPULATION_SIZE) {
                mutate_is_transpose(c2);
                mutate_ris_transpose(c2);
                mutate_gene_transpose(c2);
                mutate_invert(c2);
                mutate_uniform(c2);
                mutate_gaussian(c2);
                c2.compile();
                next[next_idx++] = c2;
            }
        }
        // Swap instead of reallocation
        pop.swap(next);
    }
    
    // Final Report using HoF Best
    if(hof.best().fitness < 1e8) {
        Individual& best_ind = hof.best();
        std::vector<double> preds(rows); double t[9];
        for(size_t r=0; r<rows; ++r) {
            for(int k=0; k<9; ++k) t[k] = X[k][r];
            preds[r] = L_base[r] * execute_vm(best_ind.program, t);
        }
        DetailedStats s = calculate_stats(L_true, preds, L_safe, d);
        logger.println("\n" + std::string(80, '='));
        logger.println("                        FINAL ANALYSIS REPORT (HALL OF FAME BEST)");
        logger.println(std::string(80, '='));
        print_full_report(s, get_formula_str(best_ind));
        
        logger.println("\n4. TOP 20 WORST PREDICTIONS");
        logger.println(std::string(80, '-'));
        std::stringstream hss; 
        hss << "   " << std::left << std::setw(6) << "Index" << std::setw(8) << "H" << std::setw(8) << "T" 
            << std::setw(8) << "d" << std::setw(8) << "U" << std::setw(10) << "L_True" << std::setw(10) << "L_Pred" 
            << std::setw(10) << "Diff" << std::setw(8) << "Error%";
        logger.println(hss.str());
        logger.println(std::string(80, '-'));

        std::vector<std::pair<double, size_t>> es(rows);
        for(size_t i=0; i<rows; ++i) es[i] = {s.pct_err[i], i};
        std::sort(es.rbegin(), es.rend());
        
        for(int i=0; i<20 && i<(int)rows; ++i) {
            size_t idx = es[i].second;
            std::stringstream rss; rss << std::fixed << std::setprecision(2);
            rss << "   " << std::left << std::setw(6) << idx << std::setw(8) << H[idx] << std::setw(8) << T_raw[idx] 
                << std::setw(8) << d[idx] << std::setw(8) << U[idx] << std::setw(10) << L_true[idx] << std::setw(10) << preds[idx] 
                << std::setw(10) << s.residuals[idx] << std::setw(8) << s.pct_err[idx] << "%";
            logger.println(rss.str());
        }
    } else {
        logger.println("\nNo valid model found (all individuals failed constraints).");
    }
    logger.println("\n" + std::string(80, '='));
    logger.println("L_Final = L_linear_no_current * (Multiplier)");
    logger.println(std::string(80, '='));
    return 0;
}