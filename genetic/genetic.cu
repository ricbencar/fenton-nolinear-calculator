/*
====================================================================================================
|                                     WAVE REGRESSION SYSTEM                                       |
|                               (RAW CUDA + SHARED MEMORY EDITION)                                 |
====================================================================================================

1. SYSTEM OVERVIEW
----------------------------------------------------------------------------------------------------
This software is a high-performance, industrial-grade C++/CUDA implementation of a Hybrid Physics-AI 
Symbolic Regression system. It is designed to solve a critical problem in coastal engineering and 
oceanography: accurately predicting wave characteristics (specifically wavelength) in the presence 
of strong currents and non-linear effects.

Standard linear wave theory often fails when waves interact with currents or become steep 
(non-linear). This system uses Gene Expression Programming (GEP) to discover an explicit 
mathematical "Correction Factor" that bridges the gap between theoretical linear predictions 
and observed reality.

2. GPU ARCHITECTURE
----------------------------------------------------------------------------------------------------
Unlike standard genetic algorithms that run on CPU, this system offloads the heavy mathematical 
evaluation to the GPU using a custom "Virtual Machine" kernel.

    A. RAW CUDA KERNEL & SHARED MEMORY
    ----------------------------------
    We utilize a "One Block per Individual" topology.
    - **Instruction Caching**: At the start of execution, the block cooperatively loads the 
      individual's bytecode into **Shared Memory (L1 Cache)**. This reduces instruction fetch 
      latency from ~400 cycles (Global Memory) to ~2 cycles (Shared Memory).
    - **Structure of Arrays (SoA)**: Feature data is stored in column-major format in VRAM, 
      allowing threads to perform coalesced memory reads (32 threads reading continuous bytes), 
      saturating the memory bandwidth.
    - **Lock-Step Execution**: Threads in a warp execute the same instructions on different 
      data rows, minimizing divergence.

    B. PRECISION
    ------------
    The GPU kernel operates in Single Precision (float). This provides a 2x-32x speedup over 
    Double Precision on consumer GPUs (GeForce/RTX), which is sufficient for the chaotic 
    regime of wave mechanics.

3. PHYSICS BACKGROUND & METHODOLOGY
----------------------------------------------------------------------------------------------------
The core philosophy is a "Grey-Box" model: combining a White-Box physics kernel with a 
Black-Box AI corrector.

    A. THE DISPERSION RELATION (The "Scaffold")
    -------------------------------------------
    Linear Wave Theory (Airy Wave Theory) describes the relationship between wave period (T), 
    wavelength (L), and water depth (d) via the transcendental Dispersion Equation:
    
        omega^2 = g * k * tanh(k * d)
    
    This equation is implicit for 'k' (and thus 'L'). The system includes a robust Newton-Raphson 
    numerical solver to find the exact root of this equation for every data point. This calculated 
    'L_linear' serves as the "Physics Baseline."

    B. THE DOPPLER EFFECT & NON-LINEARITY
    -------------------------------------
    When waves propagate on a current (U), their apparent frequency shifts (Doppler Effect). 
    Additionally, as waves become steep (large Height H), non-linear effects (Stokes waves) 
    alter the wavelength.
    
    The AI searches for the function 'Correction_Factor(X)' that minimizes error:
    
        L_final = L_linear_no_current * Correction_Factor(X)

    C. DIMENSIONLESS PARAMETERS (The Inputs)
    ----------------------------------------
    The AI sees 9 dimensionless groups representing physical regimes:
    
    [Regime & Nonlinearity]
    - x0 = ln(d/L): Relative Depth.
    - x1 = ln(H/L): Wave Steepness.
    - x2 = ln(H/d): Relative Height.
    - x3 = ln(Ur):  Ursell Number (H*L^2 / d^3).
    
    [Wave-Current Interaction]
    - x4 = Fr:      Froude Number.
    - x5 = Doppler: (U*T)/L.
    - x6 = U/C0:    Velocity Ratio.
    
    [Proxies]
    - x7 = ln(H/L0): Deep water steepness proxy.
    - x8 = ln(T*sqrt(g/d)): Dimensionless Period.

4. COMPILATION INSTRUCTIONS
----------------------------------------------------------------------------------------------------
    This code requires the NVIDIA CUDA Toolkit (nvcc).
    
    COMPILE COMMAND:
    nvcc -O3 -Xptxas -O3,-v -Xcompiler "-O3 -march=native -DNDEBUG" \
    -arch=native --use_fast_math -extra-device-vectorization \
    --cudart static -std=c++17 genetic.cu -o genetic.exe

====================================================================================================
*/

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
#include <chrono>
#include <limits>
#include <cstring>

// CUDA Runtime
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// ------------------------------------------------------------------------------------------------
//                                     CUDA ERROR HANDLING
// ------------------------------------------------------------------------------------------------

#define checkCuda(call) { \
    const cudaError_t error = call; \
    if (error != cudaSuccess) { \
        std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << ", " \
                  << cudaGetErrorString(error) << std::endl; \
        exit(1); \
    } \
}

// ------------------------------------------------------------------------------------------------
//                                     GLOBAL CONFIGURATION
// ------------------------------------------------------------------------------------------------

const std::string FILENAME = "list.txt";
const std::string OUTPUT_FILE = "output.txt";
const int POPULATION_SIZE = 1000;
const int N_GENERATIONS = 1000000;

// GEP Gene Configuration (Non-const to allow user input)
int N_GENES = 3;       // Default, will be updated by user
int HEAD_LENGTH = 9;   // Default, will be updated by user
int GENE_LENGTH = 0;   // Calculated at runtime: HEAD_LENGTH + (HEAD_LENGTH * (MAX_ARITY - 1) + 1)

const int MAX_ARITY = 2; // Binary operators (+, -, *, /) have arity 2. Sqrt/Exp have arity 1.

// VM Constants
// Stack size must be sufficient for the deepest possible tree
const int MAX_STACK_SIZE = 64;
// Increased program size buffer to accommodate potentially large user-defined HEAD_LENGTH
const int MAX_PROGRAM_SIZE = 2048; 

// Physics Constants (Float for GPU, Double for CPU setup)
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
//                                UTILITIES: FAST BITWISE MATH (Host Side)
// ------------------------------------------------------------------------------------------------

// Standard std::isnan() can be optimized away by -ffast-math. 
inline bool host_safe_isnan(double d) {
    uint64_t u;
    std::memcpy(&u, &d, sizeof(d));
    // IEEE 754: NaN has exponent bits all 1, and mantissa non-zero
    return (u & 0x7ff0000000000000ULL) == 0x7ff0000000000000ULL && (u & 0x000fffffffffffffULL);
}

inline bool host_safe_isinf(double d) {
    uint64_t u;
    std::memcpy(&u, &d, sizeof(d));
    return (u & 0x7fffffffffffffffULL) == 0x7ff0000000000000ULL;
}

// ------------------------------------------------------------------------------------------------
//                                UTILITIES: FAST RNG (XorShift256)
// ------------------------------------------------------------------------------------------------

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

XorShift256 rng;

double random_double(double min, double max) {
    return min + (rng.next_double() * (max - min));
}

int random_int(int min, int max) {
    return rng.next_int(min, max);
}

// ================================================================================================
//  SECTION 1: PHYSICS KERNEL & FEATURE EXTRACTION (HOST)
// ================================================================================================

// Solves the transcendental Dispersion Relation for Linear Wave Theory
// Equation: omega^2 = g * k * tanh(k * d)
// Method: Newton-Raphson iteration
std::vector<double> solve_linear_no_current(const std::vector<double>& T_arr, const std::vector<double>& d_arr) {
    size_t n = T_arr.size();
    std::vector<double> L_res(n);
    const double PI2 = 2.0 * PI;
    const double PI4 = 4.0 * PI * PI;

    for (size_t i = 0; i < n; i++) {
        double T = T_arr[i];
        double d = d_arr[i];
        
        // Safety check for T=0
        if(std::abs(T) < 1e-5) { L_res[i] = 0.0; continue; }

        double omega = PI2 / T;
        // Initial guess: Deep water approximation
        double k = PI4 / (G * T * T);
        if (k == 0) k = 1e-4;

        // Newton-Raphson Iteration (Max 40 steps)
        for (int j = 0; j < 40; j++) {
            if (k < 1e-4) k = 1e-4; 
            if (k > 200.0) k = 200.0; 

            double kd = k * d;
            double th = std::tanh(kd);
            double sigma = std::sqrt(G * k * th);
            double f = sigma - omega;
            
            double df = (sigma > 1e-9) ? (G * th + G * kd * (1.0 - th * th)) / (2 * sigma) : 0.0;
            
            if (std::abs(df) < 1e-9) break;
            double k_new = k - f / df;
            
            k = 0.8 * k + 0.2 * k_new;
            if (std::abs(k_new - k) < 1e-7) break; 
        }
        L_res[i] = PI2 / k;
    }
    return L_res;
}

// Builds Feature Matrix (transposed for GPU SoA layout later)
// Returns a flat vector in preparation for GPU transfer (Features * Rows)
std::vector<double> build_soa_features_float(
    const std::vector<double>& H, const std::vector<double>& T, 
    const std::vector<double>& d, const std::vector<double>& U, 
    const std::vector<double>& L_base) 
{
    size_t n = H.size();
    std::vector<double> soa(n * 9); 
    const double PI2 = 2.0 * PI;

    for (size_t i = 0; i < n; i++) {
        double L = std::max(L_base[i], 0.1);
        double d_safe = std::max(d[i], 0.1);
        double g_d = G * d_safe;
        double sqrt_gd = std::sqrt(g_d);
        
        // Logarithmic transforms help the AI find power laws easily
        soa[0 * n + i] = std::log(std::max(d[i] / L, 1e-7));
        soa[1 * n + i] = std::log(std::max(H[i] / L, 1e-7));
        soa[2 * n + i] = std::log(std::max(H[i] / d[i], 1e-7));
        soa[3 * n + i] = std::log(std::max((H[i] * L * L) / (d_safe * d_safe * d_safe), 1e-7));
        
        // Linear terms for current interaction
        soa[4 * n + i] = U[i] / sqrt_gd; // Froude
        soa[5 * n + i] = (U[i] * T[i]) / L; // Doppler proxy
        
        double L0 = (G * T[i] * T[i]) / PI2;
        double C0 = L0 / T[i];
        soa[6 * n + i] = U[i] / C0; // Velocity ratio
        
        // Proxies
        soa[7 * n + i] = std::log(std::max(H[i] / L0, 1e-7));
        soa[8 * n + i] = std::log(std::max(T[i] * std::sqrt(G / d_safe), 1e-7));
    }
    return soa;
}

// ================================================================================================
//  SECTION 2: GPU DEVICE FUNCTIONS & KERNEL
// ================================================================================================

// Shared OpCodes
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
    float val; // GPU uses float
};

// Device Math Intrinsics (Float optimized)
// Using native CUDA intrinsics where possible
__device__ inline float d_protected_div(float a, float b) { 
    return fabsf(b) < 1e-6f ? NAN : __fdividef(a, b); 
}
__device__ inline float d_protected_sqrt(float a) { return sqrtf(fabsf(a)); }
__device__ inline float d_protected_sq(float a) { return (fabsf(a) > 1e4f) ? 1e8f : a * a; }
__device__ inline float d_protected_tanh(float a) { return tanhf(a); }
__device__ inline float d_protected_exp(float a) { 
    if (a > 20.0f) return 4.85e8f; if (a < -20.0f) return 0.0f; return __expf(a); 
}
__device__ inline float d_protected_ln(float a) { 
    float v = fabsf(a); return (v < 1e-6f) ? -13.8f : logf(v); 
}

// ------------------------------------------------------------------------------------------------
//  KERNEL: Evaluate Population with Shared Memory Caching
// ------------------------------------------------------------------------------------------------
// Grid: Blocks = Population Size.
// Threads: Process data rows in parallel.
// Shared Memory: Stores the bytecode for the block's individual.

__global__ void kernel_evaluate_vm(
    const Instruction* __restrict__ all_codes, 
    const int* __restrict__ offsets, 
    const int* __restrict__ lengths,
    const float* __restrict__ soa_features, 
    const float* __restrict__ L_true,
    const float* __restrict__ L_base,
    float* __restrict__ out_errors, // Returns raw absolute errors
    int num_rows,
    int soa_stride
) {
    // 1. Identify Individual (One Block per Individual)
    int ind_idx = blockIdx.x;
    
    // 2. Load Program into Shared Memory
    // The shared memory size is defined dynamically at launch
    extern __shared__ Instruction s_code[];
    
    int my_offset = offsets[ind_idx];
    int my_len = lengths[ind_idx];

    // Cooperative loading: threads load instructions in parallel
    for (int i = threadIdx.x; i < my_len; i += blockDim.x) {
        s_code[i] = all_codes[my_offset + i];
    }
    
    // Barrier to ensure program is fully loaded in L1 cache
    __syncthreads();

    // 3. Process Data Rows (Grid Stride Loop)
    // Even if rows > blockDim, the loop handles it.
    for (int r = threadIdx.x; r < num_rows; r += blockDim.x) {
        
        // VM Stack (Local Registers)
        float stack[MAX_STACK_SIZE];
        int sp = 0;
        bool error_flag = false;

        // Execute Bytecode from L1 Shared Memory
        for (int pc = 0; pc < my_len; ++pc) {
            const Instruction instr = s_code[pc];
            
            switch(instr.op) {
                case OP_VAR: 
                    // Coalesced Read from SoA Global Memory
                    stack[sp++] = soa_features[instr.idx * soa_stride + r]; 
                    break;
                case OP_CONST: 
                    stack[sp++] = instr.val; 
                    break;
                case OP_ADD: sp--; stack[sp-1] += stack[sp]; break;
                case OP_MUL: sp--; stack[sp-1] *= stack[sp]; break;
                case OP_SUB: sp--; stack[sp-1] -= stack[sp]; break;
                case OP_DIV: 
                    sp--; 
                    stack[sp-1] = d_protected_div(stack[sp-1], stack[sp]); 
                    break;
                case OP_SQRT: stack[sp-1] = d_protected_sqrt(stack[sp-1]); break;
                case OP_SQ:   stack[sp-1] = d_protected_sq(stack[sp-1]); break;
                case OP_TANH: stack[sp-1] = d_protected_tanh(stack[sp-1]); break;
                case OP_EXP:  stack[sp-1] = d_protected_exp(stack[sp-1]); break;
                case OP_LN:   stack[sp-1] = d_protected_ln(stack[sp-1]); break;
            }
            
            // Stack overflow check
            if (sp >= MAX_STACK_SIZE) { error_flag = true; break; }
        }

        float result = (sp > 0) ? stack[0] : 0.0f;
        
        // Error handling & Difference Calculation
        if (error_flag || isnan(result) || isinf(result) || fabsf(result) > 100.0f) {
            // Signal failure (-1.0)
            out_errors[ind_idx * num_rows + r] = -1.0f; 
        } else {
            // Calculate Absolute Error: abs(True - (Base * Multiplier))
            float pred = L_base[r] * result;
            out_errors[ind_idx * num_rows + r] = fabsf(L_true[r] - pred);
        }
    }
}

// ================================================================================================
//  SECTION 3: CPU HOST CLASSES & GEP ENGINE
// ================================================================================================

struct CompiledProgram {
    Instruction code[MAX_PROGRAM_SIZE];
    int length;
};

struct Node { int id; double value; bool is_int; };
struct Gene { std::vector<Node> sequence; }; 

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
            for (int i = 0; i < GENE_LENGTH; i++) {
                int min_id = (i < HEAD_LENGTH) ? 0 : 9;
                
                // 30% Chance for Constant (ID 18)
                if(random_double(0, 1) < 0.3) {
                    genes[g].sequence[i].id = 18; 
                } else {
                    genes[g].sequence[i].id = random_int(min_id, 17);
                }
                
                genes[g].sequence[i].is_int = false; 
                
                if (genes[g].sequence[i].id == 18) {
                    if (random_double(0, 1) < 0.5) { 
                        genes[g].sequence[i].value = (double)random_int(-10, 10);
                        genes[g].sequence[i].is_int = true;
                    } else {
                         if (random_double(0, 1) < 0.3) {
                            double magics[] = { 0.0, 1.0, 2.0, 0.5, PI, G, 2*PI };
                            genes[g].sequence[i].value = magics[random_int(0, 6)];
                        } else {
                            genes[g].sequence[i].value = random_double(-5, 5);
                        }
                    }
                }
            }
        }
    }

    // Compiles Gene to Bytecode (Host Side)
    void compile() {
        program.length = 0;
        int prog_idx = 0;

        for (int g = 0; g < N_GENES; g++) {
            struct TreeNode { int id; double val; TreeNode *l=0, *r=0; };
            TreeNode pool[256]; int ptr = 0; 
            TreeNode* q[256]; int head=0, tail=0;
            
            TreeNode* root = &pool[ptr++];
            root->id = genes[g].sequence[0].id; root->val = genes[g].sequence[0].value;
            q[tail++] = root;
            
            int seq_idx = 1;
            while(head < tail) {
                TreeNode* curr = q[head++];
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

            std::function<void(TreeNode*)> emit = [&](TreeNode* n) {
                if(!n) return;
                emit(n->l);
                emit(n->r);
                if(prog_idx >= MAX_PROGRAM_SIZE - 2) return; 
                Instruction instr;
                if(n->id >= 9) {
                    if(n->id == 18) { instr.op = OP_CONST; instr.val = (float)n->val; } // Cast to float
                    else { instr.op = OP_VAR; instr.idx = (uint8_t)(n->id - 9); }
                } else {
                    switch(n->id) {
                        case 0: instr.op = OP_ADD; break; case 1: instr.op = OP_SUB; break;
                        case 2: instr.op = OP_MUL; break; case 3: instr.op = OP_DIV; break;
                        case 4: instr.op = OP_SQRT; break; case 5: instr.op = OP_SQ; break;
                        case 6: instr.op = OP_TANH; break; case 7: instr.op = OP_EXP; break;
                        case 8: instr.op = OP_LN; break;
                    }
                }
                program.code[prog_idx++] = instr;
            };
            emit(root);
            if (g > 0 && prog_idx < MAX_PROGRAM_SIZE - 1) {
                Instruction link; link.op = OP_ADD; program.code[prog_idx++] = link;
            }
        }
        program.length = prog_idx;
    }
};

// ================================================================================================
//  SECTION 4: SYMBOLIC MATH & REPORTING (CPU HOST)
// ================================================================================================

std::string to_fraction(double val, bool is_int = false) {
    if (is_int || (std::abs(val - std::round(val)) < 1e-9)) return std::to_string((int)std::round(val));
    if (std::abs(val) < 1e-9) return "0";
    double v = std::abs(val); int sign = (val < 0) ? -1 : 1;
    int best_n = 1, best_d = 1; double min_err = 1e9;
    for (int d = 1; d <= 1000; ++d) {
        int n = (int)std::round(v * d); double err = std::abs(v - (double)n / d);
        if (err < min_err) { min_err = err; best_n = n; best_d = d; }
        if (err < 1e-9) break;
    }
    if (best_d == 1) return std::to_string(sign * best_n);
    return std::to_string(sign * best_n) + "/" + std::to_string(best_d);
}

struct SymNode {
    int id; double val; bool is_int;
    SymNode *l = nullptr, *r = nullptr; double power = 1.0; 
    bool is_const() const { return id == 18; }
    bool equals(const SymNode* other) const {
        if (!other) return false; if (id != other->id) return false;
        if (is_const()) return std::abs(val - other->val) < 1e-6;
        return (l && other->l ? l->equals(other->l) : l==other->l) && (r && other->r ? r->equals(other->r) : r==other->r);
    }
    bool is_non_negative() const { return (id==5 || id==4 || id==7) || (is_const() && val >= 0); }
};

struct SymArena {
    std::vector<SymNode*> nodes;
    ~SymArena() { for(auto n : nodes) delete n; }
    SymNode* create(int id, double val = 0.0, bool is_int = false) {
        SymNode* n = new SymNode{id, val, is_int}; nodes.push_back(n); return n;
    }
};

SymNode* simplify_ast(SymNode* n, SymArena& arena) {
    if (!n) return nullptr; if (n->id >= 9) return n; 
    n->l = simplify_ast(n->l, arena); n->r = simplify_ast(n->r, arena);
    if ((n->l && n->l->is_const()) && (!n->r || (n->r && n->r->is_const()))) {
        double l_val = n->l->val; double r_val = (n->r) ? n->r->val : 0.0; double res = 0.0; bool ok = true;
        switch(n->id) {
            case 0: res = l_val + r_val; break; case 1: res = l_val - r_val; break;
            case 2: res = l_val * r_val; break; case 3: if(std::abs(r_val)>1e-9) res=l_val/r_val; else ok=false; break;
            case 4: res = std::sqrt(std::abs(l_val)); break; case 5: res = l_val*l_val; break;
            case 6: res = std::tanh(l_val); break; case 7: res = std::exp(l_val); break;
            case 8: res = std::log(std::abs(l_val)+1e-9); break;
        }
        if (ok) { n->id = 18; n->val = res; n->is_int = (std::abs(res - std::round(res)) < 1e-9); n->l = 0; n->r = 0; return n; }
    }
    return n;
}

int get_precedence(int id) { if(id==0||id==1)return 1; if(id==2||id==3)return 2; if(id==5)return 3; return 4; }

std::string print_ast(SymNode* n, int parent_prec = 0, double accumulated_pow = 1.0) {
    if (!n) return "";
    double current_pow = accumulated_pow * n->power;
    if (n->id == 4) { 
        if (n->l && (n->l->id == 4 || n->l->id == 5)) return print_ast(n->l, parent_prec, current_pow * 0.5);
        double final_p = current_pow * 0.5; std::string inner = print_ast(n->l, 0, 1.0);
        if (!n->l->is_non_negative()) inner = "abs(" + inner + ")";
        if (std::abs(final_p - 0.5) < 1e-5) return "sqrt(" + inner + ")";
        return "(" + inner + "**" + to_fraction(final_p) + ")";
    }
    if (n->id == 18) return to_fraction(n->val, n->is_int);
    if (n->id >= 9) return "x" + std::to_string(n->id - 9);
    if (n->id == 6) return "tanh(" + print_ast(n->l, 0, 1.0) + ")";
    if (n->id == 7) return "exp(" + print_ast(n->l, 0, 1.0) + ")";
    if (n->id == 8) return "ln(" + print_ast(n->l, 0, 1.0) + ")";
    
    int my_prec = get_precedence(n->id); bool parens = (my_prec < parent_prec); std::string s;
    switch(n->id) {
        case 0: s = print_ast(n->l, my_prec) + " + " + print_ast(n->r, my_prec); break;
        case 1: s = print_ast(n->l, my_prec) + " - " + print_ast(n->r, my_prec + 1); break;
        case 2: s = print_ast(n->l, my_prec) + " * " + print_ast(n->r, my_prec); break;
        case 3: s = print_ast(n->l, my_prec) + " / " + print_ast(n->r, my_prec + 1); break;
    }
    return parens ? "(" + s + ")" : s;
}

std::string get_formula_str(const Individual& ind) {
    std::string s = ""; SymArena arena;
    for(int i=0; i<N_GENES; ++i) {
        if(i > 0) s += " + ";
        const Gene& gene = ind.genes[i];
        std::vector<SymNode*> q; 
        SymNode* root = arena.create(gene.sequence[0].id, gene.sequence[0].value, gene.sequence[0].is_int);
        q.push_back(root); int idx = 1, cur = 0;
        while(cur < q.size() && idx < GENE_LENGTH) {
            SymNode* node = q[cur++];
            if(node->id <= 8) {
                int arity = (node->id <= 3) ? 2 : 1;
                for(int k=0; k<arity; ++k) {
                    if(idx>=GENE_LENGTH) break;
                    SymNode* c = arena.create(gene.sequence[idx].id, gene.sequence[idx].value, gene.sequence[idx].is_int);
                    idx++; if(k==0) node->l=c; else node->r=c; q.push_back(c);
                }
            }
        }
        simplify_ast(root, arena); simplify_ast(root, arena); s += print_ast(root);
    }
    return s;
}

struct DetailedStats {
    double mape, max_err, rmse, bias, p50, p90, p99;
    std::vector<double> residuals, pct_err;
    struct Regime { double mape, max; int count; };
    Regime shallow, inter, deep;
};

DetailedStats calculate_stats(const std::vector<double>& true_val, const std::vector<double>& pred_val, 
                              const std::vector<double>& safe_val, const std::vector<double>& d_arr) {
    DetailedStats s; size_t n = true_val.size();
    s.residuals.resize(n); s.pct_err.resize(n);
    double sum_pct=0, sum_sq=0, sum_res=0;
    for(size_t i=0; i<n; ++i) {
        double p = pred_val[i]; if (host_safe_isnan(p) || host_safe_isinf(p)) p = 0.0; 
        s.residuals[i] = true_val[i] - p;
        s.pct_err[i] = (std::abs(s.residuals[i]) / std::abs(safe_val[i])) * 100.0;
        sum_pct += s.pct_err[i]; sum_sq += s.residuals[i]*s.residuals[i]; sum_res += s.residuals[i];
    }
    s.mape = sum_pct/n; s.rmse = std::sqrt(sum_sq/n); s.bias = sum_res/n;
    std::vector<double> srtd = s.pct_err; std::sort(srtd.begin(), srtd.end());
    s.max_err = srtd.back(); s.p50 = srtd[size_t(0.5*n)]; s.p90 = srtd[size_t(0.9*n)]; 
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
    std::vector<double> t = v; 
    for(auto& val : t) if(host_safe_isnan(val) || host_safe_isinf(val)) val = 1e9;
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
    auto pr = [&](std::string n, auto r){ ss<<"   "<<std::left<<std::setw(15)<<n<<" | "<<std::setw(6)<<r.count<<" | "<<std::setw(9)<<r.mape<<"% | "<<r.max<<"%\n"; };
    pr("Shallow",s.shallow); pr("Intermediate",s.inter); pr("Deep",s.deep);
    ss<<std::string(80,'='); logger.println(ss.str());
}

// ==============================================================================
//  SECTION 5: EVOLUTION OPERATORS & HALL OF FAME
// ==============================================================================

class HallOfFame {
public:
    std::vector<Individual> members; size_t capacity;
    HallOfFame(size_t cap) : capacity(cap) { members.reserve(cap); }
    void update(const std::vector<Individual>& pop) {
        size_t candidates_count = std::min(pop.size(), capacity * 2); 
        std::vector<Individual> candidates = members;
        candidates.insert(candidates.end(), pop.begin(), pop.begin() + candidates_count);
        std::sort(candidates.begin(), candidates.end(), [](const Individual& a, const Individual& b){
            bool a_nan = host_safe_isnan(a.fitness); bool b_nan = host_safe_isnan(b.fitness);
            if (a_nan && b_nan) return false; if (a_nan) return false; if (b_nan) return true;
            return a.fitness < b.fitness;
        });
        auto is_duplicate = [](const Individual& a, const Individual& b) {
            return (std::abs(a.fitness - b.fitness) <= 1e-9) && (a.program.length == b.program.length);
        };
        auto last = std::unique(candidates.begin(), candidates.end(), is_duplicate);
        candidates.erase(last, candidates.end());
        if (candidates.size() > capacity) candidates.resize(capacity);
        members = candidates;
    }
    Individual& best() { return members[0]; }
};

void crossover_one_point(Individual& p1, Individual& p2) {
    if (random_double(0, 1) < 0.4) {
        int total = N_GENES * GENE_LENGTH; int pt = random_int(1, total - 1);
        for (int i = pt; i < total; i++) {
            int g = i / GENE_LENGTH, n = i % GENE_LENGTH;
            std::swap(p1.genes[g].sequence[n], p2.genes[g].sequence[n]);
        }
    }
}
void crossover_two_point(Individual& p1, Individual& p2) {
    if (random_double(0, 1) < 0.2) {
        int total = N_GENES * GENE_LENGTH; int pt1 = random_int(0, total - 2); int pt2 = random_int(pt1 + 1, total - 1);
        for(int i = pt1; i <= pt2; ++i) {
            int g = i / GENE_LENGTH, n = i % GENE_LENGTH;
            std::swap(p1.genes[g].sequence[n], p2.genes[g].sequence[n]);
        }
    }
}
void crossover_gene(Individual& p1, Individual& p2) {
    if (N_GENES > 1 && random_double(0, 1) < 0.1) {
        int g = random_int(0, N_GENES - 1); std::swap(p1.genes[g], p2.genes[g]);
    }
}
void mutate_uniform(Individual& ind) {
    const double MUTATION_RATE = 0.044;
    for (int g = 0; g < N_GENES; g++) {
        for (int i = 0; i < GENE_LENGTH; i++) {
            if (random_double(0, 1) < MUTATION_RATE) {
                int min_id = (i < HEAD_LENGTH) ? 0 : 9;
                if (random_double(0, 1) < 0.3) ind.genes[g].sequence[i].id = 18;
                else ind.genes[g].sequence[i].id = random_int(min_id, 17);
                ind.genes[g].sequence[i].is_int = false;
                if (ind.genes[g].sequence[i].id == 18) {
                     if(random_double(0, 1) < 0.5) { ind.genes[g].sequence[i].value = (double)random_int(-10, 10); ind.genes[g].sequence[i].is_int = true; }
                     else { ind.genes[g].sequence[i].value += random_double(-1.0, 1.0); }
                }
            }
        }
    }
}
void mutate_gaussian(Individual& ind) {
    for (int g = 0; g < N_GENES; g++) {
        for (int i = 0; i < GENE_LENGTH; i++) {
            if (ind.genes[g].sequence[i].id == 18) {
                if (random_double(0, 1) < 0.1) {
                    if (ind.genes[g].sequence[i].is_int) { int step = (random_double(0, 1) < 0.5) ? 1 : -1; ind.genes[g].sequence[i].value += step; }
                    else { double perturbation = random_double(-0.05, 0.05); ind.genes[g].sequence[i].value *= (1.0 + perturbation); }
                }
            }
        }
    }
}
void mutate_is_transpose(Individual& ind) {
    if (random_double(0, 1) < 0.1) {
        int g = random_int(0, N_GENES - 1); int max_len = (HEAD_LENGTH > 3) ? (HEAD_LENGTH / 3) : 1; 
        int len = random_int(1, std::max(1, max_len)); int src = random_int(0, GENE_LENGTH - len); int tgt = random_int(0, HEAD_LENGTH - 1); 
        std::vector<Node> chunk(len); for(int i=0; i<len; ++i) chunk[i] = ind.genes[g].sequence[src+i];
        std::vector<Node> new_head = ind.genes[g].sequence; for(int i=0; i<len; ++i) if(tgt+i < HEAD_LENGTH) new_head[tgt+i] = chunk[i];
        int r = tgt; int w = tgt + len; while(w < HEAD_LENGTH) new_head[w++] = ind.genes[g].sequence[r++];
        for(int i=0; i<HEAD_LENGTH; ++i) ind.genes[g].sequence[i] = new_head[i];
    }
}
void mutate_ris_transpose(Individual& ind) {
    if (random_double(0, 1) < 0.1) {
        int g = random_int(0, N_GENES - 1); std::vector<int> func_indices;
        for(int i=0; i<HEAD_LENGTH; ++i) if(ind.genes[g].sequence[i].id <= 8) func_indices.push_back(i);
        if(func_indices.empty()) return;
        int start_idx = func_indices[random_int(0, func_indices.size()-1)];
        int max_len = (HEAD_LENGTH > 3) ? (HEAD_LENGTH / 3) : 1; int len = random_int(1, std::max(1, max_len));
        if(start_idx + len > GENE_LENGTH) len = GENE_LENGTH - start_idx;
        std::vector<Node> ris(len); for(int i=0; i<len; ++i) ris[i] = ind.genes[g].sequence[start_idx+i];
        std::vector<Node> new_head(HEAD_LENGTH); for(int i=0; i<len; ++i) if(i<HEAD_LENGTH) new_head[i] = ris[i];
        int r = 0; int w = len; while(w < HEAD_LENGTH) new_head[w++] = ind.genes[g].sequence[r++];
        for(int i=0; i<HEAD_LENGTH; ++i) ind.genes[g].sequence[i] = new_head[i];
    }
}
void mutate_gene_transpose(Individual& ind) {
    if (N_GENES > 1 && random_double(0, 1) < 0.1) {
        int src_g = random_int(0, N_GENES - 1); Gene target_gene = ind.genes[src_g];
        for(int g = src_g; g > 0; --g) ind.genes[g] = ind.genes[g-1];
        ind.genes[0] = target_gene;
    }
}
void mutate_invert(Individual& ind) {
    if (random_double(0, 1) < 0.1) {
        int g = random_int(0, N_GENES - 1); int p1 = random_int(0, HEAD_LENGTH - 1), p2 = random_int(0, HEAD_LENGTH - 1);
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

    // Gene Config
    std::cout << "\n--- GENE CONFIGURATION ---\n";
    std::cout << "Enter Number of Genes: "; std::cin >> N_GENES;
    std::cout << "Enter Head Length: "; std::cin >> HEAD_LENGTH;
    GENE_LENGTH = HEAD_LENGTH + (HEAD_LENGTH * (MAX_ARITY - 1) + 1);
    
    // Load Data
    std::ifstream f(FILENAME);
    if(!f.good()) return 1;
    std::vector<double> H, T_raw, d_arr, U, L_true;
    std::string line; std::getline(f, line);
    while(std::getline(f, line)) {
        std::stringstream ss(line);
        double h, t, dv, u, l;
        if(ss >> h >> t >> dv >> u >> l) {
            H.push_back(h); T_raw.push_back(t); d_arr.push_back(dv); U.push_back(u); L_true.push_back(l);
        }
    }
    size_t rows = L_true.size();
    if(rows == 0) return 1;
    logger.println("Loaded " + std::to_string(rows) + " rows.");

    // Physics Baseline
    std::vector<double> L_base = solve_linear_no_current(T_raw, d_arr);
    
    // Build SoA (Float) for GPU
    // This converts the 9 features calculated on CPU into a flat float array
    std::vector<double> soa_flat = build_soa_features_float(H, T_raw, d_arr, U, L_base);
    
    // Convert Host Doubles to Float for Device
    std::vector<float> soa_f(soa_flat.begin(), soa_flat.end());
    std::vector<float> L_true_f(L_true.begin(), L_true.end());
    std::vector<float> L_base_f(L_base.begin(), L_base.end());
    
    std::vector<double> L_safe = L_true; for(auto& v : L_safe) if(v==0) v=1e-6;

    // ----------------------------------------------------------------
    // GPU MEMORY ALLOCATION (PERSISTENT)
    // ----------------------------------------------------------------
    float *d_soa, *d_L_true, *d_L_base, *d_errors;
    Instruction *d_codes; 
    int *d_offsets, *d_lengths;

    size_t code_buf_size = POPULATION_SIZE * MAX_PROGRAM_SIZE * sizeof(Instruction);
    
    checkCuda(cudaMalloc(&d_soa, soa_f.size() * sizeof(float)));
    checkCuda(cudaMalloc(&d_L_true, rows * sizeof(float)));
    checkCuda(cudaMalloc(&d_L_base, rows * sizeof(float)));
    checkCuda(cudaMalloc(&d_errors, POPULATION_SIZE * rows * sizeof(float)));
    checkCuda(cudaMalloc(&d_codes, code_buf_size));
    checkCuda(cudaMalloc(&d_offsets, POPULATION_SIZE * sizeof(int)));
    checkCuda(cudaMalloc(&d_lengths, POPULATION_SIZE * sizeof(int)));

    // Copy Static Data Once
    checkCuda(cudaMemcpy(d_soa, soa_f.data(), soa_f.size() * sizeof(float), cudaMemcpyHostToDevice));
    checkCuda(cudaMemcpy(d_L_true, L_true_f.data(), rows * sizeof(float), cudaMemcpyHostToDevice));
    checkCuda(cudaMemcpy(d_L_base, L_base_f.data(), rows * sizeof(float), cudaMemcpyHostToDevice));

    // Population
    std::vector<Individual> pop(POPULATION_SIZE);
    std::vector<Individual> next(POPULATION_SIZE);
    for(auto& ind : pop) { ind.randomize(); ind.compile(); }
    
    HallOfFame hof(POPULATION_SIZE / 10);
    double global_best_score = DBL_MAX;

    std::vector<Instruction> host_code_flat; host_code_flat.reserve(POPULATION_SIZE * 50);
    std::vector<int> host_offsets(POPULATION_SIZE);
    std::vector<int> host_lengths(POPULATION_SIZE);
    std::vector<float> host_errors(POPULATION_SIZE * rows);

    logger.println("Starting Evolution (" + std::to_string(N_GENERATIONS) + " gens) on GPU... Press Ctrl+C to stop.\n");

    for(int gen = 1; gen <= N_GENERATIONS; ++gen) {
        
        // 1. Flatten Population
        host_code_flat.clear();
        for(int i=0; i<POPULATION_SIZE; ++i) {
            host_offsets[i] = host_code_flat.size();
            host_lengths[i] = pop[i].program.length;
            for(int k=0; k<pop[i].program.length; ++k) 
                host_code_flat.push_back(pop[i].program.code[k]);
        }

        // 2. Transfer Code to GPU
        checkCuda(cudaMemcpy(d_codes, host_code_flat.data(), host_code_flat.size() * sizeof(Instruction), cudaMemcpyHostToDevice));
        checkCuda(cudaMemcpy(d_offsets, host_offsets.data(), POPULATION_SIZE * sizeof(int), cudaMemcpyHostToDevice));
        checkCuda(cudaMemcpy(d_lengths, host_lengths.data(), POPULATION_SIZE * sizeof(int), cudaMemcpyHostToDevice));

        // 3. Launch Kernel
        // 1 Block per Individual. Threads per block = 256. 
        // Shared Mem size = MAX_PROGRAM_SIZE * sizeof(Instruction).
        int threads = 256;
        int blocks = POPULATION_SIZE;
        size_t shared_mem = MAX_PROGRAM_SIZE * sizeof(Instruction); 

        kernel_evaluate_vm<<<blocks, threads, shared_mem>>>(
            d_codes, d_offsets, d_lengths, 
            d_soa, d_L_true, d_L_base, d_errors, 
            (int)rows, (int)rows
        );
        checkCuda(cudaGetLastError());

        // 4. Retrieve Errors
        checkCuda(cudaMemcpy(host_errors.data(), d_errors, host_errors.size() * sizeof(float), cudaMemcpyDeviceToHost));

        // 5. Aggregate Stats (Host Reduction)
        // We do reduction on CPU to support the complex metrics (P99, Median) easily.
        // Transfer overhead is minimal compared to CPU evaluation time.
        #pragma omp parallel for schedule(static)
        for(int i=0; i<POPULATION_SIZE; ++i) {
            float* start = &host_errors[i * rows];
            double sum_pct = 0.0;
            double max_pct = 0.0;
            double sum_sq = 0.0;
            std::vector<double> p99_buf; 
            bool need_p99 = (choice == 4 || choice == 6 || choice == 7);
            if(need_p99) p99_buf.reserve(rows);

            bool failed = false;

            for(int r=0; r<rows; ++r) {
                float diff = start[r]; // This is abs(True - Pred) from GPU
                if(diff < 0.0f) { failed = true; break; }

                double safe = std::abs(L_safe[r]);
                double pct = (double)diff / safe * 100.0;
                
                if (choice == 1) sum_pct += pct;
                else if (choice == 2) max_pct = std::max(max_pct, pct);
                else if (choice == 3) sum_sq += (double)diff * diff;
                else {
                    if (need_p99) p99_buf.push_back(pct);
                    if (choice >= 5) {
                        sum_sq += (double)diff * diff;
                        max_pct = std::max(max_pct, pct);
                    }
                }
            }

            double penalty = pop[i].program.length * 0.05;

            if(failed) pop[i].fitness = 1e9;
            else {
                if(choice == 1) pop[i].fitness = (sum_pct / rows) + penalty;
                else if(choice == 2) pop[i].fitness = max_pct + penalty;
                else if(choice == 3) pop[i].fitness = std::sqrt(sum_sq/rows) + penalty;
                else if(choice == 4) pop[i].fitness = calc_p99(p99_buf) + penalty;
                else if(choice == 6) pop[i].fitness = (calc_p99(p99_buf) * max_pct) + penalty;
                else if(choice == 7) pop[i].fitness = (std::sqrt(sum_sq/rows) * calc_p99(p99_buf) * max_pct) + penalty;
                else if(choice == 5) pop[i].fitness = (std::sqrt(sum_sq/rows) * max_pct) + penalty;
            }
        }

        // 6. Evolution Logic
        for(auto& ind : pop) {
            if (host_safe_isnan(ind.fitness) || host_safe_isinf(ind.fitness) || ind.fitness < 0) ind.fitness = 1e9;
        }

        int elites = 5;
        // Partial sort to find elites
        std::partial_sort(pop.begin(), pop.begin() + std::max(elites, POPULATION_SIZE/10), pop.end(), 
            [](const Individual& a, const Individual& b){ return a.fitness < b.fitness; });
        
        hof.update(pop);

        if(hof.best().fitness < global_best_score - 1e-5 && hof.best().fitness < 1e8) {
            global_best_score = hof.best().fitness;
            
            // Reconstruct predictions on CPU for full reporting accuracy
            // We use the double-precision CPU math here to generate the report
            std::vector<double> preds(rows); 
            for(size_t r=0; r<rows; ++r) {
                // Manually execute VM on CPU for report
                double stack[64]; int sp=0;
                const Instruction* code = hof.best().program.code;
                int len = hof.best().program.length;
                for(int pc=0; pc<len; ++pc) {
                     Instruction in = code[pc];
                     switch(in.op) {
                         // Note: soa_flat is organized as [Feature0...Feature1...]
                         case OP_VAR: stack[sp++] = soa_flat[in.idx*rows + r]; break;
                         case OP_CONST: stack[sp++] = (double)in.val; break;
                         case OP_ADD: sp--; stack[sp-1]+=stack[sp]; break;
                         case OP_MUL: sp--; stack[sp-1]*=stack[sp]; break;
                         case OP_SUB: sp--; stack[sp-1]-=stack[sp]; break;
                         case OP_DIV: sp--; stack[sp-1] = fabsf(stack[sp])<1e-9?NAN:stack[sp-1]/stack[sp]; break;
                         case OP_SQRT: stack[sp-1] = sqrt(fabs(stack[sp-1])); break;
                         case OP_SQ: stack[sp-1] = stack[sp-1]*stack[sp-1]; break;
                         case OP_TANH: stack[sp-1] = tanh(stack[sp-1]); break;
                         case OP_EXP: stack[sp-1] = (stack[sp-1]>20)?4.85e8:exp(stack[sp-1]); break;
                         case OP_LN: stack[sp-1] = log(fabs(stack[sp-1])+1e-9); break;
                     }
                }
                double res = (sp>0)?stack[0]:0.0;
                if(host_safe_isnan(res)) res=0.0;
                preds[r] = L_base[r] * res;
            }
            print_full_report(calculate_stats(L_true, preds, L_safe, d_arr), get_formula_str(hof.best()), gen);
            logger.flush(); 
        }

        int next_idx = 0;
        for(int k=0; k<elites; ++k) next[next_idx++] = pop[k];
        
        while(next_idx < POPULATION_SIZE) {
            auto tourn = [&](){ 
                int b = random_int(0, POPULATION_SIZE-1);
                for(int k=0; k<2; k++) { 
                    int c = random_int(0, POPULATION_SIZE-1); 
                    if(pop[c].fitness < pop[b].fitness) b = c; 
                }
                return pop[b];
            };

            Individual c1 = tourn(), c2 = tourn();
            crossover_one_point(c1, c2); crossover_two_point(c1, c2); crossover_gene(c1, c2); 
            mutate_is_transpose(c1); mutate_ris_transpose(c1); mutate_gene_transpose(c1); mutate_invert(c1);
            mutate_uniform(c1); mutate_gaussian(c1);      
            c1.compile(); next[next_idx++] = c1;
            
            if(next_idx < POPULATION_SIZE) {
                mutate_is_transpose(c2); mutate_ris_transpose(c2); mutate_gene_transpose(c2); mutate_invert(c2);
                mutate_uniform(c2); mutate_gaussian(c2);
                c2.compile(); next[next_idx++] = c2;
            }
        }
        pop.swap(next);
    }
    
    // Final Report
    if(hof.best().fitness < 1e8) {
        Individual& best_ind = hof.best();
        std::vector<double> preds(rows);
        for(size_t r=0; r<rows; ++r) {
            double stack[64]; int sp=0;
            const Instruction* code = best_ind.program.code;
            int len = best_ind.program.length;
            for(int pc=0; pc<len; ++pc) {
                 Instruction in = code[pc];
                 switch(in.op) {
                     case OP_VAR: stack[sp++] = soa_flat[in.idx*rows + r]; break;
                     case OP_CONST: stack[sp++] = (double)in.val; break;
                     case OP_ADD: sp--; stack[sp-1]+=stack[sp]; break;
                     case OP_MUL: sp--; stack[sp-1]*=stack[sp]; break;
                     case OP_SUB: sp--; stack[sp-1]-=stack[sp]; break;
                     case OP_DIV: sp--; stack[sp-1] = fabsf(stack[sp])<1e-9?NAN:stack[sp-1]/stack[sp]; break;
                     case OP_SQRT: stack[sp-1] = sqrt(fabs(stack[sp-1])); break;
                     case OP_SQ: stack[sp-1] = stack[sp-1]*stack[sp-1]; break;
                     case OP_TANH: stack[sp-1] = tanh(stack[sp-1]); break;
                     case OP_EXP: stack[sp-1] = (stack[sp-1]>20)?4.85e8:exp(stack[sp-1]); break;
                     case OP_LN: stack[sp-1] = log(fabs(stack[sp-1])+1e-9); break;
                 }
            }
            double res = (sp>0)?stack[0]:0.0;
            if(host_safe_isnan(res)) res=0.0;
            preds[r] = L_base[r] * res;
        }
        DetailedStats s = calculate_stats(L_true, preds, L_safe, d_arr);
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
                << std::setw(8) << d_arr[idx] << std::setw(8) << U[idx] << std::setw(10) << L_true[idx] << std::setw(10) << preds[idx] 
                << std::setw(10) << s.residuals[idx] << std::setw(8) << s.pct_err[idx] << "%";
            logger.println(rss.str());
        }
    } else {
        logger.println("\nNo valid model found (all individuals failed constraints).");
    }
    logger.println("\n" + std::string(80, '='));
    logger.println("L_Final = L_linear_no_current * (Multiplier)");
    logger.println(std::string(80, '='));
    
    // Cleanup
    cudaFree(d_soa); cudaFree(d_L_true); cudaFree(d_L_base);
    cudaFree(d_errors); cudaFree(d_codes); cudaFree(d_offsets); cudaFree(d_lengths);
    return 0;
}