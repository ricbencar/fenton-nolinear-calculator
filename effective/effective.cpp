/*
====================================================================================================
|                                     EVOLUTIONARY FENTON SOLVER                                   |
|                                     (DUAL-FUNCTION EVOLUTION)                                    |
====================================================================================================

PROBLEM DEFINITION:
----------------------------------------------------------------------------------------------------
Target: Predict Wavelength (L) by evolving TWO simultaneous functions:
  1. Effective Period: T_eff = g(H, T, U, d)
  2. Effective Depth:  d_eff = h(H, T, U, d)

The Solver then computes: 
  L_pred = Linear_Dispersion_Relation(T_eff, d_eff)

INPUTS MAPPING (Unified for both functions):
----------------------------------------------------------------------------------------------------
  Var 0: H (Wave Height)
  Var 1: T (Wave Period)
  Var 2: U (Current Velocity)
  Var 3: d (Water Depth)

COMPILATION (Windows MinGW64):
----------------------------------------------------------------------------------------------------
g++ -O3 -std=c++17 -static -fopenmp -ffast-math -funroll-loops ^
-march=native -fstack-protector-strong -Wl,--stack,67108864 ^
-Wl,--allow-multiple-definition effective.cpp -o effective.exe
====================================================================================================
*/

#ifdef __MINGW32__
#define NOMINMAX
#include <stdio.h>
#include <time.h>
#include <windows.h>
extern "C" {
    #pragma GCC diagnostic ignored "-Wattributes"
    int fseeko64(FILE* stream, long long offset, int origin) { return _fseeki64(stream, offset, origin); }
    long long ftello64(FILE* stream) { return _ftelli64(stream); }
    int (*__imp_fseeko64)(FILE*, long long, int) = fseeko64;
    long long (*__imp_ftello64)(FILE*) = ftello64;
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
#include <omp.h>
#include <chrono>
#include <limits>
#include <cstring>

// ------------------------------------------------------------------------------------------------
//                                     GLOBAL CONFIGURATION
// ------------------------------------------------------------------------------------------------

const std::string FILENAME = "list.txt";
const std::string OUTPUT_FILE = "output.txt";
const int POPULATION_SIZE = 1000;
const int N_GENERATIONS = 100000;

// GEP Gene Configuration 
int GENES_PER_FUNC = 2; // Genes per function (Total genes = 2 * this)
int HEAD_LENGTH = 8;   
int GENE_LENGTH = 0;    // Calculated at runtime

const int MAX_ARITY = 2; 
const int MAX_STACK_SIZE = 64;
const int MAX_PROGRAM_SIZE = 2048; 
const int INPUT_VARS = 4; // H, T, U, d

// Physics Constants
const double PI = 3.14159265359;
const double G = 9.80665;

// ------------------------------------------------------------------------------------------------
//                                     UTILITIES
// ------------------------------------------------------------------------------------------------

class DualLogger {
    std::ofstream file;
public:
    DualLogger(const std::string& filename) { file.open(filename); }
    template <typename T> void print(const T& msg) { std::cout << msg; if (file.is_open()) file << msg; }
    template <typename T> void println(const T& msg) { std::cout << msg << std::endl; if (file.is_open()) file << msg << std::endl; }
    void flush() { std::cout.flush(); if (file.is_open()) file.flush(); }
};
DualLogger logger(OUTPUT_FILE);

inline bool safe_isnan(double d) {
    uint64_t u; std::memcpy(&u, &d, sizeof(d));
    return (u & 0x7ff0000000000000ULL) == 0x7ff0000000000000ULL && (u & 0x000fffffffffffffULL);
}
inline bool safe_isinf(double d) {
    uint64_t u; std::memcpy(&u, &d, sizeof(d));
    return (u & 0x7fffffffffffffffULL) == 0x7ff0000000000000ULL;
}

class XorShift256 {
    uint64_t s[4];
    static uint64_t splitmix64(uint64_t& x) {
        uint64_t z = (x += 0x9e3779b97f4a7c15);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
        z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
        return z ^ (z >> 31);
    }
public:
    XorShift256(uint64_t seed = 0) {
        if(seed == 0) seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        s[0] = splitmix64(seed); s[1] = splitmix64(seed); s[2] = splitmix64(seed); s[3] = splitmix64(seed);
    }
    uint64_t next() {
        const uint64_t result = s[0] + s[3];
        const uint64_t t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t; s[3] = (s[3] << 45) | (s[3] >> 19);
        return result;
    }
    double next_double() { return (next() >> 11) * 0x1.0p-53; }
    int next_int(int min, int max) { return min + (next() % (max - min + 1)); }
};
XorShift256 rng;
double random_double(double min, double max) { return min + (rng.next_double() * (max - min)); }
int random_int(int min, int max) { return rng.next_int(min, max); }

// ================================================================================================
//  SECTION 1: PHYSICS KERNEL
// ================================================================================================

inline double protected_div(double a, double b) { return std::abs(b) < 1e-9 ? std::numeric_limits<double>::quiet_NaN() : a / b; }
inline double protected_sqrt(double a) { return std::sqrt(std::abs(a)); }
inline double protected_sq(double a) { return (std::abs(a) > 1e5) ? 1e10 : a * a; }
inline double protected_tanh(double a) { return std::tanh(a); }
inline double protected_exp(double a) { if (a > 20.0) return 4.85e8; if (a < -20.0) return 0.0; return std::exp(a); }
inline double protected_ln(double a) { double v = std::abs(a); return (v < 1e-6) ? -13.8 : std::log(v); }

// SOLVER: L = Linear(T_eff, d_eff)
inline double solve_wavelength_dual(double T, double d) {
    if (std::abs(T) < 1e-4) return 0.0; 
    
    // Physical constraints: T and d must be positive
    double T_abs = std::abs(T);
    double d_abs = std::abs(d);
    if(T_abs < 0.1) T_abs = 0.1;
    if(d_abs < 0.01) d_abs = 0.01;

    double omega = (2.0 * PI) / T_abs;
    double k = (omega * omega) / G; // Deep water guess
    if (k < 1e-6) k = 1e-6;

    // Newton-Raphson
    for (int j = 0; j < 5; j++) {
        if (k > 500.0) k = 500.0; 
        double kd = k * d_abs;
        double th = std::tanh(kd);
        double f = G * k * th - (omega * omega);
        double df = G * th + G * kd * (1.0 - th * th);
        
        if (std::abs(df) < 1e-9) break;
        double diff = f / df;
        k = k - diff;
        if (std::abs(diff) < 1e-5) break;
    }
    return (2.0 * PI) / k;
}

// ================================================================================================
//  SECTION 2: JIT BYTECODE ENGINE
// ================================================================================================

enum OpCode : uint8_t {
    OP_ADD, OP_SUB, OP_MUL, OP_DIV, OP_SQRT, OP_SQ, OP_TANH, OP_EXP, OP_LN,
    OP_VAR, OP_CONST, OP_END
};

struct Instruction { uint8_t op; uint8_t idx; double val; };
struct CompiledProgram { Instruction code[MAX_PROGRAM_SIZE]; int length; };
struct Node { int id; double value; bool is_int; };
struct Gene { std::vector<Node> sequence; }; 

class Individual {
public:
    std::vector<Gene> genes; // Holds 2 * GENES_PER_FUNC
    double fitness;
    
    // TWO compiled programs per individual
    CompiledProgram prog_T; 
    CompiledProgram prog_d;

    Individual() : fitness(DBL_MAX) { 
        prog_T.length = 0; prog_d.length = 0;
        genes.resize(GENES_PER_FUNC * 2);
        for(auto& g : genes) g.sequence.resize(GENE_LENGTH);
    }

    void randomize() {
        for (auto& g : genes) {
            for (int i = 0; i < GENE_LENGTH; i++) {
                int min_id = (i < HEAD_LENGTH) ? 0 : 9;
                if(random_double(0, 1) < 0.3) {
                    g.sequence[i].id = 18; 
                } else {
                    g.sequence[i].id = random_int(min_id, 8 + INPUT_VARS); 
                }
                g.sequence[i].is_int = false; 
                if (g.sequence[i].id == 18) {
                    if (random_double(0, 1) < 0.5) { 
                        g.sequence[i].value = (double)random_int(-10, 10);
                        g.sequence[i].is_int = true;
                    } else {
                         if (random_double(0, 1) < 0.3) {
                            double magics[] = { 0.0, 1.0, 2.0, 0.5, PI, G, 2*PI };
                            g.sequence[i].value = magics[random_int(0, 6)];
                        } else {
                            g.sequence[i].value = random_double(-5, 5);
                        }
                    }
                }
            }
        }
    }

    // Compiles a range of genes into a target program
    void compile_range(int start_gene, int end_gene, CompiledProgram& out_prog) {
        out_prog.length = 0;
        int prog_idx = 0;
        for (int g = start_gene; g < end_gene; g++) {
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
                emit(n->l); emit(n->r);
                if(prog_idx >= MAX_PROGRAM_SIZE - 2) return; 
                Instruction instr;
                if(n->id >= 9) {
                    if(n->id == 18) { instr.op = OP_CONST; instr.val = n->val; }
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
                out_prog.code[prog_idx++] = instr;
            };
            emit(root);
            if (g > start_gene && prog_idx < MAX_PROGRAM_SIZE - 1) {
                Instruction link; link.op = OP_ADD; out_prog.code[prog_idx++] = link;
            }
        }
        out_prog.length = prog_idx;
    }

    void compile() {
        compile_range(0, GENES_PER_FUNC, prog_T);
        compile_range(GENES_PER_FUNC, GENES_PER_FUNC * 2, prog_d);
    }
};

inline double execute_vm(const CompiledProgram& prog, const double* vars) {
    double stack[MAX_STACK_SIZE];
    int sp = 0;
    for(int i = 0; i < prog.length; ++i) {
        const Instruction& instr = prog.code[i];
        switch(instr.op) {
            case OP_VAR: stack[sp++] = vars[instr.idx]; break;
            case OP_CONST: stack[sp++] = instr.val; break;
            case OP_ADD: sp--; stack[sp-1] += stack[sp]; break;
            case OP_MUL: sp--; stack[sp-1] *= stack[sp]; break;
            case OP_SUB: sp--; stack[sp-1] -= stack[sp]; break;
            case OP_DIV: sp--; stack[sp-1] = protected_div(stack[sp-1], stack[sp]); break;
            case OP_SQRT: stack[sp-1] = protected_sqrt(stack[sp-1]); break;
            case OP_SQ:   stack[sp-1] = protected_sq(stack[sp-1]); break;
            case OP_TANH: stack[sp-1] = protected_tanh(stack[sp-1]); break;
            case OP_EXP:  stack[sp-1] = protected_exp(stack[sp-1]); break;
            case OP_LN:   stack[sp-1] = protected_ln(stack[sp-1]); break;
        }
    }
    return (sp > 0) ? stack[0] : 0.0;
}

// ================================================================================================
//  SECTION 3: SYMBOLIC MATH ENGINE (UNIFIED 4-VAR SUPPORT)
// ================================================================================================

std::string to_fraction(double val, bool is_int = false) {
    if (is_int || (std::abs(val - std::round(val)) < 1e-9)) return std::to_string((int)std::round(val));
    if (std::abs(val) < 1e-9) return "0";
    double v = std::abs(val); int sign = (val < 0) ? -1 : 1;
    double min_err = 1e9; int best_n = 1, best_d = 1;
    for (int d = 1; d <= 1000; ++d) {
        int n = (int)std::round(v * d);
        double err = std::abs(v - (double)n / d);
        if (err < min_err) { min_err = err; best_n = n; best_d = d; }
        if (err < 1e-9) break;
    }
    if (best_d == 1) return std::to_string(sign * best_n);
    return std::to_string(sign * best_n) + "/" + std::to_string(best_d);
}

int get_precedence(int id) {
    if (id == 0 || id == 1) return 1; if (id == 2 || id == 3) return 2; if (id == 5) return 3; return 4;
}

struct SymNode {
    int id; double val; bool is_int; SymNode *l = nullptr, *r = nullptr; double power = 1.0; 
    bool is_const() const { return id == 18; }
};

struct SymArena {
    std::vector<SymNode*> nodes;
    ~SymArena() { for(auto n : nodes) delete n; }
    SymNode* create(int id, double val = 0.0, bool is_int = false) {
        SymNode* n = new SymNode{id, val, is_int}; nodes.push_back(n); return n;
    }
};

SymNode* simplify_ast(SymNode* n, SymArena& arena) {
    if (!n) return nullptr;
    if (n->id >= 9) return n; 
    n->l = simplify_ast(n->l, arena); n->r = simplify_ast(n->r, arena);

    if ((n->l && n->l->is_const()) && (!n->r || (n->r && n->r->is_const()))) {
        double lv = n->l->val; double rv = (n->r) ? n->r->val : 0.0;
        double res = 0.0; bool ok = true;
        switch(n->id) {
            case 0: res = lv + rv; break; case 1: res = lv - rv; break;
            case 2: res = lv * rv; break; case 3: if(std::abs(rv) > 1e-9) res = lv/rv; else ok = false; break;
            case 4: res = protected_sqrt(lv); break; case 5: res = protected_sq(lv); break;
            case 6: res = protected_tanh(lv); break; case 7: res = protected_exp(lv); break;
            case 8: res = protected_ln(lv); break;
        }
        if (ok) { n->id = 18; n->val = res; n->is_int = (std::abs(res - std::round(res)) < 1e-9); n->l = 0; n->r = 0; return n; }
    }
    if (n->id == 0) { if (n->r && n->r->is_const() && std::abs(n->r->val) < 1e-9) return n->l; if (n->l && n->l->is_const() && std::abs(n->l->val) < 1e-9) return n->r; }
    if (n->id == 2) { if (n->r && n->r->is_const() && std::abs(n->r->val - 1.0) < 1e-9) return n->l; if (n->l && n->l->is_const() && std::abs(n->l->val - 1.0) < 1e-9) return n->r; }
    return n;
}

std::string print_ast(SymNode* n, int parent_prec = 0, double accumulated_pow = 1.0) {
    if (!n) return "";
    double current_pow = accumulated_pow * n->power;
    if (n->id == 4) { 
        if (n->l && (n->l->id == 4 || n->l->id == 5)) return print_ast(n->l, parent_prec, current_pow * 0.5);
        double fp = current_pow * 0.5; std::string in = print_ast(n->l, 0, 1.0);
        if (std::abs(fp - 0.5) < 1e-5) return "sqrt(" + in + ")";
        if (std::abs(fp - 1.0) < 1e-5) return in;
        return "(" + in + "**" + to_fraction(fp) + ")";
    }
    if (n->id == 5) {
        if (n->l && (n->l->id == 4 || n->l->id == 5)) return print_ast(n->l, parent_prec, current_pow * 2.0);
        double fp = current_pow * 2.0; std::string in = print_ast(n->l, 3, 1.0);
        if (std::abs(fp - 1.0) < 1e-5) return in;
        return "(" + in + "**" + to_fraction(fp) + ")";
    }
    if (n->id == 18) return to_fraction(n->val, n->is_int);
    
    // UNIFIED VARIABLE MAPPING: 0=H, 1=T, 2=U, 3=d
    if (n->id >= 9) {
        int idx = n->id - 9;
        if(idx == 0) return "H";
        if(idx == 1) return "T";
        if(idx == 2) return "U";
        if(idx == 3) return "d";
        return "?";
    }
    
    if (n->id == 6) return "tanh(" + print_ast(n->l, 0, 1.0) + ")";
    if (n->id == 7) return "exp(" + print_ast(n->l, 0, 1.0) + ")";
    if (n->id == 8) return "ln(" + print_ast(n->l, 0, 1.0) + ")";

    int mp = get_precedence(n->id); bool par = (mp < parent_prec); 
    std::string s;
    switch(n->id) {
        case 0: s = print_ast(n->l, mp) + " + " + print_ast(n->r, mp); break;
        case 1: s = print_ast(n->l, mp) + " - " + print_ast(n->r, mp + 1); break;
        case 2: s = print_ast(n->l, mp) + " * " + print_ast(n->r, mp); break;
        case 3: s = print_ast(n->l, mp) + " / " + print_ast(n->r, mp + 1); break;
    }
    return par ? "(" + s + ")" : s;
}

std::string simplify_gene_str(const Gene& gene) {
    SymArena arena; std::vector<SymNode*> q; 
    SymNode* root = arena.create(gene.sequence[0].id, gene.sequence[0].value, gene.sequence[0].is_int); q.push_back(root);
    int idx = 1, cur = 0;
    while(cur < q.size() && idx < GENE_LENGTH) {
        SymNode* node = q[cur++];
        if(node->id <= 8) {
            int arity = (node->id <= 3) ? 2 : 1;
            for(int k=0; k<arity; ++k) {
                if(idx >= GENE_LENGTH) break;
                SymNode* c = arena.create(gene.sequence[idx].id, gene.sequence[idx].value, gene.sequence[idx].is_int);
                idx++; if(k==0) node->l = c; else node->r = c; q.push_back(c);
            }
        }
    }
    simplify_ast(root, arena); simplify_ast(root, arena);
    return print_ast(root);
}

// ==============================================================================
//  SECTION 4: EVOLUTION LOGIC & REPORTING
// ==============================================================================

struct DetailedStats {
    double mape, max_err, rmse, bias, p50, p90, p99;
    std::vector<double> residuals, pct_err;
    struct Regime { double mape, max; int count; };
    Regime shallow, inter, deep;
};

DetailedStats calculate_stats(const std::vector<double>& true_val, const std::vector<double>& pred_val, 
                              const std::vector<double>& d_arr) {
    DetailedStats s;
    size_t n = true_val.size();
    s.residuals.resize(n); s.pct_err.resize(n);
    double sum_pct=0, sum_sq=0, sum_res=0;
    
    for(size_t i=0; i<n; ++i) {
        double p = pred_val[i];
        if (safe_isnan(p) || safe_isinf(p)) p = 0.0; 
        
        s.residuals[i] = true_val[i] - p;
        s.pct_err[i] = (std::abs(s.residuals[i]) / std::abs(true_val[i])) * 100.0;
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
            double rel = d_arr[i]/true_val[i];
            if(rel>=l && rel<h) { r_sum+=s.pct_err[i]; if(s.pct_err[i]>r_max)r_max=s.pct_err[i]; c++; }
        }
        return { c?r_sum/c:0, r_max, c };
    };
    s.shallow = reg(0,0.05); s.inter = reg(0.05,0.5); s.deep = reg(0.5,1e9);
    return s;
}

double calc_p99(const std::vector<double>& v) {
    std::vector<double> t = v; 
    for(auto& val : t) if(safe_isnan(val) || safe_isinf(val)) val = 1e9;
    size_t idx = size_t(std::min(size_t(0.99*t.size()), t.size()-1));
    std::nth_element(t.begin(), t.begin() + idx, t.end());
    return t[idx];
}

void print_full_report(const DetailedStats& s, const Individual& ind, int gen=-1) {
    std::stringstream ss;
    if(gen!=-1) ss<<"\n["<<std::setw(6)<<gen<<"] >>> NEW BEST FOUND <<<\n";
    ss<<std::string(80,'-')<<"\n";
    
    // Print T_eff
    ss << "T_eff(H,T,U,d) = ";
    for(int i=0; i<GENES_PER_FUNC; ++i) {
        if(i>0) ss << " + ";
        ss << simplify_gene_str(ind.genes[i]);
    }
    ss << "\n";
    
    // Print d_eff
    ss << "d_eff(H,T,U,d) = ";
    for(int i=GENES_PER_FUNC; i<GENES_PER_FUNC*2; ++i) {
        if(i>GENES_PER_FUNC) ss << " + ";
        ss << simplify_gene_str(ind.genes[i]);
    }
    ss << "\n";
    
    ss<<std::string(80,'-')<<"\n";
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

// Hall of Fame
class HallOfFame {
public:
    std::vector<Individual> members;
    size_t capacity;
    HallOfFame(size_t cap) : capacity(cap) { members.reserve(cap); }
    void update(const std::vector<Individual>& pop) {
        size_t candidates_count = std::min(pop.size(), capacity * 2); 
        std::vector<Individual> candidates = members;
        candidates.insert(candidates.end(), pop.begin(), pop.begin() + candidates_count);
        std::sort(candidates.begin(), candidates.end(), [](const Individual& a, const Individual& b){
            bool a_nan = safe_isnan(a.fitness); bool b_nan = safe_isnan(b.fitness);
            if (a_nan && b_nan) return false; if (a_nan) return false; if (b_nan) return true;
            return a.fitness < b.fitness;
        });
        auto is_duplicate = [](const Individual& a, const Individual& b) {
            if (std::abs(a.fitness - b.fitness) > 1e-9) return false;
            return a.prog_T.length == b.prog_T.length && a.prog_d.length == b.prog_d.length;
        };
        auto last = std::unique(candidates.begin(), candidates.end(), is_duplicate);
        candidates.erase(last, candidates.end());
        if (candidates.size() > capacity) candidates.resize(capacity);
        members = candidates;
    }
    Individual& best() { return members[0]; }
};

// ------------------------------------------------------------------------------------------------
//  SECTION 5: FULL SUITE OF GENETIC OPERATORS
// ------------------------------------------------------------------------------------------------

// 1. One-Point Crossover (40%)
void crossover_one_point(Individual& p1, Individual& p2) {
    if (random_double(0, 1) < 0.4) {
        int total_genes = GENES_PER_FUNC * 2 * GENE_LENGTH;
        int pt = random_int(1, total_genes - 1);
        for (int i = pt; i < total_genes; i++) {
            int g = i / GENE_LENGTH, n = i % GENE_LENGTH;
            std::swap(p1.genes[g].sequence[n], p2.genes[g].sequence[n]);
        }
    }
}

// 2. Two-Point Crossover (20%)
void crossover_two_point(Individual& p1, Individual& p2) {
    if (random_double(0, 1) < 0.2) {
        int total = GENES_PER_FUNC * 2 * GENE_LENGTH;
        int pt1 = random_int(0, total - 2);
        int pt2 = random_int(pt1 + 1, total - 1);
        for(int i = pt1; i <= pt2; ++i) {
            int g = i / GENE_LENGTH, n = i % GENE_LENGTH;
            std::swap(p1.genes[g].sequence[n], p2.genes[g].sequence[n]);
        }
    }
}

// 3. Gene Crossover (10%)
void crossover_gene(Individual& p1, Individual& p2) {
    if (random_double(0, 1) < 0.1) {
        int g = random_int(0, (GENES_PER_FUNC * 2) - 1);
        std::swap(p1.genes[g], p2.genes[g]);
    }
}

// 4. Uniform Mutation (Rate varies)
void mutate_uniform(Individual& ind) {
    for (auto& gene : ind.genes) {
        for (int i = 0; i < GENE_LENGTH; i++) {
            if (random_double(0, 1) < 0.044) {
                int min_id = (i < HEAD_LENGTH) ? 0 : 9;
                if (random_double(0, 1) < 0.3) gene.sequence[i].id = 18;
                else gene.sequence[i].id = random_int(min_id, 8 + INPUT_VARS);
                gene.sequence[i].is_int = false;
                if (gene.sequence[i].id == 18) {
                    if(random_double(0, 1) < 0.5) {
                         gene.sequence[i].value = (double)random_int(-10, 10);
                         gene.sequence[i].is_int = true;
                    } else {
                        gene.sequence[i].value += random_double(-1.0, 1.0);
                    }
                }
            }
        }
    }
}

// 5. Gaussian / Step Mutation (Fine Tuning)
void mutate_gaussian(Individual& ind) {
    for (auto& gene : ind.genes) {
        for (int i = 0; i < GENE_LENGTH; i++) {
            if (gene.sequence[i].id == 18) {
                if (random_double(0, 1) < 0.1) {
                    if (gene.sequence[i].is_int) {
                        int step = (random_double(0, 1) < 0.5) ? 1 : -1;
                        gene.sequence[i].value += step;
                    } else {
                        double val = gene.sequence[i].value;
                        gene.sequence[i].value = val + (val * random_double(-0.05, 0.05)); 
                    }
                }
            }
        }
    }
}

// 6. IS Transposition
void mutate_is_transpose(Individual& ind) {
    if (random_double(0, 1) < 0.1) {
        int g = random_int(0, (GENES_PER_FUNC * 2) - 1);
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

// 7. RIS Transposition
void mutate_ris_transpose(Individual& ind) {
    if (random_double(0, 1) < 0.1) {
        int g = random_int(0, (GENES_PER_FUNC * 2) - 1);
        std::vector<int> func_indices;
        for(int i=0; i<HEAD_LENGTH; ++i) if(ind.genes[g].sequence[i].id <= 8) func_indices.push_back(i);
        if(func_indices.empty()) return;
        int start_idx = func_indices[random_int(0, func_indices.size()-1)];
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

// 8. Gene Transposition
void mutate_gene_transpose(Individual& ind) {
    if (random_double(0, 1) < 0.1) {
        int src_g = random_int(0, (GENES_PER_FUNC * 2) - 1);
        Gene target_gene = ind.genes[src_g];
        for(int g = src_g; g > 0; --g) ind.genes[g] = ind.genes[g-1];
        ind.genes[0] = target_gene;
    }
}

// 9. Inversion
void mutate_invert(Individual& ind) {
    if (random_double(0, 1) < 0.1) {
        int g = random_int(0, (GENES_PER_FUNC * 2) - 1);
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
    logger.println("\n--- FENTON SOLVER CONFIGURATION (DUAL-FUNCTION V5) ---");
    logger.println("Optimizing T_eff(H,T,U,d) and d_eff(H,T,U,d)");
    
    // Metric Selection (Full suite)
    logger.println("\nSelect Optimization Metric:");
    logger.println("  1) MAPE");
    logger.println("  2) MAX ERROR");
    logger.println("  3) RMSE");
    logger.println("  4) P99 ERROR");
    logger.println("  5) RMSE * MAX ERROR");
    logger.println("  6) P99 ERROR * MAX ERROR");
    logger.println("  7) RMSE * P99 ERROR * MAX ERROR");
    
    int choice = 1;
    std::string in;
    while(true) {
        std::cout << "Enter choice (1-7): ";
        std::getline(std::cin, in);
        if(!in.empty()){ try{choice=stoi(in);}catch(...){choice=0;} if(choice>=1 && choice<=7) break;}
        std::cout << "Invalid selection.\n";
    }
    
    // Config
    std::cout << "\nEnter GENES_PER_FUNC (e.g. 1): "; while(!(std::cin >> GENES_PER_FUNC) || GENES_PER_FUNC<=0);
    std::cout << "Enter HEAD_LENGTH (e.g. 8): "; while(!(std::cin >> HEAD_LENGTH) || HEAD_LENGTH<=0);
    GENE_LENGTH = HEAD_LENGTH + (HEAD_LENGTH * (MAX_ARITY - 1) + 1);
    
    // Data Loading
    std::ifstream f(FILENAME);
    if(!f.good()) { logger.println("Error: File not found."); return 1; }
    std::vector<double> H_raw, T_raw, d_raw, U_raw, L_true;
    std::string line; std::getline(f, line);
    while(std::getline(f, line)) {
        std::stringstream ss(line);
        double h, t, dv, u, l;
        if(ss >> h >> t >> dv >> u >> l) {
            H_raw.push_back(h); T_raw.push_back(t); d_raw.push_back(dv); U_raw.push_back(u); L_true.push_back(l);
        }
    }
    size_t rows = L_true.size();
    if(rows == 0) return 1;
    logger.print("Loaded "); logger.print(rows); logger.println(" rows.");

    // Population Init
    std::vector<Individual> pop(POPULATION_SIZE);
    for(auto& ind : pop) { ind.randomize(); ind.compile(); }
    std::vector<Individual> next(POPULATION_SIZE);
    HallOfFame hof(POPULATION_SIZE / 10);
    double global_best_score = DBL_MAX;

    // Evolution
    for(int gen = 1; gen <= N_GENERATIONS; ++gen) {
        
        #pragma omp parallel for schedule(static)
        for(int i = 0; i < POPULATION_SIZE; ++i) {
            double vars[4]; // 0=H, 1=T, 2=U, 3=d
            double fit=0, max_e=0, sq_e=0;
            std::vector<double> p99_vec; if(choice>=4) p99_vec.resize(rows);
            bool err = false;
            
            for(size_t r = 0; r < rows; ++r) {
                // UNIFIED MAPPING
                vars[0] = H_raw[r];
                vars[1] = T_raw[r];
                vars[2] = U_raw[r];
                vars[3] = d_raw[r];
                
                double T_eff = execute_vm(pop[i].prog_T, vars);
                double d_eff = execute_vm(pop[i].prog_d, vars);
                
                // Solve
                double L_pred = solve_wavelength_dual(T_eff, d_eff);

                if(safe_isnan(L_pred) || safe_isinf(L_pred)) { err = true; break; }
                
                double diff = std::abs(L_true[r] - L_pred);
                double pct = (L_true[r] > 1e-6) ? diff/std::abs(L_true[r])*100.0 : 0.0;
                
                if(choice==1) fit += pct;
                else if(choice==2) max_e = std::max(max_e, pct);
                else if(choice==3) sq_e += diff*diff;
                else if(choice==4) p99_vec[r] = pct;
                else if(choice==6) { p99_vec[r] = pct; max_e = std::max(max_e, pct); }
                else if(choice==7) { p99_vec[r] = pct; max_e = std::max(max_e, pct); sq_e += diff*diff; }
                else { sq_e += diff*diff; max_e = std::max(max_e, pct); }
            }
            
            double penalty = (pop[i].prog_T.length + pop[i].prog_d.length) * 0.01;
            if(err) pop[i].fitness = 1e9;
            else if(choice==1) pop[i].fitness = (fit/rows) + penalty;
            else if(choice==2) pop[i].fitness = max_e + penalty;
            else if(choice==3) pop[i].fitness = std::sqrt(sq_e/rows) + penalty;
            else if(choice==4) pop[i].fitness = calc_p99(p99_vec) + penalty;
            else if(choice==6) pop[i].fitness = (calc_p99(p99_vec) * max_e) + penalty;
            else if(choice==7) pop[i].fitness = (std::sqrt(sq_e/rows) * calc_p99(p99_vec) * max_e) + penalty;
            else pop[i].fitness = (std::sqrt(sq_e/rows) * max_e) + penalty;
        }

        // Selection & Elitism
        for(auto& ind : pop) if(safe_isnan(ind.fitness) || ind.fitness < 0) ind.fitness = 1e9;
        std::partial_sort(pop.begin(), pop.begin() + 10, pop.end(), [](auto& a, auto& b){ return a.fitness < b.fitness; });
        
        hof.update(pop);

        if(hof.best().fitness < global_best_score - 1e-5 && hof.best().fitness < 1e8) {
            global_best_score = hof.best().fitness;
            std::vector<double> preds(rows); double v[4];
            for(size_t r=0; r<rows; ++r) {
                v[0]=H_raw[r]; v[1]=T_raw[r]; v[2]=U_raw[r]; v[3]=d_raw[r];
                double te = execute_vm(hof.best().prog_T, v);
                double de = execute_vm(hof.best().prog_d, v);
                preds[r] = solve_wavelength_dual(te, de);
            }
            print_full_report(calculate_stats(L_true, preds, d_raw), hof.best(), gen);
            logger.flush();
        }

        int elites = 5;
        for(int k=0; k<elites; ++k) next[k] = pop[k];
        int next_idx = elites;
        while(next_idx < POPULATION_SIZE) {
            int b = random_int(0, POPULATION_SIZE-1);
            for(int k=0; k<2; k++) { int c = random_int(0, POPULATION_SIZE-1); if(pop[c].fitness < pop[b].fitness) b = c; }
            Individual c1 = pop[b];

            b = random_int(0, POPULATION_SIZE-1);
            for(int k=0; k<2; k++) { int c = random_int(0, POPULATION_SIZE-1); if(pop[c].fitness < pop[b].fitness) b = c; }
            Individual c2 = pop[b];

            crossover_one_point(c1, c2);
            crossover_two_point(c1, c2);
            crossover_gene(c1, c2);

            mutate_is_transpose(c1);
            mutate_ris_transpose(c1);
            mutate_gene_transpose(c1);
            mutate_invert(c1);
            mutate_uniform(c1); 
            mutate_gaussian(c1);
            c1.compile();
            next[next_idx++] = c1;

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
        pop = next;
    }
    
    // Final Report
    if(hof.best().fitness < 1e8) {
        Individual& best_ind = hof.best();
        std::vector<double> preds(rows); double v[4];
        for(size_t r=0; r<rows; ++r) {
            v[0]=H_raw[r]; v[1]=T_raw[r]; v[2]=U_raw[r]; v[3]=d_raw[r];
            double te = execute_vm(best_ind.prog_T, v);
            double de = execute_vm(best_ind.prog_d, v);
            preds[r] = solve_wavelength_dual(te, de);
        }
        DetailedStats s = calculate_stats(L_true, preds, d_raw);
        logger.println("\n" + std::string(80, '='));
        logger.println("                        FINAL ANALYSIS REPORT (HALL OF FAME BEST)");
        logger.println(std::string(80, '='));
        print_full_report(s, best_ind);
        
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
            rss << "   " << std::left << std::setw(6) << idx << std::setw(8) << H_raw[idx] << std::setw(8) << T_raw[idx] 
                << std::setw(8) << d_raw[idx] << std::setw(8) << U_raw[idx] << std::setw(10) << L_true[idx] << std::setw(10) << preds[idx] 
                << std::setw(10) << s.residuals[idx] << std::setw(8) << s.pct_err[idx] << "%";
            logger.println(rss.str());
        }
    } else {
        logger.println("\nNo valid model found.");
    }

    return 0;
}