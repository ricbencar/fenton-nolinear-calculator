# -*- coding: utf-8 -*-
"""
pade_nn.py
==========
Single (global) neural-network surrogate for wavelength L.

Core workflow is intentionally kept the same as the previous Padé script:
- Reads list.txt (H, T, d, Uc, L).
- Builds a physics baseline L_lin from linear dispersion with current.
- Builds a bank of 5 physics-inspired candidate features from (H, T, d, Uc, L_lin).
- For each modelling strategy (x/y transforms + optional tail reweighting):
    - Selects the best 3 features via Spearman relevance + mRMR (relevance vs redundancy).
    - Fits a small neural network on those 3 features.
    - Scans DEGREE_GRID and selects best by validation MAPE.
- Retrains the best configuration on all data.
- Generates:
    function.py (self-contained: includes all NN coefficients + preprocessing; no external files)
    diagnostics.pdf
    plot_all.png

Neural network choice (for speed + determinism):
- ELM-style MLP: random hidden layer (tanh) + trained linear output (ridge).
  This is still a neural network (nonlinear hidden units), but training reduces to a fast
  weighted ridge solve. This keeps the "fit/score/select" loop practical even with many
  strategies and sizes.

Input:
    list.txt with columns: H T d Uc L
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from numba import njit

# ==============================================================================
# CONSTANTS / CONFIG
# ==============================================================================
G = 9.80665
PI = 3.141592653589793

# Kept name for continuity with the previous Padé script: now used as NN size grid.
DEGREE_GRID = [10, 25, 50]

# Hidden units = HIDDEN_MULT * degree
HIDDEN_MULT = 10

LAMBDA_RIDGE = 1e-6
TRAIN_FRACTION = 0.80
MAX_TRAIN_SAMPLES = None  # inferred from list.txt at runtime

Y_CLIP = 2.5
REWEIGHT_BETA = 5.0
WEIGHT_CLIP_MIN = 1e-3
WEIGHT_CLIP_MAX = 100.0

RNG_SEED = 12345

FEATURE_NAMES = [
    "WaveSteepness (H/L_lin)",
    "RelativeDepth (d/L_lin)",
    "DopplerFactor (Uc*T/L_lin)",
    "UrsellNumber (H*L_lin^2/d^3)",
    "CurrentFroude (Uc/sqrt(g d))",
]

# ==============================================================================
# IO
# ==============================================================================
def load_data(path: str | Path = "list.txt") -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Missing input file: {path.resolve()}")

    try:
        df = pd.read_csv(path, sep="\t")
        if df.shape[1] < 5:
            df = pd.read_csv(path, sep=r"\s+")
    except Exception:
        df = pd.read_csv(path, sep=r"\s+")

    if "s" in df.columns and "d" not in df.columns:
        df = df.rename(columns={"s": "d"})

    required = {"H", "T", "d", "Uc", "L"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"list.txt missing required columns: {sorted(missing)}")

    return df


def infer_max_train_samples(df: pd.DataFrame) -> int:
    """Infer the effective MAX_TRAIN_SAMPLES directly from list.txt.

    The training cap is taken as the number of valid rows that were actually
    loaded from list.txt. This preserves the legacy variable name while making
    the limit data-driven rather than hard-coded.
    """
    n_rows = int(len(df))
    if n_rows <= 0:
        raise ValueError("list.txt contains no valid data rows.")
    return n_rows


# ==============================================================================
# PHYSICS BASELINE (linear dispersion + current)
# ==============================================================================
@njit(cache=True, fastmath=True)
def _solve_k_no_current(omega: float, d: float) -> float:
    k_deep = (omega * omega) / G
    c_shallow = math.sqrt(G * max(d, 1e-12))
    k_shallow = omega / max(c_shallow, 1e-12)
    k = k_deep if (k_deep * d) > 1.0 else k_shallow

    for _ in range(8):
        k = max(1e-9, min(k, 400.0))
        kd = k * d
        th = math.tanh(kd)
        sigma = math.sqrt(max(G * k * th, 1e-16))
        f = sigma - omega
        sech2 = 1.0 - th * th
        if sigma > 1e-12:
            d_sigma = (G * th + G * kd * sech2) / (2.0 * sigma)
        else:
            d_sigma = 0.0
        if abs(d_sigma) < 1e-12:
            break
        k_new = k - f / d_sigma
        if k_new <= 1e-9:
            k_new = 1e-9
        k = 0.7 * k + 0.3 * k_new

    return max(k, 1e-9)


@njit(cache=True, fastmath=True)
def _solve_k_with_current(omega: float, d: float, Uc: float) -> float:
    k = _solve_k_no_current(omega, d)

    for _ in range(40):
        k = max(1e-9, min(k, 400.0))
        kd = k * d
        th = math.tanh(kd)
        sigma = math.sqrt(max(G * k * th, 1e-16))
        f = sigma + k * Uc - omega

        sech2 = 1.0 - th * th
        if sigma > 1e-12:
            d_sigma = (G * th + G * kd * sech2) / (2.0 * sigma)
        else:
            d_sigma = 0.0

        df = d_sigma + Uc
        if abs(df) < 1e-12:
            break

        k_new = k - f / df
        if k_new <= 1e-9:
            k_new = 1e-9

        k = 0.8 * k + 0.2 * k_new
        if abs(k_new - k) / max(k, 1e-9) < 1e-10:
            break

    return max(k, 1e-9)


@njit(cache=True, fastmath=True)
def solve_linear_wavelength(T_arr: np.ndarray, d_arr: np.ndarray, U_arr: np.ndarray) -> np.ndarray:
    n = len(T_arr)
    L = np.empty(n, dtype=np.float64)
    for i in range(n):
        omega = 2.0 * PI / T_arr[i]
        k = _solve_k_with_current(omega, d_arr[i], U_arr[i])
        L[i] = 2.0 * PI / k
    return L


# ==============================================================================
# STRATEGIES (x/y transforms + robust reweighting)
# ==============================================================================
@dataclass(frozen=True)
class Strategy:
    name: str
    y_mode: str
    x_mode: str
    lam: float = LAMBDA_RIDGE
    reweight_mode: str = "mape"   # none | mape | tail | p95 | topq | huber_down
    reweight_beta: float = REWEIGHT_BETA
    reweight_gamma: float = 1.0
    reweight_iters: int = 1


STRATEGIES = [
    Strategy("A_logratio_loggroups", "log_ratio", "log_groups"),
    Strategy("B_ratiominus1_rawgroups", "ratio_minus1", "raw_groups"),

    Strategy("C_logratio_rawgroups", "log_ratio", "raw_groups"),
    Strategy("D_ratiominus1_loggroups", "ratio_minus1", "log_groups"),

    Strategy("E_asinhLratio_loggroups", "asinh_ratio_minus1", "log_groups"),
    Strategy("F_asinhLratio_logasinh", "asinh_ratio_minus1", "log_asinh_groups"),
    Strategy("G_signedlogLratio_signedlog", "signedlog_ratio_minus1", "signedlog_groups"),
    Strategy("H_signedlogLratio_rawsignedlog", "signedlog_ratio_minus1", "raw_signedlog_groups"),

    Strategy("I_kratio_rawgroups", "k_ratio_minus1", "raw_groups"),
    Strategy("J_kratio_loggroups", "k_ratio_minus1", "log_groups"),
    Strategy("K_asinhKratio_logasinh", "asinh_k_ratio_minus1", "log_asinh_groups"),
    Strategy("L_signedlogKratio_signedlog", "signedlog_k_ratio_minus1", "signedlog_groups"),

    Strategy("M_logratio_log_tail", "log_ratio", "log_groups",
             lam=1e-6, reweight_mode="tail", reweight_beta=8.0, reweight_gamma=2.0, reweight_iters=2),
    Strategy("N_kratio_log_p95", "k_ratio_minus1", "log_groups",
             lam=1e-6, reweight_mode="p95", reweight_beta=6.0, reweight_gamma=2.0, reweight_iters=2),
]


def build_feature_bank(
    H: np.ndarray, T: np.ndarray, d: np.ndarray, Uc: np.ndarray, L_lin: np.ndarray, x_mode: str
) -> np.ndarray:
    eps = 1e-12
    Ls = np.maximum(L_lin, eps)
    ds = np.maximum(d, eps)

    S = H / Ls
    R = d / Ls
    D = Uc * T / Ls
    Ur = H * (L_lin * L_lin) / (ds * ds * ds)
    Fr = Uc / np.sqrt(G * ds)

    def _signed_log1p(z: np.ndarray) -> np.ndarray:
        return np.sign(z) * np.log1p(np.abs(z))

    if x_mode == "log_groups":
        return np.column_stack([
            np.log(np.maximum(S, eps)),
            np.log(np.maximum(R, eps)),
            D,
            np.log(np.maximum(Ur, eps)),
            Fr,
        ])

    if x_mode == "raw_groups":
        return np.column_stack([S, R, D, Ur, Fr])

    if x_mode == "asinh_groups":
        return np.column_stack([np.arcsinh(S), np.arcsinh(R), np.arcsinh(D), np.arcsinh(Ur), np.arcsinh(Fr)])

    if x_mode == "log_asinh_groups":
        return np.column_stack([
            np.log(np.maximum(S, eps)),
            np.log(np.maximum(R, eps)),
            np.arcsinh(D),
            np.log(np.maximum(Ur, eps)),
            np.arcsinh(Fr),
        ])

    if x_mode == "signedlog_groups":
        return np.column_stack([
            np.log(np.maximum(S, eps)),
            np.log(np.maximum(R, eps)),
            _signed_log1p(D),
            np.log(np.maximum(Ur, eps)),
            _signed_log1p(Fr),
        ])

    if x_mode == "raw_signedlog_groups":
        return np.column_stack([S, R, _signed_log1p(D), Ur, _signed_log1p(Fr)])

    raise ValueError(f"Unknown x_mode: {x_mode}")


def build_target(L_true: np.ndarray, L_lin: np.ndarray, y_mode: str) -> np.ndarray:
    eps = 1e-12
    Ls = np.maximum(L_lin, eps)
    Lt = np.maximum(L_true, eps)

    if y_mode == "log_ratio":
        return np.log(Lt / Ls)

    if y_mode == "ratio_minus1":
        return (Lt / Ls) - 1.0

    if y_mode == "asinh_ratio_minus1":
        return np.arcsinh((Lt / Ls) - 1.0)

    if y_mode == "signedlog_ratio_minus1":
        r = (Lt / Ls) - 1.0
        return np.sign(r) * np.log1p(np.abs(r))

    if y_mode == "k_ratio_minus1":
        return (Ls / Lt) - 1.0

    if y_mode == "asinh_k_ratio_minus1":
        return np.arcsinh((Ls / Lt) - 1.0)

    if y_mode == "signedlog_k_ratio_minus1":
        r = (Ls / Lt) - 1.0
        return np.sign(r) * np.log1p(np.abs(r))

    raise ValueError(f"Unknown y_mode: {y_mode}")


def reconstruct_L(L_lin: np.ndarray, y_pred: np.ndarray, strat: Strategy) -> np.ndarray:
    if strat.y_mode == "log_ratio":
        y = np.clip(y_pred, -Y_CLIP, Y_CLIP)
        return L_lin * np.exp(y)

    if strat.y_mode == "ratio_minus1":
        corr = 1.0 + y_pred
        corr = np.maximum(corr, 1e-12)
        return L_lin * corr

    if strat.y_mode == "asinh_ratio_minus1":
        corr = 1.0 + np.sinh(y_pred)
        corr = np.maximum(corr, 1e-12)
        return L_lin * corr

    if strat.y_mode == "signedlog_ratio_minus1":
        r = np.sign(y_pred) * (np.exp(np.abs(y_pred)) - 1.0)
        corr = 1.0 + r
        corr = np.maximum(corr, 1e-12)
        return L_lin * corr

    if strat.y_mode == "k_ratio_minus1":
        corr_k = 1.0 + y_pred
        corr_k = np.maximum(corr_k, 1e-12)
        return L_lin / corr_k

    if strat.y_mode == "asinh_k_ratio_minus1":
        corr_k = 1.0 + np.sinh(y_pred)
        corr_k = np.maximum(corr_k, 1e-12)
        return L_lin / corr_k

    if strat.y_mode == "signedlog_k_ratio_minus1":
        r = np.sign(y_pred) * (np.exp(np.abs(y_pred)) - 1.0)
        corr_k = 1.0 + r
        corr_k = np.maximum(corr_k, 1e-12)
        return L_lin / corr_k

    raise ValueError(f"Unknown y_mode: {strat.y_mode}")


# ==============================================================================
# FEATURE SELECTION (Spearman + mRMR) - unchanged
# ==============================================================================
def _rankdata(a: np.ndarray) -> np.ndarray:
    temp = np.argsort(a, kind="mergesort")
    ranks = np.empty_like(temp, dtype=np.float64)
    ranks[temp] = np.arange(len(a), dtype=np.float64)

    sa = a[temp]
    i = 0
    n = len(a)
    while i < n:
        j = i
        while (j + 1) < n and sa[j + 1] == sa[i]:
            j += 1
        if j > i:
            avg = 0.5 * (i + j)
            ranks[temp[i:j + 1]] = avg
        i = j + 1
    return ranks


def spearman_corr(x: np.ndarray, y: np.ndarray) -> float:
    rx = _rankdata(x)
    ry = _rankdata(y)
    rxm = rx - rx.mean()
    rym = ry - ry.mean()
    denom = math.sqrt(float((rxm * rxm).sum() * (rym * rym).sum()))
    if denom <= 0.0:
        return 0.0
    return float((rxm * rym).sum() / denom)


def mrmr_select_3(X: np.ndarray, y: np.ndarray, alpha: float = 0.5) -> Tuple[List[int], np.ndarray]:
    p = X.shape[1]
    corr_y = np.array([abs(spearman_corr(X[:, j], y)) for j in range(p)], dtype=np.float64)

    sel = [int(np.argmax(corr_y))]
    while len(sel) < 3:
        best_score = -1e99
        best_j = -1
        for j in range(p):
            if j in sel:
                continue
            redundancy = 0.0
            for s in sel:
                redundancy = max(redundancy, abs(spearman_corr(X[:, j], X[:, s])))
            score = float(corr_y[j] - alpha * redundancy)
            if score > best_score:
                best_score = score
                best_j = j
        sel.append(best_j)

    return sel, corr_y


# ==============================================================================
# ELM-style MLP: random hidden layer (tanh) + trained linear output
# ==============================================================================
@dataclass
class ELMParams:
    x_mean: np.ndarray
    x_std: np.ndarray
    y_mean: float
    y_std: float
    W: np.ndarray     # (3, H)
    b: np.ndarray     # (H,)
    beta: np.ndarray  # (H+1,) includes bias term at the end
    hidden: int


def _hidden_act(Xn: np.ndarray, W: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.tanh(Xn @ W + b)


def _solve_weighted_ridge(Phi: np.ndarray, y: np.ndarray, lam: float, w: np.ndarray | None) -> np.ndarray:
    # Phi: (n, m), y: (n,)
    if w is None:
        A = Phi
        b = y
    else:
        sw = np.sqrt(np.maximum(w, 1e-12)).reshape(-1, 1)
        A = Phi * sw
        b = y * sw[:, 0]

    ATA = A.T @ A
    ATA[np.diag_indices_from(ATA)] += lam
    ATb = A.T @ b
    return np.linalg.solve(ATA, ATb)


def predict_y_elm(p: ELMParams, X3: np.ndarray) -> np.ndarray:
    Xn = (X3 - p.x_mean) / p.x_std
    Hh = _hidden_act(Xn, p.W, p.b)
    Phi = np.concatenate([Hh, np.ones((Hh.shape[0], 1), dtype=np.float64)], axis=1)
    y_scaled = Phi @ p.beta
    return p.y_mean + p.y_std * y_scaled


def _compute_reweight(L_lin: np.ndarray, L_true: np.ndarray, y_pred: np.ndarray, strat: Strategy) -> np.ndarray:
    L_pred = reconstruct_L(L_lin, y_pred, strat)
    ape = np.abs(L_true - L_pred) / np.maximum(L_true, 1e-12)

    if strat.reweight_mode == "none":
        w = np.ones_like(ape)
    elif strat.reweight_mode == "mape":
        w = 1.0 + strat.reweight_beta * ape
    elif strat.reweight_mode == "tail":
        w = 1.0 + strat.reweight_beta * np.power(ape, strat.reweight_gamma)
    elif strat.reweight_mode == "p95":
        p95 = float(np.quantile(ape, 0.95))
        p95 = max(p95, 1e-6)
        w = 1.0 + strat.reweight_beta * np.power(ape / p95, strat.reweight_gamma)
    elif strat.reweight_mode == "topq":
        tq = float(np.quantile(ape, 0.90))
        tq = max(tq, 1e-6)
        boost = (ape >= tq).astype(np.float64)
        w = 1.0 + strat.reweight_beta * boost * np.power(ape / tq, strat.reweight_gamma)
    elif strat.reweight_mode == "huber_down":
        c = float(np.quantile(ape, 0.90))
        c = max(c, 1e-6)
        w = 1.0 / (1.0 + np.power(ape / c, 2.0))
    else:
        raise ValueError(f"Unknown reweight_mode: {strat.reweight_mode}")

    return np.clip(w, WEIGHT_CLIP_MIN, WEIGHT_CLIP_MAX)


def fit_elm(
    X3: np.ndarray,
    y: np.ndarray,
    L_lin: np.ndarray,
    L_true: np.ndarray,
    strat: Strategy,
    degree: int,
) -> ELMParams:
    # Standardize input and target for conditioning
    x_mean = X3.mean(axis=0)
    x_std = X3.std(axis=0)
    x_std[x_std == 0.0] = 1.0

    y_mean = float(np.mean(y))
    y_std = float(np.std(y))
    if y_std <= 0.0:
        y_std = 1.0

    Xn = (X3 - x_mean) / x_std
    ys = (y - y_mean) / y_std

    hidden = int(HIDDEN_MULT * degree)
    rng = np.random.default_rng(RNG_SEED + 1000 * degree)

    # Random hidden layer parameters (fixed after init)
    W = rng.normal(scale=math.sqrt(1.0 / 3.0), size=(3, hidden))
    b = rng.normal(scale=1.0, size=(hidden,))

    Hh = _hidden_act(Xn, W, b)
    Phi = np.concatenate([Hh, np.ones((Hh.shape[0], 1), dtype=np.float64)], axis=1)

    beta = _solve_weighted_ridge(Phi, ys, lam=strat.lam, w=None)

    # Optional L-space reweight iterations: keep (W,b) fixed, refit beta
    if strat.reweight_mode != "none" and strat.reweight_iters > 0:
        for _ in range(int(strat.reweight_iters)):
            y_pred = y_mean + y_std * (Phi @ beta)
            w = _compute_reweight(L_lin, L_true, y_pred, strat)
            beta = _solve_weighted_ridge(Phi, ys, lam=strat.lam, w=w)

    return ELMParams(
        x_mean=x_mean,
        x_std=x_std,
        y_mean=y_mean,
        y_std=y_std,
        W=W,
        b=b,
        beta=beta,
        hidden=hidden,
    )


# ==============================================================================
# METRICS + PLOTS
# ==============================================================================
def metrics(L_true: np.ndarray, L_pred: np.ndarray) -> Dict[str, float]:
    err = L_true - L_pred
    ape = np.abs(err) / np.maximum(L_true, 1e-12)
    return {
        "mse": float(np.mean(err * err)),
        "mape_pct": float(np.mean(ape) * 100.0),
        "p95_ape_pct": float(np.quantile(ape, 0.95) * 100.0),
        "max_ape_pct": float(np.max(ape) * 100.0),
    }


def save_plots(df: pd.DataFrame, pdf_path: str = "diagnostics.pdf", png_path: str = "plot_all.png") -> None:
    with PdfPages(pdf_path) as pdf:
        plt.figure(figsize=(10, 8))
        plt.scatter(df["L_true"], df["L_pred"], s=6, alpha=0.25)
        mx = max(df["L_true"].max(), df["L_pred"].max())
        plt.plot([0, mx], [0, mx], "k--", lw=2)
        plt.xlabel("True L (m)")
        plt.ylabel("Predicted L (m)")
        plt.title("Parity plot")
        plt.grid(True, alpha=0.3)
        pdf.savefig()
        plt.savefig(png_path, dpi=200, bbox_inches="tight")
        plt.close()

        plt.figure(figsize=(10, 6))
        plt.hist(df["APE_pct"], bins=60, alpha=0.8, edgecolor="black")
        plt.xlabel("Absolute Percentage Error (%)")
        plt.ylabel("Count")
        plt.title("Error distribution")
        plt.grid(True, alpha=0.3)
        pdf.savefig()
        plt.close()

        for xcol, title in [("Uc", "Error vs current Uc"), ("d", "Error vs depth d"), ("T", "Error vs period T")]:
            plt.figure(figsize=(10, 7))
            plt.scatter(df[xcol], df["APE_pct"], s=6, alpha=0.25)
            plt.xlabel(xcol)
            plt.ylabel("Absolute Percentage Error (%)")
            plt.title(title)
            plt.grid(True, alpha=0.3)
            pdf.savefig()
            plt.close()

    print(f"[SUCCESS] Saved {pdf_path} and {png_path}")


# ==============================================================================
# function.py GENERATION (self-contained ELM)
# ==============================================================================
def generate_function_py(out_path: str, strat: Strategy, feature_idx: List[int], p: ELMParams) -> None:
    # Scalar expressions mirroring build_feature_bank
    bank_expr = {
        "log_groups": [
            "math.log(max(H/L_lin, 1e-12))",
            "math.log(max(d/L_lin, 1e-12))",
            "(Uc*T)/max(L_lin, 1e-12)",
            "math.log(max((H*(L_lin**2))/max(d,1e-12)**3, 1e-12))",
            "Uc/math.sqrt(G*max(d,1e-12))",
        ],
        "raw_groups": [
            "H/max(L_lin, 1e-12)",
            "d/max(L_lin, 1e-12)",
            "(Uc*T)/max(L_lin, 1e-12)",
            "(H*(L_lin**2))/max(d,1e-12)**3",
            "Uc/math.sqrt(G*max(d,1e-12))",
        ],
        "asinh_groups": [
            "math.asinh(H/max(L_lin, 1e-12))",
            "math.asinh(d/max(L_lin, 1e-12))",
            "math.asinh((Uc*T)/max(L_lin, 1e-12))",
            "math.asinh((H*(L_lin**2))/max(d,1e-12)**3)",
            "math.asinh(Uc/math.sqrt(G*max(d,1e-12)))",
        ],
        "log_asinh_groups": [
            "math.log(max(H/L_lin, 1e-12))",
            "math.log(max(d/L_lin, 1e-12))",
            "math.asinh((Uc*T)/max(L_lin, 1e-12))",
            "math.log(max((H*(L_lin**2))/max(d,1e-12)**3, 1e-12))",
            "math.asinh(Uc/math.sqrt(G*max(d,1e-12)))",
        ],
        "signedlog_groups": [
            "math.log(max(H/L_lin, 1e-12))",
            "math.log(max(d/L_lin, 1e-12))",
            "_signed_log1p((Uc*T)/max(L_lin, 1e-12))",
            "math.log(max((H*(L_lin**2))/max(d,1e-12)**3, 1e-12))",
            "_signed_log1p(Uc/math.sqrt(G*max(d,1e-12)))",
        ],
        "raw_signedlog_groups": [
            "H/max(L_lin, 1e-12)",
            "d/max(L_lin, 1e-12)",
            "_signed_log1p((Uc*T)/max(L_lin, 1e-12))",
            "(H*(L_lin**2))/max(d,1e-12)**3",
            "_signed_log1p(Uc/math.sqrt(G*max(d,1e-12)))",
        ],
    }[strat.x_mode]

    x0_expr = bank_expr[feature_idx[0]]
    x1_expr = bank_expr[feature_idx[1]]
    x2_expr = bank_expr[feature_idx[2]]

    m0, m1, m2 = p.x_mean.astype(np.float64).tolist()
    s0, s1, s2 = p.x_std.astype(np.float64).tolist()

    # For scalar evaluation store per-neuron weights
    W = p.W.astype(np.float64).T.tolist()   # (H, 3)
    b = p.b.astype(np.float64).tolist()
    beta = p.beta.astype(np.float64).tolist()

    code = f'''# -*- coding: utf-8 -*-
"""
function.py
===========
Generated by pade_nn.py (single global neural network model: tanh hidden layer + linear output)

Inputs:
  H  (m)   wave height
  T  (s)   wave period (observed / lab frame)
  d  (m)   mean water depth
  Uc (m/s) Eulerian current (+ with waves)

Output:
  L  (m)   wavelength
"""
from __future__ import annotations

import math
import sys

G = {G:.16e}
PI = {PI:.16e}

Y_MODE = "{strat.y_mode}"
X_MODE = "{strat.x_mode}"
FEATURE_INDEX = {feature_idx}

M0, M1, M2 = {m0:.16e}, {m1:.16e}, {m2:.16e}
S0, S1, S2 = {s0:.16e}, {s1:.16e}, {s2:.16e}

Y_MEAN = {p.y_mean:.16e}
Y_STD  = {p.y_std:.16e}

Y_CLIP = {Y_CLIP:.16e}

# Hidden layer parameters (H neurons)
W = {W}          # list of neurons; each neuron has 3 weights
b = {b}          # list length H
beta = {beta}    # list length H+1 (last element is linear output bias term)

def _signed_log1p(x: float) -> float:
    return math.copysign(math.log1p(abs(x)), x)

def _signed_expm1(x: float) -> float:
    return math.copysign(math.expm1(abs(x)), x)

def _solve_k_no_current(omega: float, d: float) -> float:
    k_deep = (omega * omega) / G
    c_shallow = math.sqrt(G * max(d, 1e-12))
    k_shallow = omega / max(c_shallow, 1e-12)
    k = k_deep if (k_deep * d) > 1.0 else k_shallow
    for _ in range(8):
        k = max(1e-9, min(k, 400.0))
        kd = k * d
        th = math.tanh(kd)
        sigma = math.sqrt(max(G * k * th, 1e-16))
        f = sigma - omega
        sech2 = 1.0 - th * th
        if sigma > 1e-12:
            d_sigma = (G * th + G * kd * sech2) / (2.0 * sigma)
        else:
            d_sigma = 0.0
        if abs(d_sigma) < 1e-12:
            break
        k_new = k - f / d_sigma
        if k_new <= 1e-9:
            k_new = 1e-9
        k = 0.7 * k + 0.3 * k_new
    return max(k, 1e-9)

def _solve_k_with_current(omega: float, d: float, Uc: float) -> float:
    k = _solve_k_no_current(omega, d)
    for _ in range(40):
        k = max(1e-9, min(k, 400.0))
        kd = k * d
        th = math.tanh(kd)
        sigma = math.sqrt(max(G * k * th, 1e-16))
        f = sigma + k * Uc - omega
        sech2 = 1.0 - th * th
        if sigma > 1e-12:
            d_sigma = (G * th + G * kd * sech2) / (2.0 * sigma)
        else:
            d_sigma = 0.0
        df = d_sigma + Uc
        if abs(df) < 1e-12:
            break
        k_new = k - f / df
        if k_new <= 1e-9:
            k_new = 1e-9
        k = 0.8 * k + 0.2 * k_new
        if abs(k_new - k) / max(k, 1e-9) < 1e-10:
            break
    return max(k, 1e-9)

def _elm_y_scaled(x0: float, x1: float, x2: float) -> float:
    # Hidden activations
    acc = beta[-1]  # output bias term
    for wi, bi, betai in zip(W, b, beta[:-1]):
        z = bi + wi[0]*x0 + wi[1]*x1 + wi[2]*x2
        acc += math.tanh(z) * betai
    return acc

def L(H: float, T: float, d: float, Uc: float = 0.0) -> float:
    omega = 2.0 * PI / T
    k = _solve_k_with_current(omega, d, Uc)
    L_lin = 2.0 * PI / k

    x0_raw = {x0_expr}
    x1_raw = {x1_expr}
    x2_raw = {x2_expr}

    x0 = (x0_raw - M0) / S0
    x1 = (x1_raw - M1) / S1
    x2 = (x2_raw - M2) / S2

    y_scaled = _elm_y_scaled(x0, x1, x2)
    y = Y_MEAN + Y_STD * y_scaled

    if Y_MODE == "log_ratio":
        y = max(-Y_CLIP, min(Y_CLIP, y))
        return L_lin * math.exp(y)

    if Y_MODE == "ratio_minus1":
        corr = 1.0 + y
        if corr <= 1e-12:
            corr = 1e-12
        return L_lin * corr

    if Y_MODE == "asinh_ratio_minus1":
        corr = 1.0 + math.sinh(y)
        if corr <= 1e-12:
            corr = 1e-12
        return L_lin * corr

    if Y_MODE == "signedlog_ratio_minus1":
        r = _signed_expm1(y)
        corr = 1.0 + r
        if corr <= 1e-12:
            corr = 1e-12
        return L_lin * corr

    if Y_MODE == "k_ratio_minus1":
        corr_k = 1.0 + y
        if corr_k <= 1e-12:
            corr_k = 1e-12
        return L_lin / corr_k

    if Y_MODE == "asinh_k_ratio_minus1":
        corr_k = 1.0 + math.sinh(y)
        if corr_k <= 1e-12:
            corr_k = 1e-12
        return L_lin / corr_k

    if Y_MODE == "signedlog_k_ratio_minus1":
        r = _signed_expm1(y)
        corr_k = 1.0 + r
        if corr_k <= 1e-12:
            corr_k = 1e-12
        return L_lin / corr_k

    raise ValueError("Unknown Y_MODE: " + str(Y_MODE))

def _usage() -> None:
    print("Usage:")
    print("  python function.py H T d Uc")
    print("Example:")
    print("  python function.py 3 9 5 1")

if __name__ == "__main__":
    if len(sys.argv) == 5:
        H = float(sys.argv[1]); T = float(sys.argv[2]); d = float(sys.argv[3]); Uc = float(sys.argv[4])
        print("L = %.10f m" % (L(H, T, d, Uc),))
    elif len(sys.argv) == 1:
        H = float(input("Wave height H (m) [3.0]: ") or "3.0")
        T = float(input("Wave period T (s) [9.0]: ") or "9.0")
        d = float(input("Water depth d (m) [5.0]: ") or "5.0")
        Uc = float(input("Eulerian current Uc (m/s, + with wave) [1.0]: ") or "1.0")
        print("L = %.10f m" % (L(H, T, d, Uc),))
    else:
        _usage()
'''
    Path(out_path).write_text(code, encoding="utf-8")
    print(f"[SUCCESS] Wrote {out_path}")


# ==============================================================================
# neural.bas GENERATION (self-contained VBA module)
# ==============================================================================
def _vba_num(x: float) -> str:
    xf = float(x)
    if not np.isfinite(xf):
        raise ValueError(f"Cannot export non-finite VBA coefficient: {x!r}")
    return f"{xf:.16E}#"


def _vba_assign_array(lines: List[str], name: str, values: np.ndarray) -> None:
    vals = np.asarray(values, dtype=np.float64).ravel()
    for i, v in enumerate(vals):
        lines.append(f"    {name}({i}) = {_vba_num(v)}")


def _vba_feature_lines(x_mode: str, feature_idx: List[int]) -> List[str]:
    lines: List[str] = [
        "    Dim L_feat As Double",
        "    Dim d_safe As Double",
        "    Dim f(0 To 4) As Double",
        "",
        "    L_feat = MaxVal(L_base, 0.1#)",
        "    d_safe = MaxVal(d, 0.1#)",
    ]

    if x_mode == "log_groups":
        body = [
            "    f(0) = Log(MaxVal(H / L_feat, EPS_LOG))",
            "    f(1) = Log(MaxVal(d / L_feat, EPS_LOG))",
            "    f(2) = (Uc * T) / L_feat",
            "    f(3) = Log(MaxVal((H * L_feat * L_feat) / _",
            "                    (d_safe * d_safe * d_safe), EPS_LOG))",
            "    f(4) = Uc / Sqr(G * d_safe)",
        ]
    elif x_mode == "raw_groups":
        body = [
            "    f(0) = H / L_feat",
            "    f(1) = d / L_feat",
            "    f(2) = (Uc * T) / L_feat",
            "    f(3) = (H * L_feat * L_feat) / (d_safe * d_safe * d_safe)",
            "    f(4) = Uc / Sqr(G * d_safe)",
        ]
    elif x_mode == "asinh_groups":
        body = [
            "    f(0) = FnAsinh(H / L_feat)",
            "    f(1) = FnAsinh(d / L_feat)",
            "    f(2) = FnAsinh((Uc * T) / L_feat)",
            "    f(3) = FnAsinh((H * L_feat * L_feat) / _",
            "                    (d_safe * d_safe * d_safe))",
            "    f(4) = FnAsinh(Uc / Sqr(G * d_safe))",
        ]
    elif x_mode == "log_asinh_groups":
        body = [
            "    f(0) = Log(MaxVal(H / L_feat, EPS_LOG))",
            "    f(1) = Log(MaxVal(d / L_feat, EPS_LOG))",
            "    f(2) = FnAsinh((Uc * T) / L_feat)",
            "    f(3) = Log(MaxVal((H * L_feat * L_feat) / _",
            "                    (d_safe * d_safe * d_safe), EPS_LOG))",
            "    f(4) = FnAsinh(Uc / Sqr(G * d_safe))",
        ]
    elif x_mode == "signedlog_groups":
        body = [
            "    f(0) = Log(MaxVal(H / L_feat, EPS_LOG))",
            "    f(1) = Log(MaxVal(d / L_feat, EPS_LOG))",
            "    f(2) = SignedLog1p((Uc * T) / L_feat)",
            "    f(3) = Log(MaxVal((H * L_feat * L_feat) / _",
            "                    (d_safe * d_safe * d_safe), EPS_LOG))",
            "    f(4) = SignedLog1p(Uc / Sqr(G * d_safe))",
        ]
    elif x_mode == "raw_signedlog_groups":
        body = [
            "    f(0) = H / L_feat",
            "    f(1) = d / L_feat",
            "    f(2) = SignedLog1p((Uc * T) / L_feat)",
            "    f(3) = (H * L_feat * L_feat) / (d_safe * d_safe * d_safe)",
            "    f(4) = SignedLog1p(Uc / Sqr(G * d_safe))",
        ]
    else:
        raise ValueError(f"Unsupported x_mode for VBA export: {x_mode}")

    lines.extend(body)
    lines.extend([
        "",
        f"    x0_raw = f({int(feature_idx[0])})",
        f"    x1_raw = f({int(feature_idx[1])})",
        f"    x2_raw = f({int(feature_idx[2])})",
    ])
    return lines


def generate_neural_bas(out_path: str | Path, strat: Strategy, feature_idx: List[int], p: ELMParams) -> None:
    hidden = int(p.hidden)
    w = np.asarray(p.W, dtype=np.float64)
    b = np.asarray(p.b, dtype=np.float64)
    beta = np.asarray(p.beta, dtype=np.float64)

    if w.shape != (3, hidden):
        raise ValueError(f"Unexpected W shape for VBA export: {w.shape}")
    if b.shape != (hidden,):
        raise ValueError(f"Unexpected b shape for VBA export: {b.shape}")
    if beta.shape != (hidden + 1,):
        raise ValueError(f"Unexpected beta shape for VBA export: {beta.shape}")

    lines: List[str] = []
    lines.append('Attribute VB_Name = "neural"')
    lines.append('Option Explicit')
    lines.append('')
    lines.append(f'Private Const G As Double = {_vba_num(G)}')
    lines.append(f'Private Const PI As Double = {_vba_num(PI)}')
    lines.append('Private Const EPS_K As Double = 0.000000001#')
    lines.append('Private Const EPS_LOG As Double = 0.000000000001#')
    lines.append(f'Private Const Y_CLIP As Double = {_vba_num(Y_CLIP)}')
    lines.append(f'Private Const HIDDEN_NEURONS As Long = {hidden}')
    lines.append(f'Private Const Y_MODE As String = "{strat.y_mode}"')
    lines.append('')
    lines.append(f'Private Const X0_MEAN As Double = {_vba_num(p.x_mean[0])}')
    lines.append(f'Private Const X1_MEAN As Double = {_vba_num(p.x_mean[1])}')
    lines.append(f'Private Const X2_MEAN As Double = {_vba_num(p.x_mean[2])}')
    lines.append(f'Private Const X0_STD As Double = {_vba_num(p.x_std[0])}')
    lines.append(f'Private Const X1_STD As Double = {_vba_num(p.x_std[1])}')
    lines.append(f'Private Const X2_STD As Double = {_vba_num(p.x_std[2])}')
    lines.append(f'Private Const Y_MEAN As Double = {_vba_num(p.y_mean)}')
    lines.append(f'Private Const Y_STD As Double = {_vba_num(p.y_std)}')
    lines.append('')
    lines.append('Private ModelReady As Boolean')
    lines.append('Private W0() As Double')
    lines.append('Private W1() As Double')
    lines.append('Private W2() As Double')
    lines.append('Private BHid() As Double')
    lines.append('Private BetaH() As Double')
    lines.append('Private BetaBias As Double')
    lines.append('')
    lines.append('Public Function L_neural(ByVal H As Double, ByVal T As Double, _')
    lines.append('                         ByVal d As Double, ByVal Uc As Double) As Double')
    lines.append('    Dim L_lin As Double')
    lines.append('')
    lines.append('    If T <= 0# Or d <= 0# Then')
    lines.append('        L_neural = 0#')
    lines.append('        Exit Function')
    lines.append('    End If')
    lines.append('')
    lines.append('    L_lin = SolveLinearDoppler(T, d, Uc)')
    lines.append('    L_neural = PredictNN(H, T, d, Uc, L_lin)')
    lines.append('End Function')
    lines.append('')
    lines.append('Public Function L(ByVal H As Double, ByVal T As Double, _')
    lines.append('                  ByVal d As Double, ByVal Uc As Double) As Double')
    lines.append('    L = L_neural(H, T, d, Uc)')
    lines.append('End Function')
    lines.append('')
    lines.append('Private Sub InitModel()')
    lines.append('    If ModelReady Then Exit Sub')
    lines.append('')
    lines.append('    ReDim W0(0 To HIDDEN_NEURONS - 1)')
    lines.append('    ReDim W1(0 To HIDDEN_NEURONS - 1)')
    lines.append('    ReDim W2(0 To HIDDEN_NEURONS - 1)')
    lines.append('    ReDim BHid(0 To HIDDEN_NEURONS - 1)')
    lines.append('    ReDim BetaH(0 To HIDDEN_NEURONS - 1)')
    lines.append('')
    _vba_assign_array(lines, 'W0', w[0, :])
    _vba_assign_array(lines, 'W1', w[1, :])
    _vba_assign_array(lines, 'W2', w[2, :])
    _vba_assign_array(lines, 'BHid', b)
    _vba_assign_array(lines, 'BetaH', beta[:-1])
    lines.append(f'    BetaBias = {_vba_num(beta[-1])}')
    lines.append('')
    lines.append('    ModelReady = True')
    lines.append('End Sub')
    lines.append('')
    lines.append('Private Function PredictNN(ByVal H As Double, ByVal T As Double, _')
    lines.append('                           ByVal d As Double, ByVal Uc As Double, _')
    lines.append('                           ByVal L_base As Double) As Double')
    lines.append('    Dim x0_raw As Double, x1_raw As Double, x2_raw As Double')
    lines.append('    Dim x0 As Double, x1 As Double, x2 As Double')
    lines.append('    Dim acc As Double, y As Double, corr As Double')
    lines.append('    Dim i As Long')
    lines.append('')
    lines.append('    InitModel')
    lines.append('    Call BuildSelectedFeatures(H, T, d, Uc, L_base, x0_raw, x1_raw, x2_raw)')
    lines.append('')
    lines.append('    x0 = (x0_raw - X0_MEAN) / X0_STD')
    lines.append('    x1 = (x1_raw - X1_MEAN) / X1_STD')
    lines.append('    x2 = (x2_raw - X2_MEAN) / X2_STD')
    lines.append('')
    lines.append('    acc = BetaBias')
    lines.append('    For i = 0 To HIDDEN_NEURONS - 1')
    lines.append('        acc = acc + BetaH(i) * FnTanh(BHid(i) + W0(i) * x0 + _')
    lines.append('                                      W1(i) * x1 + W2(i) * x2)')
    lines.append('    Next i')
    lines.append('')
    lines.append('    y = Y_MEAN + Y_STD * acc')
    lines.append('')
    lines.append('    Select Case Y_MODE')
    lines.append('    Case "log_ratio"')
    lines.append('        If y < -Y_CLIP Then y = -Y_CLIP')
    lines.append('        If y > Y_CLIP Then y = Y_CLIP')
    lines.append('        PredictNN = MaxVal(L_base, 0.1#) * Exp(y)')
    lines.append('')
    lines.append('    Case "ratio_minus1"')
    lines.append('        corr = 1# + y')
    lines.append('        If corr <= 0.000000000001# Then corr = 0.000000000001#')
    lines.append('        PredictNN = MaxVal(L_base, 0.1#) * corr')
    lines.append('')
    lines.append('    Case "asinh_ratio_minus1"')
    lines.append('        corr = 1# + FnSinh(y)')
    lines.append('        If corr <= 0.000000000001# Then corr = 0.000000000001#')
    lines.append('        PredictNN = MaxVal(L_base, 0.1#) * corr')
    lines.append('')
    lines.append('    Case "signedlog_ratio_minus1"')
    lines.append('        corr = 1# + SignedExpm1(y)')
    lines.append('        If corr <= 0.000000000001# Then corr = 0.000000000001#')
    lines.append('        PredictNN = MaxVal(L_base, 0.1#) * corr')
    lines.append('')
    lines.append('    Case "k_ratio_minus1"')
    lines.append('        corr = 1# + y')
    lines.append('        If corr <= 0.000000000001# Then corr = 0.000000000001#')
    lines.append('        PredictNN = MaxVal(L_base, 0.1#) / corr')
    lines.append('')
    lines.append('    Case "asinh_k_ratio_minus1"')
    lines.append('        corr = 1# + FnSinh(y)')
    lines.append('        If corr <= 0.000000000001# Then corr = 0.000000000001#')
    lines.append('        PredictNN = MaxVal(L_base, 0.1#) / corr')
    lines.append('')
    lines.append('    Case "signedlog_k_ratio_minus1"')
    lines.append('        corr = 1# + SignedExpm1(y)')
    lines.append('        If corr <= 0.000000000001# Then corr = 0.000000000001#')
    lines.append('        PredictNN = MaxVal(L_base, 0.1#) / corr')
    lines.append('')
    lines.append('    Case Else')
    lines.append('        PredictNN = MaxVal(L_base, 0.1#)')
    lines.append('    End Select')
    lines.append('End Function')
    lines.append('')
    lines.append('Private Sub BuildSelectedFeatures(ByVal H As Double, ByVal T As Double, _')
    lines.append('                                  ByVal d As Double, ByVal Uc As Double, _')
    lines.append('                                  ByVal L_base As Double, _')
    lines.append('                                  ByRef x0_raw As Double, _')
    lines.append('                                  ByRef x1_raw As Double, _')
    lines.append('                                  ByRef x2_raw As Double)')
    lines.extend(_vba_feature_lines(strat.x_mode, feature_idx))
    lines.append('End Sub')
    lines.append('')
    lines.append('Private Function SolveLinearDoppler(ByVal T As Double, ByVal d As Double, _')
    lines.append('                                    ByVal Uc As Double) As Double')
    lines.append('    Dim omega As Double')
    lines.append('    Dim k As Double')
    lines.append('    Dim k_deep As Double')
    lines.append('    Dim c_shallow As Double')
    lines.append('    Dim k_shallow As Double')
    lines.append('    Dim kd As Double')
    lines.append('    Dim th As Double')
    lines.append('    Dim sigma As Double')
    lines.append('    Dim f As Double')
    lines.append('    Dim sech2 As Double')
    lines.append('    Dim d_sigma As Double')
    lines.append('    Dim df As Double')
    lines.append('    Dim k_new As Double')
    lines.append('    Dim it As Long')
    lines.append('')
    lines.append('    omega = 2# * PI / T')
    lines.append('    k_deep = (omega * omega) / G')
    lines.append('    c_shallow = Sqr(G * MaxVal(d, 0.000000000001#))')
    lines.append('    k_shallow = omega / MaxVal(c_shallow, 0.000000000001#)')
    lines.append('')
    lines.append('    If k_deep * d > 1# Then')
    lines.append('        k = k_deep')
    lines.append('    Else')
    lines.append('        k = k_shallow')
    lines.append('    End If')
    lines.append('')
    lines.append('    For it = 1 To 8')
    lines.append('        If k < EPS_K Then k = EPS_K')
    lines.append('        If k > 400# Then k = 400#')
    lines.append('        kd = k * d')
    lines.append('        th = FnTanh(kd)')
    lines.append('        sigma = Sqr(MaxVal(G * k * th, 0.0000000000000001#))')
    lines.append('        f = sigma - omega')
    lines.append('        sech2 = 1# - th * th')
    lines.append('        If sigma > 0.000000000001# Then')
    lines.append('            d_sigma = (G * th + G * kd * sech2) / (2# * sigma)')
    lines.append('        Else')
    lines.append('            d_sigma = 0#')
    lines.append('        End If')
    lines.append('        If Abs(d_sigma) < 0.000000000001# Then Exit For')
    lines.append('        k_new = k - f / d_sigma')
    lines.append('        If k_new <= EPS_K Then k_new = EPS_K')
    lines.append('        k = 0.7# * k + 0.3# * k_new')
    lines.append('    Next it')
    lines.append('')
    lines.append('    For it = 1 To 40')
    lines.append('        If k < EPS_K Then k = EPS_K')
    lines.append('        If k > 400# Then k = 400#')
    lines.append('        kd = k * d')
    lines.append('        th = FnTanh(kd)')
    lines.append('        sigma = Sqr(MaxVal(G * k * th, 0.0000000000000001#))')
    lines.append('        f = sigma + k * Uc - omega')
    lines.append('        sech2 = 1# - th * th')
    lines.append('        If sigma > 0.000000000001# Then')
    lines.append('            d_sigma = (G * th + G * kd * sech2) / (2# * sigma)')
    lines.append('        Else')
    lines.append('            d_sigma = 0#')
    lines.append('        End If')
    lines.append('        df = d_sigma + Uc')
    lines.append('        If Abs(df) < 0.000000000001# Then Exit For')
    lines.append('        k_new = k - f / df')
    lines.append('        If k_new <= EPS_K Then k_new = EPS_K')
    lines.append('        k = 0.8# * k + 0.2# * k_new')
    lines.append('        If Abs(k_new - k) / MaxVal(k, EPS_K) < 0.0000000001# Then Exit For')
    lines.append('    Next it')
    lines.append('')
    lines.append('    If k <= EPS_K Then k = EPS_K')
    lines.append('    SolveLinearDoppler = 2# * PI / k')
    lines.append('End Function')
    lines.append('')
    lines.append('Private Function MaxVal(ByVal a As Double, ByVal b As Double) As Double')
    lines.append('    If a >= b Then')
    lines.append('        MaxVal = a')
    lines.append('    Else')
    lines.append('        MaxVal = b')
    lines.append('    End If')
    lines.append('End Function')
    lines.append('')
    lines.append('Private Function FnTanh(ByVal x As Double) As Double')
    lines.append('    If x > 20# Then')
    lines.append('        FnTanh = 1#')
    lines.append('    ElseIf x < -20# Then')
    lines.append('        FnTanh = -1#')
    lines.append('    Else')
    lines.append('        FnTanh = (Exp(x) - Exp(-x)) / (Exp(x) + Exp(-x))')
    lines.append('    End If')
    lines.append('End Function')
    lines.append('')
    lines.append('Private Function FnSinh(ByVal x As Double) As Double')
    lines.append('    FnSinh = 0.5# * (Exp(x) - Exp(-x))')
    lines.append('End Function')
    lines.append('')
    lines.append('Private Function FnAsinh(ByVal x As Double) As Double')
    lines.append('    FnAsinh = Log(x + Sqr(x * x + 1#))')
    lines.append('End Function')
    lines.append('')
    lines.append('Private Function SignedLog1p(ByVal x As Double) As Double')
    lines.append('    If x >= 0# Then')
    lines.append('        SignedLog1p = Log(1# + x)')
    lines.append('    Else')
    lines.append('        SignedLog1p = -Log(1# + Abs(x))')
    lines.append('    End If')
    lines.append('End Function')
    lines.append('')
    lines.append('Private Function SignedExpm1(ByVal x As Double) As Double')
    lines.append('    If x >= 0# Then')
    lines.append('        SignedExpm1 = Exp(x) - 1#')
    lines.append('    Else')
    lines.append('        SignedExpm1 = -(Exp(-x) - 1#)')
    lines.append('    End If')
    lines.append('End Function')

    Path(out_path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[SUCCESS] Wrote {out_path}")


# ==============================================================================
# MAIN
# ==============================================================================
def main() -> None:
    print("--- Loading Data ---")
    df = load_data("list.txt")

    max_train_samples = infer_max_train_samples(df)
    print(f"[INFO] Inferred MAX_TRAIN_SAMPLES from list.txt: {max_train_samples}")

    H = df["H"].to_numpy(np.float64)
    T = df["T"].to_numpy(np.float64)
    d = df["d"].to_numpy(np.float64)
    Uc = df["Uc"].to_numpy(np.float64)
    L_true = df["L"].to_numpy(np.float64)

    print("--- Physics Baseline (linear dispersion + current) ---")
    L_lin = solve_linear_wavelength(T, d, Uc)

    rng = np.random.default_rng(123)
    n = len(L_true)
    idx_all = rng.permutation(n)
    n_train = int(TRAIN_FRACTION * n)
    train_idx = idx_all[:n_train]
    val_idx = idx_all[n_train:]

    if len(train_idx) > max_train_samples:
        train_idx = rng.choice(train_idx, size=max_train_samples, replace=False)

    best = None
    best_val_mape = 1e99

    for strat in STRATEGIES:
        print(f"\n--- Strategy: {strat.name} ---")
        X_bank = build_feature_bank(H, T, d, Uc, L_lin, strat.x_mode)
        y_all = build_target(L_true, L_lin, strat.y_mode)

        X_train_bank = X_bank[train_idx]
        y_train = y_all[train_idx]

        sel3, corr_y = mrmr_select_3(X_train_bank, y_train, alpha=0.5)
        print(
            "Feature relevance |spearman|:",
            ", ".join([f"{FEATURE_NAMES[i]}={corr_y[i]:.3f}" for i in range(5)]),
        )
        print("Selected 3 features:", [FEATURE_NAMES[i] for i in sel3])

        X3_train = X_train_bank[:, sel3]
        X3_val = X_bank[val_idx][:, sel3]
        L_val = L_lin[val_idx]
        L_true_val = L_true[val_idx]

        for deg in DEGREE_GRID:
            p = fit_elm(X3_train, y_train, L_lin[train_idx], L_true[train_idx], strat, deg)
            y_pred_val = predict_y_elm(p, X3_val)
            L_pred_val = reconstruct_L(L_val, y_pred_val, strat)
            m = metrics(L_true_val, L_pred_val)

            hidden = int(HIDDEN_MULT * deg)
            print(
                f"  Size={deg} (hidden={hidden}): Val MAPE={m['mape_pct']:.3f}%  "
                f"(p95={m['p95_ape_pct']:.2f}%, max={m['max_ape_pct']:.2f}%)"
            )

            if m["mape_pct"] < best_val_mape:
                best_val_mape = m["mape_pct"]
                best = (strat, deg, sel3)

    assert best is not None
    strat_best, deg_best, sel3_best = best

    print("\n==========================================")
    print(" BEST MODEL (by validation MAPE)")
    print("==========================================")
    print(f"Strategy : {strat_best.name}")
    print(f"Size     : {deg_best} (hidden={int(HIDDEN_MULT * deg_best)})")
    print("Features :", [FEATURE_NAMES[i] for i in sel3_best])
    print(f"Val MAPE : {best_val_mape:.3f}%")
    print("==========================================\n")

    X_bank = build_feature_bank(H, T, d, Uc, L_lin, strat_best.x_mode)
    y_all = build_target(L_true, L_lin, strat_best.y_mode)
    X3 = X_bank[:, sel3_best]

    print("--- Final training on ALL data ---")
    p = fit_elm(X3, y_all, L_lin, L_true, strat_best, deg_best)

    y_pred = predict_y_elm(p, X3)
    L_pred = reconstruct_L(L_lin, y_pred, strat_best)

    m = metrics(L_true, L_pred)
    print("==========================================")
    print(" FINAL STATISTICS (all data)")
    print("==========================================")
    print(f"Global MSE      : {m['mse']:.6f}")
    print(f"Global MAPE     : {m['mape_pct']:.3f}%")
    print(f"95% APE         : {m['p95_ape_pct']:.2f}%")
    print(f"Max APE         : {m['max_ape_pct']:.2f}%")
    print("==========================================")

    out_df = df.copy()
    out_df["L_lin"] = L_lin
    out_df["L_true"] = L_true
    out_df["L_pred"] = L_pred
    out_df["APE_pct"] = np.abs(out_df["L_true"] - out_df["L_pred"]) / np.maximum(out_df["L_true"], 1e-12) * 100.0

    print("\n--- Worst 20 cases ---")
    worst = out_df.sort_values("APE_pct", ascending=False).head(20)
    cols = ["H", "T", "d", "Uc", "L_lin", "L_true", "L_pred", "APE_pct"]
    print(worst[cols].to_string(index=False, float_format=lambda v: f"{v:.6g}"))

    generate_function_py("function.py", strat_best, sel3_best, p)
    generate_neural_bas("neural.bas", strat_best, sel3_best, p)
    save_plots(out_df, pdf_path="diagnostics.pdf", png_path="plot_all.png")


if __name__ == "__main__":
    main()
