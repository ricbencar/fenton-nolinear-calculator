#!/usr/bin/env python3
"""
High-quality plotting utility for steady nonlinear-wave result files.

The script reads the four files written by the current fourier.cpp solver:
    solution.res   - wave summary, integral quantities and Fourier coefficients
    solver.res     - continuation history, convergence proof and residual vector
    surface.res    - free-surface coordinates and pressure-boundary check
    flowfield.res  - vertical profiles of velocity, acceleration and diagnostics

Typical use:
    python plots.py

Dependencies:
    pip install numpy matplotlib

Each chart is written as an independent figure. By default, both 300 dpi PNG and
vector SVG versions are generated in the output directory.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator


FLOAT_PATTERN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
FLOAT_RE = re.compile(FLOAT_PATTERN)
NUMBER_PATTERN = rf"(?:{FLOAT_PATTERN}|[+-]?(?:inf(?:inity)?|nan))"
PROFILE_RE = re.compile(
    rf"^#\s*X/d\s*=\s*({FLOAT_PATTERN})\s*,\s*Phase\s*=\s*({FLOAT_PATTERN})",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class IntegralQuantity:
    name: str
    symbol: str
    k_based: float
    d_based: float


@dataclass(frozen=True)
class SolutionData:
    title: str
    method: str
    height_depth: float | None
    maximum_height_depth: float | None
    fraction_of_maximum: float | None
    length_depth: float | None
    period_depth: float | None
    current_criterion: str | None
    current_value: float | None
    stokes_ursell: float | None
    quantities: tuple[IntegralQuantity, ...]
    modes: np.ndarray
    stream_coefficients: np.ndarray
    surface_coefficients: np.ndarray
    fourier_order: int | None
    residual_rms: float | None
    residual_max: float | None




@dataclass(frozen=True)
class ContinuationStep:
    trial: int
    fraction: float
    increment: float
    accepted: bool
    residual_rms: float
    residual_max: float
    next_increment: float
    tangent_available: bool


@dataclass(frozen=True)
class ResidualGroup:
    name: str
    first_index: int
    last_index: int
    count: int
    rms: float
    maximum: float
    maximum_index: int


@dataclass(frozen=True)
class ResidualEntry:
    index: int
    value: float
    absolute_value: float
    fraction_of_max_limit: float
    description: str


@dataclass(frozen=True)
class SolverData:
    fourier_order: int
    unknown_count: int
    equation_count: int
    intermediate_tolerance: float
    residual_rms_tolerance: float
    residual_max_tolerance: float
    spectral_defect_tolerance: float
    target_height_depth: float
    current_definition: str
    prescribed_current: float
    continuation: tuple[ContinuationStep, ...]
    accepted_continuation_steps: int
    total_continuation_trials: int
    final_sum_squares: float
    final_l2_norm: float
    final_residual_rms: float
    final_residual_max: float
    final_newton_correction_rms: float
    final_linear_solver: str
    spectral_defect_max: float
    residual_test_passed: bool
    spectral_test_passed: bool
    overall_convergence_passed: bool
    residual_groups: tuple[ResidualGroup, ...]
    residuals: tuple[ResidualEntry, ...]
    state_count: int

@dataclass(frozen=True)
class SurfaceData:
    x_over_d: np.ndarray
    y_over_d: np.ndarray
    pressure_check: np.ndarray


@dataclass(frozen=True)
class FlowProfile:
    x_over_d: float
    phase_deg: float
    values: np.ndarray


@dataclass(frozen=True)
class FlowFieldData:
    profiles: tuple[FlowProfile, ...]
    x: np.ndarray
    y: np.ndarray
    values: np.ndarray
    triangles: np.ndarray
    profile_index: np.ndarray
    vertical_index: np.ndarray


FLOW_COLUMNS = (
    "y_over_d",
    "u_over_sqrt_gd",
    "v_over_sqrt_gd",
    "dphi_dt_over_gd",
    "du_dt_over_g",
    "dv_dt_over_g",
    "du_dx_over_sqrt_g_over_d",
    "du_dy_over_sqrt_g_over_d",
    "bernoulli_check_over_gd",
)


# -----------------------------------------------------------------------------
# Parsing
# -----------------------------------------------------------------------------


def _as_float(text: str) -> float:
    """Read standard, Fortran-D, or Fortran-d exponent notation."""
    return float(text.replace("D", "E").replace("d", "e"))


def _first_float(text: str) -> float | None:
    match = FLOAT_RE.search(text)
    return _as_float(match.group(0)) if match else None


def _all_floats(text: str) -> list[float]:
    return [_as_float(value) for value in FLOAT_RE.findall(text)]


def _strip_comment_prefix(line: str) -> str:
    return line.lstrip().removeprefix("#").strip()


def resolve_input_path(path_text: str, expected_stem: str) -> Path:
    """
    Resolve an input path robustly.

    The exact user-supplied path is preferred. If the default name is absent,
    the current directory is searched case-insensitively for a single matching
    result such as ``solution(2).res``.
    """
    path = Path(path_text).expanduser()
    if path.is_file():
        return path.resolve()

    parent = path.parent if path.parent != Path("") else Path.cwd()
    if not parent.exists():
        raise FileNotFoundError(f"Input directory does not exist: {parent}")

    expected = expected_stem.casefold()
    candidates = sorted(
        candidate
        for candidate in parent.iterdir()
        if candidate.is_file()
        and candidate.suffix.casefold() == ".res"
        and candidate.stem.casefold().startswith(expected)
    )

    if len(candidates) == 1:
        return candidates[0].resolve()

    if not candidates:
        raise FileNotFoundError(f"Required input file not found: {path}")

    joined = "\n  ".join(str(candidate) for candidate in candidates)
    raise FileNotFoundError(
        f"Input file {path} is ambiguous. Specify one of:\n  {joined}"
    )


def parse_solution(path: Path) -> SolutionData:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    comment_lines = [_strip_comment_prefix(line) for line in lines if line.lstrip().startswith("#")]
    title = comment_lines[0] if comment_lines else path.stem
    method = comment_lines[1] if len(comment_lines) > 1 else "Steady-wave solution"

    height_depth = None
    maximum_height_depth = None
    fraction_of_maximum = None
    length_depth = None
    period_depth = None
    current_criterion = None
    current_value = None
    stokes_ursell = None
    fourier_order = None
    residual_rms = None
    residual_max = None

    for raw_line in lines:
        line = _strip_comment_prefix(raw_line)
        lower = line.casefold()

        if lower.startswith("height/depth:"):
            values = _all_floats(line)
            if values:
                height_depth = values[0]
            percent_match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*%", line)
            if percent_match:
                fraction_of_maximum = _as_float(percent_match.group(1)) / 100.0
            max_match = re.search(
                rf"maximum\s+of\s+H/d\s*=\s*({FLOAT_PATTERN})",
                line,
                re.IGNORECASE,
            )
            if max_match:
                maximum_height_depth = _as_float(max_match.group(1))

        elif lower.startswith("length/depth:"):
            length_depth = _first_float(line)

        elif lower.startswith("dimensionless period"):
            period_depth = _first_float(line)

        elif lower.startswith("current criterion:"):
            match = re.search(
                rf"Current\s+criterion:\s*([^,]+),\s*Dimensionless\s+value:\s*({FLOAT_PATTERN})",
                line,
                re.IGNORECASE,
            )
            if match:
                current_criterion = match.group(1).strip()
                current_value = _as_float(match.group(2))

        elif lower.startswith("stokes-ursell number"):
            stokes_ursell = _first_float(line)

        elif lower.startswith("final fourier order:"):
            value = _first_float(line)
            fourier_order = int(value) if value is not None else None

        elif lower.startswith("collocation residual rms:"):
            residual_rms = _first_float(line)

        elif lower.startswith("collocation residual max:"):
            residual_max = _first_float(line)

    quantities: list[IntegralQuantity] = []
    quantity_re = re.compile(
        rf"^#\s*(.*?)\s+\(([^()]*)\)\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s*$"
    )
    for line in lines:
        match = quantity_re.match(line.rstrip())
        if not match:
            continue
        name = match.group(1).strip()
        symbol = match.group(2).strip()
        quantities.append(
            IntegralQuantity(
                name=name,
                symbol=symbol,
                k_based=_as_float(match.group(3)),
                d_based=_as_float(match.group(4)),
            )
        )

    modes: list[int] = []
    stream_coefficients: list[float] = []
    surface_coefficients: list[float] = []
    coefficient_re = re.compile(
        rf"^\s*(\d+)\s+({FLOAT_PATTERN})\s+({FLOAT_PATTERN})\s*$"
    )
    for line in lines:
        match = coefficient_re.match(line)
        if match:
            modes.append(int(match.group(1)))
            stream_coefficients.append(_as_float(match.group(2)))
            surface_coefficients.append(_as_float(match.group(3)))

    if not quantities:
        raise ValueError(f"No integral quantities were found in {path}")
    if not modes:
        raise ValueError(f"No Fourier coefficients were found in {path}")

    return SolutionData(
        title=title,
        method=method,
        height_depth=height_depth,
        maximum_height_depth=maximum_height_depth,
        fraction_of_maximum=fraction_of_maximum,
        length_depth=length_depth,
        period_depth=period_depth,
        current_criterion=current_criterion,
        current_value=current_value,
        stokes_ursell=stokes_ursell,
        quantities=tuple(quantities),
        modes=np.asarray(modes, dtype=int),
        stream_coefficients=np.asarray(stream_coefficients, dtype=float),
        surface_coefficients=np.asarray(surface_coefficients, dtype=float),
        fourier_order=fourier_order if fourier_order is not None else len(modes),
        residual_rms=residual_rms,
        residual_max=residual_max,
    )


def _value_after_prefix(lines: Sequence[str], prefix: str) -> str:
    """Return the text following one case-insensitive solver-report prefix."""
    expected = prefix.casefold()
    for raw_line in lines:
        line = _strip_comment_prefix(raw_line)
        if line.casefold().startswith(expected):
            return line[len(prefix) :].strip()
    raise ValueError(f"Required solver field was not found: {prefix}")


def _solver_float(lines: Sequence[str], prefix: str) -> float:
    expected = prefix.casefold()
    for raw_line in lines:
        line = _strip_comment_prefix(raw_line)
        if not line.casefold().startswith(expected):
            continue
        remainder = line[len(prefix) :].strip()
        match = re.match(NUMBER_PATTERN, remainder, re.IGNORECASE)
        if match:
            return _as_float(match.group(0))
    raise ValueError(f"Required numeric solver field was not found: {prefix}")


def _solver_int(lines: Sequence[str], prefix: str) -> int:
    return int(_solver_float(lines, prefix))


def _solver_pass(lines: Sequence[str], prefix: str) -> bool:
    text = _value_after_prefix(lines, prefix).casefold()
    if text.startswith("pass"):
        return True
    if text.startswith("fail"):
        return False
    raise ValueError(f"Solver field {prefix!r} is neither PASS nor FAIL")


def parse_solver(path: Path) -> SolverData:
    """
    Parse the current solver.res audit report.

    The parser deliberately uses the report labels written by fourier.cpp rather
    than fixed column positions for the configuration and convergence proof.
    Continuation, residual-group and residual-vector rows have stable tabular
    grammars and are parsed with dedicated regular expressions.
    """
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()

    continuation_re = re.compile(
        rf"^\s*(\d+)\s+({NUMBER_PATTERN})\s+({NUMBER_PATTERN})\s+"
        rf"(ACCEPT|REJECT)\s+({NUMBER_PATTERN})\s+({NUMBER_PATTERN})\s+"
        rf"({NUMBER_PATTERN})\s+(YES|NO)\s*$",
        re.IGNORECASE,
    )
    continuation: list[ContinuationStep] = []
    for raw_line in lines:
        match = continuation_re.match(raw_line)
        if not match:
            continue
        continuation.append(
            ContinuationStep(
                trial=int(match.group(1)),
                fraction=_as_float(match.group(2)),
                increment=_as_float(match.group(3)),
                accepted=match.group(4).casefold() == "accept",
                residual_rms=_as_float(match.group(5)),
                residual_max=_as_float(match.group(6)),
                next_increment=_as_float(match.group(7)),
                tangent_available=match.group(8).casefold() == "yes",
            )
        )

    group_re = re.compile(
        rf"^(.*?)\s+F\[\s*(\d+)\.\.\s*(\d+)\]\s+count=\s*(\d+)\s+"
        rf"RMS=\s*({NUMBER_PATTERN})\s+MAX=\s*({NUMBER_PATTERN})\s+"
        rf"at\s+F\[\s*(\d+)\]\s*$",
        re.IGNORECASE,
    )
    residual_groups: list[ResidualGroup] = []
    for raw_line in lines:
        match = group_re.match(raw_line)
        if not match:
            continue
        residual_groups.append(
            ResidualGroup(
                name=match.group(1).strip(),
                first_index=int(match.group(2)),
                last_index=int(match.group(3)),
                count=int(match.group(4)),
                rms=_as_float(match.group(5)),
                maximum=_as_float(match.group(6)),
                maximum_index=int(match.group(7)),
            )
        )

    residual_re = re.compile(
        rf"^F\[\s*(\d+)\]\s*=\s*({NUMBER_PATTERN})\s+"
        rf"({NUMBER_PATTERN})\s+({NUMBER_PATTERN})\s+(.*?)\s*$",
        re.IGNORECASE,
    )
    residuals: list[ResidualEntry] = []
    state_indices: list[int] = []
    state_re = re.compile(r"^z\[\s*(\d+)\]\s*=")
    for raw_line in lines:
        residual_match = residual_re.match(raw_line)
        if residual_match:
            residuals.append(
                ResidualEntry(
                    index=int(residual_match.group(1)),
                    value=_as_float(residual_match.group(2)),
                    absolute_value=_as_float(residual_match.group(3)),
                    fraction_of_max_limit=_as_float(residual_match.group(4)),
                    description=residual_match.group(5).strip(),
                )
            )
            continue
        state_match = state_re.match(raw_line)
        if state_match:
            state_indices.append(int(state_match.group(1)))

    if not continuation:
        raise ValueError(f"No adaptive-continuation rows were found in {path}")
    if not residual_groups:
        raise ValueError(f"No residual-group summary was found in {path}")
    if not residuals:
        raise ValueError(f"No complete residual vector was found in {path}")

    fourier_order = _solver_int(lines, "Fourier order N")
    unknown_count = _solver_int(lines, "Nonlinear unknowns")
    equation_count = _solver_int(lines, "Nonlinear equations")
    accepted_steps = _solver_int(lines, "Accepted continuation steps:")
    total_trials = _solver_int(lines, "Total continuation trials:")

    if [entry.index for entry in residuals] != list(range(1, equation_count + 1)):
        raise ValueError(
            f"Residual indices in {path} are incomplete or out of order; "
            f"expected F[1..{equation_count}]"
        )
    if state_indices != list(range(1, unknown_count + 1)):
        raise ValueError(
            f"State indices in {path} are incomplete or out of order; "
            f"expected z[1..{unknown_count}]"
        )

    return SolverData(
        fourier_order=fourier_order,
        unknown_count=unknown_count,
        equation_count=equation_count,
        intermediate_tolerance=_solver_float(
            lines, "Intermediate residual tolerance"
        ),
        residual_rms_tolerance=_solver_float(
            lines, "Final residual RMS tolerance"
        ),
        residual_max_tolerance=_solver_float(
            lines, "Final residual max tolerance"
        ),
        spectral_defect_tolerance=_solver_float(
            lines, "Off-grid spectral defect tolerance"
        ),
        target_height_depth=_solver_float(lines, "Target H/d"),
        current_definition=_value_after_prefix(lines, "Current definition"),
        prescribed_current=_solver_float(
            lines, "Prescribed current / sqrt(gd)"
        ),
        continuation=tuple(continuation),
        accepted_continuation_steps=accepted_steps,
        total_continuation_trials=total_trials,
        final_sum_squares=_solver_float(lines, "Final residual sum of squares"),
        final_l2_norm=_solver_float(lines, "Final residual L2 norm"),
        final_residual_rms=_solver_float(lines, "Final residual RMS"),
        final_residual_max=_solver_float(
            lines, "Final maximum absolute residual"
        ),
        final_newton_correction_rms=_solver_float(
            lines, "Final scaled Newton correction RMS"
        ),
        final_linear_solver=_value_after_prefix(lines, "Final linear solver"),
        spectral_defect_max=_solver_float(
            lines, "Maximum between-node spectral defect"
        ),
        residual_test_passed=_solver_pass(
            lines, "Combined residual convergence test"
        ),
        spectral_test_passed=_solver_pass(lines, "Spectral resolution check"),
        overall_convergence_passed=_solver_pass(
            lines, "Overall nonlinear convergence"
        ),
        residual_groups=tuple(residual_groups),
        residuals=tuple(residuals),
        state_count=len(state_indices),
    )


def merge_solver_metadata(
    solution: SolutionData, solver: SolverData
) -> SolutionData:
    """Use solver.res as the authoritative source for numerical diagnostics."""
    return replace(
        solution,
        fourier_order=solver.fourier_order,
        residual_rms=solver.final_residual_rms,
        residual_max=solver.final_residual_max,
    )


def validate_result_set(
    solution: SolutionData,
    solver: SolverData,
    surface: SurfaceData,
    flow: FlowFieldData,
) -> None:
    """Reject mixed, incomplete or internally inconsistent result files."""
    if len(solution.modes) != solver.fourier_order:
        raise ValueError(
            "solution.res and solver.res disagree on Fourier order: "
            f"{len(solution.modes)} coefficient rows versus N={solver.fourier_order}"
        )
    if len(solver.continuation) != solver.total_continuation_trials:
        raise ValueError(
            "solver.res continuation row count does not match its reported trial count"
        )
    if sum(step.accepted for step in solver.continuation) != solver.accepted_continuation_steps:
        raise ValueError(
            "solver.res accepted continuation count does not match its history"
        )
    if solver.unknown_count != 2 * solver.fourier_order + 10:
        raise ValueError(
            "solver.res has an inconsistent nonlinear unknown count: "
            f"{solver.unknown_count} for N={solver.fourier_order}"
        )
    if solver.equation_count != solver.unknown_count:
        raise ValueError("solver.res does not describe a square nonlinear system")
    if not math.isclose(
        solver.target_height_depth,
        solution.height_depth if solution.height_depth is not None else math.nan,
        rel_tol=5.0e-5,
        abs_tol=5.0e-7,
    ):
        raise ValueError("solution.res and solver.res describe different H/d cases")
    if solution.current_criterion and (
        solution.current_criterion.casefold()
        != solver.current_definition.casefold()
    ):
        raise ValueError(
            "solution.res and solver.res use different current definitions"
        )
    if solution.current_value is not None and not math.isclose(
        solution.current_value,
        solver.prescribed_current,
        rel_tol=5.0e-5,
        abs_tol=5.0e-7,
    ):
        raise ValueError("solution.res and solver.res use different currents")
    if len(surface.x_over_d) < 3:
        raise ValueError("surface.res must contain at least three physical points")
    if not np.all(np.diff(surface.x_over_d) > 0.0):
        raise ValueError("surface.res X/d coordinates must be strictly increasing")
    if not np.all(np.isfinite(surface.y_over_d)):
        raise ValueError("surface.res contains non-finite free-surface elevations")
    if not np.all(np.isfinite(flow.values)):
        raise ValueError("flowfield.res contains non-finite numerical values")
    phases = np.asarray([profile.phase_deg for profile in flow.profiles])
    if not np.all(np.diff(phases) > 0.0):
        raise ValueError("flowfield.res phases must be strictly increasing")
    if abs(phases[0]) > 1.0e-8 or abs(phases[-1] - 180.0) > 1.0e-6:
        raise ValueError("flowfield.res must span the crest-to-trough half-wave")


def parse_surface(path: Path) -> SurfaceData:
    rows: list[list[float]] = []
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        numeric_part = raw_line.split("#", 1)[0]
        values = _all_floats(numeric_part)
        if len(values) >= 3:
            # The historical files may contain a dummy 0,0,0 scaling point.
            if "dummy" in raw_line.casefold():
                continue
            rows.append(values[:3])

    if not rows:
        raise ValueError(f"No surface points were found in {path}")

    array = np.asarray(rows, dtype=float)
    order = np.argsort(array[:, 0])
    array = array[order]

    return SurfaceData(
        x_over_d=array[:, 0],
        y_over_d=array[:, 1],
        pressure_check=array[:, 2],
    )


def parse_flowfield(path: Path) -> FlowFieldData:
    profiles: list[FlowProfile] = []
    current_x: float | None = None
    current_phase: float | None = None
    current_rows: list[list[float]] = []

    def commit_profile() -> None:
        nonlocal current_rows
        if current_x is None or current_phase is None:
            return
        if not current_rows:
            raise ValueError(
                f"Flow profile at X/d={current_x:g} contains no numerical rows"
            )
        values = np.asarray(current_rows, dtype=float)
        if values.shape[1] != len(FLOW_COLUMNS):
            raise ValueError(
                f"Flow profile at X/d={current_x:g} has {values.shape[1]} columns; "
                f"expected {len(FLOW_COLUMNS)}"
            )
        profiles.append(FlowProfile(current_x, current_phase, values))
        current_rows = []

    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        profile_match = PROFILE_RE.match(raw_line.strip())
        if profile_match:
            commit_profile()
            current_x = _as_float(profile_match.group(1))
            current_phase = _as_float(profile_match.group(2))
            continue

        stripped = raw_line.strip()
        if current_x is None or not stripped or stripped.startswith("#"):
            continue

        values = _all_floats(stripped)
        if len(values) >= len(FLOW_COLUMNS):
            current_rows.append(values[: len(FLOW_COLUMNS)])

    commit_profile()

    if len(profiles) < 2:
        raise ValueError(f"At least two flow profiles are required in {path}")

    profiles.sort(key=lambda profile: profile.x_over_d)
    vertical_counts = {profile.values.shape[0] for profile in profiles}
    if len(vertical_counts) != 1:
        counts = ", ".join(str(value) for value in sorted(vertical_counts))
        raise ValueError(
            "All flow profiles must have the same vertical point count; "
            f"found {counts}"
        )

    n_profiles = len(profiles)
    n_vertical = profiles[0].values.shape[0]

    x_values: list[float] = []
    y_values: list[float] = []
    all_values: list[np.ndarray] = []
    profile_indices: list[int] = []
    vertical_indices: list[int] = []

    for profile_number, profile in enumerate(profiles):
        for vertical_number, row in enumerate(profile.values):
            x_values.append(profile.x_over_d)
            y_values.append(row[0])
            all_values.append(row)
            profile_indices.append(profile_number)
            vertical_indices.append(vertical_number)

    triangles: list[tuple[int, int, int]] = []
    for i in range(n_profiles - 1):
        for j in range(n_vertical - 1):
            lower_left = i * n_vertical + j
            upper_left = lower_left + 1
            lower_right = (i + 1) * n_vertical + j
            upper_right = lower_right + 1
            triangles.append((lower_left, lower_right, upper_right))
            triangles.append((lower_left, upper_right, upper_left))

    return FlowFieldData(
        profiles=tuple(profiles),
        x=np.asarray(x_values, dtype=float),
        y=np.asarray(y_values, dtype=float),
        values=np.asarray(all_values, dtype=float),
        triangles=np.asarray(triangles, dtype=int),
        profile_index=np.asarray(profile_indices, dtype=int),
        vertical_index=np.asarray(vertical_indices, dtype=int),
    )


# -----------------------------------------------------------------------------
# Plot presentation
# -----------------------------------------------------------------------------


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "figure.figsize": (11.5, 7.2),
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "font.size": 11,
            "axes.titlesize": 16,
            "axes.labelsize": 12,
            "legend.fontsize": 10,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "axes.axisbelow": True,
            "lines.linewidth": 2.0,
            "lines.markersize": 5.0,
            "figure.constrained_layout.use": False,
        }
    )


def format_optional(value: float | None, spec: str = ".5g") -> str:
    return "n/a" if value is None else format(value, spec)


def case_caption(solution: SolutionData) -> str:
    parts = [solution.title]
    if solution.height_depth is not None:
        parts.append(f"H/d={solution.height_depth:.5g}")
    if solution.length_depth is not None:
        parts.append(f"L/d={solution.length_depth:.6g}")
    if solution.period_depth is not None:
        parts.append(f"T√(g/d)={solution.period_depth:.6g}")
    if solution.current_criterion and solution.current_value is not None:
        parts.append(f"{solution.current_criterion} current={solution.current_value:.6g}√(gd)")
    if solution.stokes_ursell is not None:
        parts.append(f"SU={solution.stokes_ursell:.5g}")
    if solution.fourier_order is not None:
        parts.append(f"N={solution.fourier_order}")
    return "  |  ".join(parts)


def add_case_footer(fig: Figure, solution: SolutionData) -> None:
    fig.text(0.5, 0.012, case_caption(solution), ha="center", va="bottom", fontsize=9)


def add_stat_box(ax: Axes, lines: Sequence[str], location: str = "upper right") -> None:
    anchors = {
        "upper right": (0.985, 0.975, "right", "top"),
        "upper left": (0.015, 0.975, "left", "top"),
        "lower right": (0.985, 0.025, "right", "bottom"),
        "lower left": (0.015, 0.025, "left", "bottom"),
    }
    x, y, horizontal, vertical = anchors[location]
    ax.text(
        x,
        y,
        "\n".join(lines),
        transform=ax.transAxes,
        ha=horizontal,
        va=vertical,
        fontsize=9.5,
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "white", "alpha": 0.88},
    )


def finish_figure(fig: Figure, solution: SolutionData, bottom: float = 0.075) -> None:
    add_case_footer(fig, solution)
    fig.tight_layout(rect=(0.02, bottom, 0.99, 0.98))


def save_figure(
    fig: Figure,
    output_dir: Path,
    stem: str,
    formats: Sequence[str],
    dpi: int,
    *,
    close: bool = True,
) -> list[Path]:
    output_paths: list[Path] = []
    for extension in formats:
        output_path = output_dir / f"{stem}.{extension}"
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        output_paths.append(output_path)
    if close:
        plt.close(fig)
    return output_paths


def safe_log_values(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    finite_nonzero = np.abs(array[np.isfinite(array) & (array != 0.0)])
    floor = (
        np.finfo(float).tiny
        if finite_nonzero.size == 0
        else max(float(np.min(finite_nonzero)) * 0.1, np.finfo(float).tiny)
    )
    result = np.full(array.shape, np.nan, dtype=float)
    finite = np.isfinite(array)
    result[finite] = np.maximum(np.abs(array[finite]), floor)
    return result


# -----------------------------------------------------------------------------
# Solution and surface plots
# -----------------------------------------------------------------------------


def plot_solution_summary(solution: SolutionData) -> Figure:
    rows = [
        [
            quantity.name,
            quantity.symbol,
            f"{quantity.k_based:.10g}",
            f"{quantity.d_based:.10g}",
        ]
        for quantity in solution.quantities
    ]

    figure_height = max(9.5, 0.39 * len(rows) + 4.8)
    fig, ax = plt.subplots(figsize=(12.5, figure_height))
    ax.axis("off")
    ax.set_title("Steady-wave solution summary", pad=24, fontweight="bold")

    maximum_percentage = (
        None
        if solution.fraction_of_maximum is None
        else 100.0 * solution.fraction_of_maximum
    )
    metadata = [
        f"Method: {solution.method}",
        f"H/d: {format_optional(solution.height_depth, '.8g')}",
        f"Maximum H/d: {format_optional(solution.maximum_height_depth, '.8g')}",
        "Fraction of theoretical maximum: "
        f"{format_optional(maximum_percentage, '.5g')}%",
        f"L/d: {format_optional(solution.length_depth, '.9g')}",
        f"T√(g/d): {format_optional(solution.period_depth, '.9g')}",
        f"Current: {solution.current_criterion or 'n/a'}; "
        f"value={format_optional(solution.current_value, '.9g')}√(gd)",
        f"Stokes–Ursell number: {format_optional(solution.stokes_ursell, '.8g')}",
        "Fourier order: "
        f"{solution.fourier_order if solution.fourier_order is not None else len(solution.modes)}",
        f"Residual RMS: {format_optional(solution.residual_rms, '.5e')}",
        f"Residual max: {format_optional(solution.residual_max, '.5e')}",
    ]
    ax.text(
        0.0,
        0.995,
        "\n".join(metadata),
        transform=ax.transAxes,
        ha="left",
        va="top",
        linespacing=1.38,
    )

    table = ax.table(
        cellText=rows,
        colLabels=["Quantity", "Symbol", "k-based value", "d-based value"],
        colLoc="center",
        cellLoc="right",
        colWidths=[0.48, 0.10, 0.21, 0.21],
        bbox=(0.0, 0.015, 1.0, 0.69),
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9.3)
    table.scale(1.0, 1.23)
    for (row_index, column_index), cell in table.get_celld().items():
        if row_index == 0:
            cell.set_text_props(fontweight="bold", ha="center")
        elif column_index == 0:
            cell.set_text_props(ha="left")
        elif column_index == 1:
            cell.set_text_props(ha="center")

    add_case_footer(fig, solution)
    fig.tight_layout(rect=(0.025, 0.04, 0.975, 0.985))
    return fig


def plot_solver_summary(solution: SolutionData, solver: SolverData) -> Figure:
    """Create a compact audit page from the authoritative solver.res values."""
    fig, ax = plt.subplots(figsize=(12.5, 8.8))
    ax.axis("off")
    ax.set_title("Nonlinear solver verification summary", pad=24, fontweight="bold")

    checks = [
        ("Residual convergence", solver.residual_test_passed),
        ("Spectral resolution", solver.spectral_test_passed),
        ("Overall convergence", solver.overall_convergence_passed),
    ]
    metrics = [
        ("Fourier order", f"{solver.fourier_order}"),
        ("Unknowns / equations", f"{solver.unknown_count} / {solver.equation_count}"),
        (
            "Continuation steps / trials",
            f"{solver.accepted_continuation_steps} / {solver.total_continuation_trials}",
        ),
        ("Final residual RMS", f"{solver.final_residual_rms:.6e}"),
        ("RMS tolerance", f"{solver.residual_rms_tolerance:.6e}"),
        ("Final maximum residual", f"{solver.final_residual_max:.6e}"),
        ("Maximum-residual tolerance", f"{solver.residual_max_tolerance:.6e}"),
        ("Final Newton correction RMS", f"{solver.final_newton_correction_rms:.6e}"),
        ("Final linear solver", solver.final_linear_solver),
        ("Maximum off-grid defect", f"{solver.spectral_defect_max:.6e}"),
        ("Off-grid defect tolerance", f"{solver.spectral_defect_tolerance:.6e}"),
    ]

    ax.text(
        0.0,
        0.97,
        "\n".join(f"{label}: {'PASS' if passed else 'FAIL'}" for label, passed in checks),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        linespacing=1.5,
        bbox={"boxstyle": "round,pad=0.55", "facecolor": "white", "alpha": 0.9},
    )

    table = ax.table(
        cellText=[[name, value] for name, value in metrics],
        colLabels=["Verification quantity", "Reported value"],
        colLoc="center",
        cellLoc="left",
        colWidths=[0.58, 0.42],
        bbox=(0.0, 0.12, 1.0, 0.67),
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10.5)
    table.scale(1.0, 1.35)
    for (row, column), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(fontweight="bold", ha="center")
        elif column == 1:
            cell.set_text_props(ha="right")

    ax.text(
        0.0,
        0.045,
        "Residual convergence and spectral resolution are independent checks. "
        "The spectral tolerance is not inferred from the collocation residuals.",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
    )
    add_case_footer(fig, solution)
    fig.tight_layout(rect=(0.025, 0.05, 0.975, 0.985))
    return fig


def plot_continuation_history(solution: SolutionData, solver: SolverData) -> Figure:
    """Plot residual norms at the end of every continuation trial."""
    fig, ax = plt.subplots(figsize=(11.8, 7.2))
    trials = np.asarray([step.trial for step in solver.continuation], dtype=int)
    rms = np.asarray([step.residual_rms for step in solver.continuation], dtype=float)
    maximum = np.asarray([step.residual_max for step in solver.continuation], dtype=float)
    accepted = np.asarray([step.accepted for step in solver.continuation], dtype=bool)

    finite_rms = np.where(np.isfinite(rms), rms, np.nan)
    finite_max = np.where(np.isfinite(maximum), maximum, np.nan)
    ax.semilogy(trials, safe_log_values(finite_rms), marker="o", label="Residual RMS")
    ax.semilogy(trials, safe_log_values(finite_max), marker="s", label="Maximum residual")
    ax.axhline(
        solver.intermediate_tolerance,
        linestyle="--",
        linewidth=1.3,
        label="Intermediate RMS tolerance",
    )
    maximum_factor = (
        solver.residual_max_tolerance / solver.residual_rms_tolerance
    )
    ax.axhline(
        maximum_factor * solver.intermediate_tolerance,
        linestyle=":",
        linewidth=1.3,
        label="Intermediate maximum tolerance",
    )

    rejected_trials = trials[~accepted]
    if rejected_trials.size:
        y_reference = np.nanmax(np.concatenate((finite_rms, finite_max)))
        ax.scatter(
            rejected_trials,
            np.full(rejected_trials.shape, y_reference),
            marker="x",
            s=70,
            label="Rejected continuation trial",
            zorder=5,
        )

    add_stat_box(
        ax,
        [
            f"Accepted steps: {solver.accepted_continuation_steps}",
            f"Total trials: {solver.total_continuation_trials}",
            "Rejected trials: "
            f"{solver.total_continuation_trials - solver.accepted_continuation_steps}",
            f"Final target fraction: {solver.continuation[-1].fraction:.12g}",
        ],
        "upper right",
    )
    ax.set_title("Adaptive continuation convergence history", fontweight="bold")
    ax.set_xlabel("Continuation trial")
    ax.set_ylabel("Residual norm")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.legend(loc="best")
    finish_figure(fig, solution)
    return fig


def plot_final_residual_vector(solution: SolutionData, solver: SolverData) -> Figure:
    """Plot every final residual relative to its configured maximum limit."""
    fig, ax = plt.subplots(figsize=(12.0, 7.0))
    indices = np.asarray([entry.index for entry in solver.residuals], dtype=int)
    fractions = np.asarray(
        [entry.fraction_of_max_limit for entry in solver.residuals], dtype=float
    )
    display_values = safe_log_values(fractions)

    ax.semilogy(indices, display_values, marker="o", markersize=3.5, label="|F[i]| / maximum limit")
    ax.axhline(1.0, linestyle="--", linewidth=1.4, label="Configured limit")

    for group in solver.residual_groups:
        if group.first_index > 1:
            ax.axvline(group.first_index - 0.5, linewidth=0.9, linestyle=":")
        center = 0.5 * (group.first_index + group.last_index)
        ax.text(
            center,
            0.98,
            group.name.replace(" surface equations", ""),
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=9,
        )

    maximum_entry = max(solver.residuals, key=lambda entry: entry.absolute_value)
    add_stat_box(
        ax,
        [
            f"Residual equations: {solver.equation_count}",
            f"Largest |F|: {maximum_entry.absolute_value:.3e} at F[{maximum_entry.index}]",
            f"Largest fraction of limit: {np.max(fractions):.3e}",
            f"Residual test: {'PASS' if solver.residual_test_passed else 'FAIL'}",
        ],
        "lower right",
    )
    ax.set_title("Complete final nonlinear residual vector", fontweight="bold")
    ax.set_xlabel("Residual equation index, i")
    ax.set_ylabel("Absolute residual divided by configured maximum limit")
    ax.set_xlim(1, solver.equation_count)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=12))
    ax.legend(loc="lower left")
    finish_figure(fig, solution)
    return fig


def plot_wave_profile(solution: SolutionData, surface: SurfaceData) -> Figure:
    fig, ax = plt.subplots()
    x = surface.x_over_d
    eta = surface.y_over_d

    ax.fill_between(x, 0.0, eta, alpha=0.10, label="Water domain")
    ax.plot(x, eta, label="Free surface y/d")
    ax.axhline(1.0, linestyle="--", linewidth=1.3, label="Mean water level y/d = 1")
    ax.axhline(0.0, linewidth=1.6, label="Horizontal bed y/d = 0")

    crest_index = int(np.argmax(eta))
    trough_index = int(np.argmin(eta))
    ax.scatter([x[crest_index], x[trough_index]], [eta[crest_index], eta[trough_index]], zorder=5)
    ax.annotate(
        f"Crest: η/d={eta[crest_index]:.5f}",
        xy=(x[crest_index], eta[crest_index]),
        xytext=(24, -38),
        textcoords="offset points",
        arrowprops={"arrowstyle": "->"},
    )
    ax.annotate(
        f"Trough: η/d={eta[trough_index]:.5f}",
        xy=(x[trough_index], eta[trough_index]),
        xytext=(14, -34),
        textcoords="offset points",
        arrowprops={"arrowstyle": "->"},
    )

    measured_height = float(np.max(eta) - np.min(eta))
    stats = [
        f"Surface points: {len(x)}",
        f"Computed crest-to-trough H/d: {measured_height:.6f}",
        f"Crest level: {np.max(eta):.6f} d",
        f"Trough level: {np.min(eta):.6f} d",
        f"Horizontal range: {np.min(x):.5f} ≤ X/d ≤ {np.max(x):.5f}",
        "Vertical scale exaggerated for clarity",
    ]
    add_stat_box(ax, stats, "upper right")

    ax.set_title("Nonlinear free-surface profile", fontweight="bold")
    ax.set_xlabel("Wave coordinate, X/d")
    ax.set_ylabel("Elevation above bed, y/d")
    ax.set_xlim(float(np.min(x)), float(np.max(x)))
    ax.set_ylim(bottom=-0.04, top=max(1.62, float(np.max(eta)) * 1.07))
    ax.set_box_aspect(0.34)
    ax.legend(loc="lower center", ncols=2)
    finish_figure(fig, solution)
    return fig


def plot_surface_pressure_check(solution: SolutionData, surface: SurfaceData) -> Figure:
    fig, ax = plt.subplots()
    x = surface.x_over_d
    check = surface.pressure_check

    ax.plot(x, check, marker="o", label="Signed surface-pressure check")
    ax.axhline(0.0, linewidth=1.0)
    max_abs = float(np.max(np.abs(check)))
    rms = float(np.sqrt(np.mean(np.square(check))))
    nonzero_abs = np.abs(check[np.nonzero(check)])
    linthresh = max(1.0e-16, float(np.min(nonzero_abs)) if nonzero_abs.size else 1.0e-16)
    ax.set_yscale("symlog", linthresh=linthresh)

    add_stat_box(
        ax,
        [
            f"Maximum |check|: {max_abs:.3e}",
            f"RMS check: {rms:.3e}",
            f"Mean check: {np.mean(check):.3e}",
            "Ideal value: 0",
        ],
        "upper right",
    )

    ax.set_title("Dynamic free-surface pressure consistency", fontweight="bold")
    ax.set_xlabel("Wave coordinate, X/d")
    ax.set_ylabel("Surface-pressure check (dimensionless)")
    ax.legend(loc="best")
    finish_figure(fig, solution)
    return fig


def plot_fourier_coefficients(solution: SolutionData) -> Figure:
    fig, ax = plt.subplots()
    modes = solution.modes
    b_values = safe_log_values(solution.stream_coefficients)
    e_values = safe_log_values(solution.surface_coefficients)

    ax.semilogy(
        modes,
        b_values,
        marker="o",
        markevery=max(1, len(modes) // 20),
        label="|B[j]| potential/stream-function",
    )
    ax.semilogy(
        modes,
        e_values,
        marker="s",
        markevery=max(1, len(modes) // 20),
        label="|E[j]| surface elevation",
    )

    b_peak = float(np.max(b_values))
    e_peak = float(np.max(e_values))
    b_last = float(b_values[-1])
    e_last = float(e_values[-1])
    add_stat_box(
        ax,
        [
            f"Modes represented: {len(modes)}",
            f"|B[N]| / max|B|: {b_last / b_peak:.3e}",
            f"|E[N]| / max|E|: {e_last / e_peak:.3e}",
            f"Reported residual max: {format_optional(solution.residual_max, '.3e')}",
        ],
        "upper right",
    )

    ax.set_title("Fourier-coefficient spectral decay", fontweight="bold")
    ax.set_xlabel("Fourier mode, j")
    ax.set_ylabel("Absolute coefficient magnitude")
    ax.set_xlim(int(np.min(modes)), int(np.max(modes)))
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=11))
    ax.legend(loc="best")
    finish_figure(fig, solution)
    return fig


def plot_depth_scaled_integrals(solution: SolutionData) -> Figure:
    excluded_symbols = {"d", "lambda", "H", "tau", "c", "u1", "u2", "U", "q", "r"}
    quantities = [
        quantity for quantity in solution.quantities if quantity.symbol not in excluded_symbols
    ]

    labels = [f"{quantity.name} ({quantity.symbol})" for quantity in quantities]
    values = np.asarray([quantity.d_based for quantity in quantities], dtype=float)
    order = np.argsort(values)
    labels = [labels[index] for index in order]
    values = values[order]

    fig, ax = plt.subplots(figsize=(12.5, max(7.0, 0.58 * len(labels) + 2.4)))
    positions = np.arange(len(values))
    bars = ax.barh(positions, values)
    ax.set_yticks(positions, labels)
    if np.all(values > 0.0):
        ax.set_xscale("log")
    ax.set_xlabel("Depth-based dimensionless value")
    ax.set_title("Depth-scaled integral quantities and invariants", fontweight="bold")

    for bar, value in zip(bars, values, strict=True):
        ax.annotate(
            f"{value:.7g}",
            xy=(value, bar.get_y() + bar.get_height() / 2.0),
            xytext=(5, 0),
            textcoords="offset points",
            va="center",
            ha="left",
            fontsize=9,
        )

    ax.grid(True, axis="x")
    ax.grid(False, axis="y")
    finish_figure(fig, solution)
    return fig


# -----------------------------------------------------------------------------
# Flow-field plotting
# -----------------------------------------------------------------------------


def triangulation(flow: FlowFieldData) -> mtri.Triangulation:
    return mtri.Triangulation(flow.x, flow.y, flow.triangles)


def flow_surface(flow: FlowFieldData) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray([profile.x_over_d for profile in flow.profiles], dtype=float)
    eta = np.asarray([profile.values[-1, 0] for profile in flow.profiles], dtype=float)
    return x, eta


def add_flow_boundaries(ax: Axes, flow: FlowFieldData) -> None:
    x_surface, y_surface = flow_surface(flow)
    ax.plot(x_surface, y_surface, linewidth=2.2, label="Free surface")
    ax.plot([x_surface[0], x_surface[-1]], [0.0, 0.0], linewidth=1.7, label="Bed")
    ax.set_xlim(float(np.min(x_surface)), float(np.max(x_surface)))
    ax.set_ylim(-0.02, float(np.max(y_surface)) * 1.035)
    ax.set_box_aspect(0.36)


def field_plot(
    solution: SolutionData,
    flow: FlowFieldData,
    field: np.ndarray,
    title: str,
    colorbar_label: str,
    symmetric: bool = False,
    overlay_vectors: tuple[np.ndarray, np.ndarray, str] | None = None,
) -> Figure:
    fig, ax = plt.subplots(figsize=(12.0, 6.4))
    tri = triangulation(flow)
    values = np.asarray(field, dtype=float)

    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise ValueError(f"Field {title!r} contains no finite values")

    minimum = float(np.min(finite))
    maximum = float(np.max(finite))
    span = maximum - minimum

    if symmetric:
        limit = max(abs(minimum), abs(maximum))
        if limit > 0.0:
            levels = np.linspace(-limit, limit, 31)
        else:
            levels = None
    else:
        levels = np.linspace(minimum, maximum, 31) if span > 1.0e-15 else None

    constant_tolerance = max(1.0e-15, 1.0e-12 * max(abs(minimum), abs(maximum), 1.0))
    constant_field = span <= constant_tolerance

    if constant_field:
        artist = ax.tripcolor(tri, values, shading="flat")
        ax.text(
            0.5,
            0.52,
            f"Constant field at file precision: {minimum:.6g}",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=12,
            bbox={"boxstyle": "round,pad=0.45", "facecolor": "white", "alpha": 0.88},
        )
    else:
        artist = ax.tricontourf(tri, values, levels=levels, extend="both")
        colorbar = fig.colorbar(artist, ax=ax, pad=0.025, shrink=0.94)
        colorbar.set_label(colorbar_label)

    if overlay_vectors is not None:
        u_vector, v_vector, key_label = overlay_vectors
        mask = (flow.vertical_index % 2 == 0)
        xq = flow.x[mask]
        yq = flow.y[mask]
        uq = np.asarray(u_vector)[mask]
        vq = np.asarray(v_vector)[mask]
        vector_magnitude = np.hypot(uq, vq)
        max_vector = float(np.max(vector_magnitude))
        if max_vector > 0.0:
            target_arrow_length = 0.28
            scale = max_vector / target_arrow_length
            quiver = ax.quiver(
                xq,
                yq,
                uq,
                vq,
                angles="xy",
                scale_units="xy",
                scale=scale,
                width=0.0024,
            )
            reference = 10.0 ** math.floor(math.log10(max_vector))
            while reference * 2.5 < max_vector:
                reference *= 2.0
            ax.quiverkey(
                quiver,
                X=0.14,
                Y=0.08,
                U=reference,
                label=f"{reference:.3g} {key_label}",
                labelpos="E",
                coordinates="axes",
            )

    add_flow_boundaries(ax, flow)
    add_stat_box(
        ax,
        [
            f"Minimum: {minimum:.6g}",
            f"Maximum: {maximum:.6g}",
            f"RMS: {np.sqrt(np.mean(np.square(finite))):.6g}",
            f"Profiles: {len(flow.profiles)}",
            f"Vertical points/profile: {flow.profiles[0].values.shape[0]}",
            "Vertical scale exaggerated for clarity",
        ],
        "lower right" if overlay_vectors is not None else "upper right",
    )

    ax.set_title(title, fontweight="bold")
    ax.set_xlabel("Wave coordinate from crest to trough, X/d")
    ax.set_ylabel("Elevation above bed, y/d")
    ax.grid(False)
    finish_figure(fig, solution)
    return fig


def plot_flow_velocity_magnitude(solution: SolutionData, flow: FlowFieldData) -> Figure:
    u = flow.values[:, 1]
    v = flow.values[:, 2]
    speed = np.hypot(u, v)
    return field_plot(
        solution,
        flow,
        speed,
        "Velocity magnitude and velocity vectors",
        "|V| / √(gd)",
        symmetric=False,
        overlay_vectors=(u, v, "√(gd)"),
    )


def plot_horizontal_velocity(solution: SolutionData, flow: FlowFieldData) -> Figure:
    return field_plot(
        solution,
        flow,
        flow.values[:, 1],
        "Horizontal velocity field",
        "u / √(gd)",
        symmetric=True,
    )


def plot_vertical_velocity(solution: SolutionData, flow: FlowFieldData) -> Figure:
    return field_plot(
        solution,
        flow,
        flow.values[:, 2],
        "Vertical velocity field",
        "v / √(gd)",
        symmetric=True,
    )


def plot_local_acceleration(solution: SolutionData, flow: FlowFieldData) -> Figure:
    du_dt = flow.values[:, 4]
    dv_dt = flow.values[:, 5]
    magnitude = np.hypot(du_dt, dv_dt)
    return field_plot(
        solution,
        flow,
        magnitude,
        "Local acceleration magnitude and vectors",
        "√[(∂u/∂t)² + (∂v/∂t)²] / g",
        symmetric=False,
        overlay_vectors=(du_dt, dv_dt, "g"),
    )


def plot_potential_time_derivative(solution: SolutionData, flow: FlowFieldData) -> Figure:
    return field_plot(
        solution,
        flow,
        flow.values[:, 3],
        "Velocity-potential time derivative",
        "(∂φ/∂t) / (gd)",
        symmetric=True,
    )


def plot_velocity_gradient_magnitude(solution: SolutionData, flow: FlowFieldData) -> Figure:
    du_dx = flow.values[:, 6]
    du_dy = flow.values[:, 7]
    magnitude = np.hypot(du_dx, du_dy)
    return field_plot(
        solution,
        flow,
        magnitude,
        "Horizontal-velocity gradient magnitude",
        "√[(∂u/∂x)² + (∂u/∂y)²] / √(g/d)",
        symmetric=False,
    )


def plot_bernoulli_check(solution: SolutionData, flow: FlowFieldData) -> Figure:
    check = flow.values[:, 8]
    return field_plot(
        solution,
        flow,
        check,
        "Bernoulli-equation consistency over the flow field",
        "Bernoulli check / (gd)",
        symmetric=True,
    )


def select_profiles(flow: FlowFieldData, target_phases: Iterable[float]) -> list[FlowProfile]:
    selected: list[FlowProfile] = []
    for target in target_phases:
        nearest = min(flow.profiles, key=lambda profile: abs(profile.phase_deg - target))
        if nearest not in selected:
            selected.append(nearest)
    return selected


def plot_horizontal_velocity_profiles(solution: SolutionData, flow: FlowFieldData) -> Figure:
    fig, ax = plt.subplots()
    selected = select_profiles(flow, (0, 30, 60, 90, 120, 150, 180))
    for profile in selected:
        ax.plot(
            profile.values[:, 1],
            profile.values[:, 0],
            marker="o",
            label=f"Phase {profile.phase_deg:.0f}°, X/d={profile.x_over_d:.3f}",
        )

    ax.axvline(0.0, linewidth=1.0)
    ax.set_title("Vertical profiles of horizontal velocity", fontweight="bold")
    ax.set_xlabel("u / √(gd)")
    ax.set_ylabel("Elevation above bed, y/d")
    ax.legend(loc="best")
    add_stat_box(
        ax,
        [
            "Phase 0°: crest",
            "Phase 180°: trough",
            "Profiles extend from bed to local free surface",
        ],
        "lower right",
    )
    finish_figure(fig, solution)
    return fig


def plot_vertical_velocity_profiles(solution: SolutionData, flow: FlowFieldData) -> Figure:
    fig, ax = plt.subplots()
    selected = select_profiles(flow, (0, 30, 60, 90, 120, 150, 180))
    for profile in selected:
        ax.plot(
            profile.values[:, 2],
            profile.values[:, 0],
            marker="o",
            label=f"Phase {profile.phase_deg:.0f}°, X/d={profile.x_over_d:.3f}",
        )

    ax.axvline(0.0, linewidth=1.0)
    ax.set_title("Vertical profiles of vertical velocity", fontweight="bold")
    ax.set_xlabel("v / √(gd)")
    ax.set_ylabel("Elevation above bed, y/d")
    ax.legend(loc="best")
    add_stat_box(
        ax,
        [
            "v=0 at crest and trough by symmetry",
            "v=0 at the impermeable bed",
            "Positive v is upward",
        ],
        "lower right",
    )
    finish_figure(fig, solution)
    return fig


def plot_surface_and_bed_kinematics(solution: SolutionData, flow: FlowFieldData) -> Figure:
    fig, ax = plt.subplots()
    phase = np.asarray([profile.phase_deg for profile in flow.profiles], dtype=float)
    surface_u = np.asarray([profile.values[-1, 1] for profile in flow.profiles], dtype=float)
    surface_v = np.asarray([profile.values[-1, 2] for profile in flow.profiles], dtype=float)
    surface_speed = np.hypot(surface_u, surface_v)
    bed_u = np.asarray([profile.values[0, 1] for profile in flow.profiles], dtype=float)

    ax.plot(phase, surface_u, marker="o", label="Surface u / √(gd)")
    ax.plot(phase, surface_v, marker="s", label="Surface v / √(gd)")
    ax.plot(phase, surface_speed, marker="^", label="Surface |V| / √(gd)")
    ax.plot(phase, bed_u, marker="d", label="Bed u / √(gd)")
    ax.axhline(0.0, linewidth=1.0)

    add_stat_box(
        ax,
        [
            f"Maximum surface speed: {np.max(surface_speed):.6g} √(gd)",
            f"Minimum surface speed: {np.min(surface_speed):.6g} √(gd)",
            f"Maximum bed u: {np.max(bed_u):.6g} √(gd)",
            f"Minimum bed u: {np.min(bed_u):.6g} √(gd)",
        ],
        "upper right",
    )

    ax.set_title("Surface and bed kinematics over the half-wave", fontweight="bold")
    ax.set_xlabel("Phase from crest to trough (degrees)")
    ax.set_ylabel("Dimensionless velocity")
    ax.set_xlim(float(np.min(phase)), float(np.max(phase)))
    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=10))
    ax.legend(loc="best")
    finish_figure(fig, solution)
    return fig


# -----------------------------------------------------------------------------
# Main workflow
# -----------------------------------------------------------------------------


def normalize_formats(values: Sequence[str]) -> tuple[str, ...]:
    formats: list[str] = []
    for value in values:
        for token in value.split(","):
            extension = token.strip().lower().lstrip(".")
            if not extension:
                continue
            if extension not in {"png", "svg"}:
                raise ValueError(
                    f"Unsupported output format {extension!r}; use png or svg"
                )
            if extension not in formats:
                formats.append(extension)
    if not formats:
        raise ValueError("At least one output format is required")
    return tuple(formats)


def write_manifest(
    output_dir: Path,
    solution_path: Path,
    solver_path: Path,
    surface_path: Path,
    flowfield_path: Path,
    created: Sequence[Path],
    solution: SolutionData,
    solver: SolverData,
    surface: SurfaceData,
    flow: FlowFieldData,
) -> Path:
    manifest = output_dir / "plot_manifest.txt"
    lines = [
        "STEADY-WAVE RESULT PLOTS",
        "========================",
        "",
        f"Solution input : {solution_path}",
        f"Solver input   : {solver_path}",
        f"Surface input  : {surface_path}",
        f"Flowfield input: {flowfield_path}",
        "",
        case_caption(solution),
        "",
        f"Fourier order: {solver.fourier_order}",
        f"Nonlinear equations: {solver.equation_count}",
        f"Overall convergence: {'PASS' if solver.overall_convergence_passed else 'FAIL'}",
        f"Maximum off-grid defect: {solver.spectral_defect_max:.8e}",
        f"Surface points: {len(surface.x_over_d)}",
        f"Flow profiles: {len(flow.profiles)}",
        f"Vertical points per profile: {flow.profiles[0].values.shape[0]}",
        f"Flow-field nodes: {len(flow.x)}",
        f"Flow-field triangles: {len(flow.triangles)}",
        "",
        "Generated files:",
    ]
    lines.extend(f"  {path.name}" for path in created)
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return manifest


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Read solution.res, solver.res, surface.res and flowfield.res and "
            "create publication-quality steady-wave plots."
        )
    )
    parser.add_argument("--solution", default="solution.res", help="Path to solution.res")
    parser.add_argument("--solver", default="solver.res", help="Path to solver.res")
    parser.add_argument("--surface", default="surface.res", help="Path to surface.res")
    parser.add_argument("--flowfield", default="flowfield.res", help="Path to flowfield.res")
    parser.add_argument(
        "--output-dir",
        default="plots",
        help="Directory for generated plots (default: plots)",
    )
    parser.add_argument(
        "--format",
        dest="formats",
        action="append",
        default=None,
        help="Output format: png or svg. Repeat or use comma-separated values. Default: png,svg",
    )
    parser.add_argument("--dpi", type=int, default=300, help="Raster resolution (default: 300 dpi)")
    parser.add_argument("--show", action="store_true", help="Display figures after writing them")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    try:
        formats = normalize_formats(args.formats or ["png", "svg"])
        if args.dpi < 72:
            raise ValueError("DPI must be at least 72")

        solution_path = resolve_input_path(args.solution, "solution")
        solver_path = resolve_input_path(args.solver, "solver")
        surface_path = resolve_input_path(args.surface, "surface")
        flowfield_path = resolve_input_path(args.flowfield, "flowfield")

        solution = parse_solution(solution_path)
        solver = parse_solver(solver_path)
        solution = merge_solver_metadata(solution, solver)
        surface = parse_surface(surface_path)
        flow = parse_flowfield(flowfield_path)
        validate_result_set(solution, solver, surface, flow)

        output_dir = Path(args.output_dir).expanduser().resolve()
        output_dir.mkdir(parents=True, exist_ok=True)

        configure_matplotlib()

        plot_specs = [
            ("01_solution_summary", lambda: plot_solution_summary(solution)),
            ("02_solver_verification", lambda: plot_solver_summary(solution, solver)),
            ("03_continuation_history", lambda: plot_continuation_history(solution, solver)),
            ("04_final_residual_vector", lambda: plot_final_residual_vector(solution, solver)),
            ("05_wave_profile", lambda: plot_wave_profile(solution, surface)),
            ("06_surface_pressure_check", lambda: plot_surface_pressure_check(solution, surface)),
            ("07_fourier_coefficients", lambda: plot_fourier_coefficients(solution)),
            ("08_depth_scaled_integrals", lambda: plot_depth_scaled_integrals(solution)),
            ("09_velocity_magnitude_vectors", lambda: plot_flow_velocity_magnitude(solution, flow)),
            ("10_horizontal_velocity", lambda: plot_horizontal_velocity(solution, flow)),
            ("11_vertical_velocity", lambda: plot_vertical_velocity(solution, flow)),
            ("12_local_acceleration", lambda: plot_local_acceleration(solution, flow)),
            (
                "13_potential_time_derivative",
                lambda: plot_potential_time_derivative(solution, flow),
            ),
            (
                "14_velocity_gradient_magnitude",
                lambda: plot_velocity_gradient_magnitude(solution, flow),
            ),
            ("15_bernoulli_check", lambda: plot_bernoulli_check(solution, flow)),
            (
                "16_horizontal_velocity_profiles",
                lambda: plot_horizontal_velocity_profiles(solution, flow),
            ),
            (
                "17_vertical_velocity_profiles",
                lambda: plot_vertical_velocity_profiles(solution, flow),
            ),
            ("18_surface_bed_kinematics", lambda: plot_surface_and_bed_kinematics(solution, flow)),
        ]

        created: list[Path] = []
        for stem, factory in plot_specs:
            figure = factory()
            created.extend(
                save_figure(
                    figure, output_dir, stem, formats, args.dpi, close=not args.show
                )
            )

        manifest = write_manifest(
            output_dir,
            solution_path,
            solver_path,
            surface_path,
            flowfield_path,
            created,
            solution,
            solver,
            surface,
            flow,
        )

        print(f"Generated {len(created)} plot files in: {output_dir}")
        print(f"Manifest: {manifest}")
        for path in created:
            print(f"  {path.name}")

        if args.show:
            plt.show()

        return 0

    except (OSError, ValueError) as exc:
        parser.exit(2, f"Error: {exc}\n")


if __name__ == "__main__":
    sys.exit(main())
