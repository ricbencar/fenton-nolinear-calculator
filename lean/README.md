# Fenton production verification

The solver computes accurately for **19 of the 21 physical configurations in the benchmark matrix**, on both NumPy and Numba. Its acceptance checks identify two steep-wave configurations that require more Fourier resolution. Read `VERIFICATION.md` for the results and the exact meaning of the Lean evidence.

## Files

| File | Purpose |
| --- | --- |
| `inputs/fenton_gui.py` | Production solver and GUI, including independent interval checks. |
| `inputs/README.md` | Mathematical specification. |
| `FentonProduction.lean` | Governing formula, differential and scaling identities. |
| `ProductionAcceptance.lean` | Numerical acceptance and evaluated-maximum contracts. |
| `FiniteSystem.lean` | Enclosure arithmetic, certified π bounds, finite residual and output-budget implications. |
| `NumericalCertificates.lean` | Exact rational checks on numerical observations. |
| `FiniteCertificates.lean` | Kernel-checked integer inequalities for every recorded residual enclosure and nine output budgets per run. |
| `audit.py` | Independent 80-digit equations, kinematics, invariants, Jacobians and interval cross-checks. |
| `check_intervals.py` | Exact-rational interval arithmetic tests and deliberate state/output perturbations. |
| `check_report.py` | Actual output formatter and rejected-result handling. |
| `build_certificates.py`, `build_finite_certificates.py` | Reproduce the recorded-value Lean sources. |
| `sources.json`, `source_map.json`, `MANIFEST.json` | Source identities, reference URLs, function mapping and package hashes. |
| `results/` | Full numerical states, decimal enclosures, execution logs and displayed reports. |

## Reproduce

Use Python 3.12 with Tk bindings installed. NumPy is required by the solver; Numba is optional. The interval checker uses Python's standard `decimal` module and adds no installation dependency. The audit additionally uses mpmath. The supplied requirements reproduce the recorded dependency versions.

```sh
python -m pip install -r requirements.txt
python audit.py
python check_intervals.py
python check_report.py
python inputs/fenton_gui.py --self-test
python build_certificates.py
python build_finite_certificates.py
lake update
lake exe cache get
lake build
```

Lean is installed through Elan. `lean-toolchain` pins Lean 4.34.0; `lakefile.toml` pins Mathlib commit `5ed2965256430c3649e86755f9576b54eca72435`. All tests and all five Lean modules must complete successfully. Passing acceptance tests include deliberate rejection of insufficiently resolved cases.

The audit checks the source hashes before importing the solver. Numerical certificate generation uses exact decimal or binary64 rational values, without rounding endpoints to binary64. The finite certificates use common powers of ten to turn the recorded inequalities into exact integer checks; the scaling implications are proved in `FiniteSystem.lean`.

## Interpretation

Lean checks the mathematical model, enclosure implications and recorded endpoint inequalities. The Python construction of the intervals remains a reviewed and tested implementation, including its Decimal and Taylor calculations. The evidence establishes a stronger finite-system check than residual sampling alone; it does not assert existence or uniqueness of an exact wave, or a formal error bound for the continuous water-wave solution.
