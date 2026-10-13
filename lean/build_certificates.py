"""Prove exact inequalities on recorded binary64 observations in Lean.

Exported values are represented by their exact rational values. These checks
do not prove how Python produced those values or bound transcendental error.
"""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent


def rational(value):
    p, q = float(value).as_integer_ratio()
    return f'({p} / {q} : ℚ)'


def main():
    path = ROOT/'results/numerical.json'
    report = json.loads(path.read_text())
    lines = ['import Mathlib.Tactic.NormNum', '',
             '/- Exact inequalities on observations; no certified numerical enclosures. -/',
             'namespace NumericalCertificates', 'set_option maxHeartbeats 4000000',
             f'def observationsSHA256 : String := "{hashlib.sha256(path.read_bytes()).hexdigest()}"', '']
    names = []

    def theorem(name, conditions):
        names.append(name)
        lines.extend([f'theorem {name} :', '    ' + ' ∧\n    '.join(conditions) + ' := by norm_num', ''])

    for i, case in enumerate(report['cases']):
        lines.append(f"-- H,T,d,U,current selector: {case['case']}; backend: {case['backend']}")
        theorem(f'case_{i:02d}_recorded_collocation', [
            f"0 ≤ {rational(case['residual_80digits'])}",
            f"{rational(case['residual_80digits'])} ≤ 1/10^10",
            f"{rational(case['condition'])} ≤ 10^11",
            f"{rational(case['point_difference'])} ≤ 1/10^10",
            f"{rational(case['node_difference'])} ≤ 1/10^10"])
        metrics = case['accuracy_metrics']
        if case['converged']:
            assert case['accuracy_verified'] and all(case['checks'].values())
            theorem(f'case_{i:02d}_recorded_accuracy', [
                f"{rational(metrics['surface_residual'])} ≤ 1/10^6",
                *[f'{rational(value)} ≤ 1' for name, value in metrics.items() if name.endswith('_ratio')]])
            theorem(f'case_{i:02d}_recorded_extrema', [
                f"0 ≤ {rational(case['outputs'][quantity])}" for quantity in ('u_bed', 'u_surf', 'acc_max', 'w_max')])
            for quantity in ('Power', 'Sxx', 'KineticEnergy', 'PotentialEnergy', 'Impulse'):
                got, reference = rational(case['outputs'][quantity]), rational(case['reference'][quantity])
                theorem(f'case_{i:02d}_{quantity.lower()}_recorded_budget', [
                    f'|{got} - {reference}| ≤ (1/10^7 : ℚ) + (1/10^8 : ℚ)*|{reference}|'])
        else:
            assert not case['accuracy_verified'] and case['status'] == 'resolution_failure'
            theorem(f'case_{i:02d}_recorded_resolution_rejection', [
                f"(1/10^6 : ℚ) < {rational(metrics['surface_residual'])}",
                f"(1 : ℚ) < {rational(max(v for k,v in metrics.items() if k.endswith('_ratio')))}"])
    for i, record in enumerate(report['jacobian_checks']):
        theorem(f'jacobian_{i}_recorded_difference', [f"{rational(record['max_scaled_error'])} < 1/10^7"])
    lines += ['end NumericalCertificates', '']
    lines += [f'#print axioms NumericalCertificates.{name}' for name in names]
    (ROOT/'NumericalCertificates.lean').write_text('\n'.join(lines) + '\n')
    print(f'Generated {len(names)} recorded-value theorems.')


if __name__ == '__main__': main()
