"""Generate Lean checks of all recorded residual bounds and nine output budgets.

Integer checks use exact common decimal scales. FiniteSystem supplies the
real-number enclosure implications. Interval construction is audited Python.
"""
from pathlib import Path
from fractions import Fraction
from decimal import Decimal
import hashlib
import json

ROOT=Path(__file__).resolve().parent


def number(value, kind='ℚ'):
    f=Fraction(value)
    return f'({f.numerator}/{f.denominator} : {kind})'


def main():
    path=ROOT/'results/numerical.json'
    report=json.loads(path.read_text())
    lines=['import FiniteSystem','import Mathlib.Tactic.NormNum','',
           'namespace FiniteCertificates','set_option maxHeartbeats 0',
           'set_option maxRecDepth 10000',
           f'def observationsSHA256 : String := "{hashlib.sha256(path.read_bytes()).hexdigest()}"','']
    names=[]
    for i,case in enumerate(report['cases']):
        c=case['finite_system_certificate']
        assert c['passed']
        assert len(c['residual_intervals'])==2*c['order']+10
        bound_name=f'case_{i:02d}_bounds'
        digits=max(10,max(-Decimal(v).as_tuple().exponent for pair in c['residual_intervals'] for v in pair))
        denominator=10**digits
        integer_bounds=[tuple(int(Fraction(v)*denominator) for v in pair) for pair in c['residual_intervals']]
        assert all(Fraction(v)*denominator==integer for pair,ints in zip(c['residual_intervals'],integer_bounds) for v,integer in zip(pair,ints))
        threshold=Fraction(c['residual_tolerance'])*denominator
        assert threshold.denominator==1
        lines += [f'-- Every integer endpoint represents its value divided by {denominator}.',
                  f'def case_{i:02d}_denominator : ℕ := {denominator}',
                  f'def {bound_name} : List (ℤ × ℤ) := [',
                  ',\n'.join(f'  ({lo}, {hi})' for lo,hi in integer_bounds),']','']
        name=f'case_{i:02d}_all_residual_bounds';names.append(name)
        tolerance=f'({threshold.numerator} : ℤ)'
        lines += [f'theorem {name} : (0 : ℤ) < {denominator} ∧',
                  f'    (10^10 : ℤ)*{tolerance} ≤ {denominator} ∧',
                  f'    ∀ p ∈ {bound_name}, -{tolerance} ≤ p.1 ∧ p.1 ≤ p.2 ∧ p.2 ≤ {tolerance} := by',
                  '  decide','']
        for key,record in c['outputs'].items():
            name=f'case_{i:02d}_{key.lower()}_output_bound';names.append(name)
            values=record['interval']+[record[k] for k in ('reported','error_upper','budget_lower')]
            digits=max(0,max(-Decimal(v).as_tuple().exponent for v in values))
            scale=10**digits
            ints=[int(Fraction(v)*scale) for v in values]
            assert all(Fraction(v)*scale==a for v,a in zip(values,ints))
            lo,hi,value,error,budget=(f'({a} : ℤ)' for a in ints)
            minimum=0 if ints[0]<=0<=ints[1] else min(abs(ints[0]),abs(ints[1]))
            lines += [f'-- Common denominator: {scale}.',
                      f'theorem {name} :',f'    {lo} ≤ {hi} ∧',
                      f'    {value}-{error} ≤ {lo} ∧ {hi} ≤ {value}+{error} ∧',
                      f'    0 ≤ {error} ∧ {error} ≤ {budget} ∧ (0 : ℤ) < {scale} ∧',
                      f'    (10^8 : ℤ)*{budget} ≤ 10*{scale}+{minimum} := by decide','']
    lines += ['end FiniteCertificates','']+[f'#print axioms FiniteCertificates.{n}' for n in names]
    (ROOT/'FiniteCertificates.lean').write_text('\n'.join(lines)+'\n')
    print('Generated',len(names),'finite-system enclosure theorems.')


if __name__=='__main__':main()
