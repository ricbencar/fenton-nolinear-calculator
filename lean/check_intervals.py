"""Exercise interval containment independently and reject perturbed states."""
from pathlib import Path
import importlib.util
import os
import sys
from fractions import Fraction
import random

os.environ['OPENBLAS_NUM_THREADS']='1'
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('fenton_intervals',ROOT/'inputs/fenton_gui.py')
mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
mod.NUMBA_AVAILABLE=False
I=mod._Interval


def contains(interval, value):
    return Fraction(interval.lo)<=value<=Fraction(interval.hi)


def main():
    rng=random.Random(104729)
    for _ in range(250):
        a,b=sorted(rng.uniform(-100,100) for _ in range(2))
        c,d=sorted(rng.uniform(-100,100) for _ in range(2))
        x,y=I(a,b),I(c,d)
        for p in (Fraction(a),Fraction(b),(Fraction(a)+Fraction(b))/2):
            for q in (Fraction(c),Fraction(d),(Fraction(c)+Fraction(d))/2):
                assert contains(x+y,p+q) and contains(x-y,p-q)
                assert contains(x*y,p*q) and contains(x.square(),p*p)
                if not c<=0<=d:assert contains(x/y,p/q)
    try:
        I(1)/I(-1,1)
        raise AssertionError('Zero-containing divisor accepted')
    except FloatingPointError:pass
    print('Signed interval arithmetic: 250 exact-rational sample families passed.')
    s=mod.FentonStreamFunction(1,9,10,-.5);s.solve()
    assert s.converged and s.finite_system_verified
    original=s.z[10];s.z[10]+=1e-6
    assert not s._verify_finite_system()
    s.z[10]=original
    assert s._verify_finite_system()
    s.Power*=1.001
    assert not s._verify_finite_system()
    assert mod.Decimal(s.finite_system_certificate['outputs']['Power']['error_upper'])>mod.Decimal(s.finite_system_certificate['outputs']['Power']['budget_lower'])
    print('Finite-state and output perturbations rejected; restored state accepted.')


if __name__=='__main__':main()
