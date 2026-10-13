"""Exercise report rendering and numerical rejection without a Tk event loop."""
from pathlib import Path
import importlib.util
import os
import sys

os.environ['OPENBLAS_NUM_THREADS'] = '1'
ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('fenton_report_audit', ROOT/'inputs/fenton_gui.py')
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
module.NUMBA_AVAILABLE = False


class Widget:
    def config(self, **kwargs): pass
    def delete(self, *args): pass
    def yview_moveto(self, *args): pass


class App:
    txt_output = Widget()
    btn_calc = Widget()
    def __init__(self): self.lines = []; self.errors = []
    def log(self, text, file=None): self.lines.append(text)
    def _finalize_error(self, *args): self.errors.append(args)


def main():
    os.chdir(ROOT/'results')
    zero = module.FentonStreamFunction(1, 9, 10, 0)
    current = module.FentonStreamFunction(1, 9, 10, -.5)
    zero.solve(); current.solve()
    assert zero.converged and current.converged
    assert current.u_bed_crest < 0 < current.u_bed
    app = App()
    module.FentonApp._finalize_output(app, zero, current, 1, 9, 10, -.5)
    assert not app.errors
    text = ''.join(app.lines)
    (ROOT/'results/rendered_report.txt').write_text(text)
    row = next(line for line in text.splitlines() if 'Bed max |u| (estimated)' in line)
    assert '0.88334' in row and '-0.08487' not in row
    print('Opposing-current magnitude report: PASS')

    zero = module.FentonStreamFunction(2, 3, 5, 0)
    current = module.FentonStreamFunction(2, 3, 5, 1)
    zero.solve(); current.solve()
    assert not zero.converged and current.converged
    app = App()
    module.FentonApp._finalize_output(app, zero, current, 2, 3, 5, 1)
    assert not app.errors
    text = ''.join(app.lines)
    (ROOT/'results/partial_report.txt').write_text(text)
    row = next(line for line in text.splitlines() if 'Bed max |u| (estimated)' in line)
    assert row.split('|')[-4].strip() == '-'
    print('Rejected no-current column suppressed: PASS')

    app = App()
    module.FentonApp._finalize_output(app, zero, zero, 2, 3, 5, 0)
    assert app.errors and 'accuracy budgets' in str(app.errors)
    print('Entirely rejected result reports failure: PASS')

    current.Power *= 1.001
    current._verify_accuracy()
    assert not current.converged and not current.accuracy_verified
    assert current.solution_status == 'resolution_failure'
    assert current.accuracy_metrics['Power_integral_ratio'] > 1
    print('Perturbed energy flux rejected: PASS')


if __name__ == '__main__': main()
