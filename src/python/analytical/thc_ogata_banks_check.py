"""Compare a THC ogata_banks.in run against the Ogata-Banks (1961) solution.

Reuses AnalyticalSolution.ogata_banks from src/python/analytical_solutions.py,
whose hydrodynamic dispersion convention matches THC exactly:
  D = dispersivity*U + porosity*tortuosity*saturation*D_mol,  then D/porosity
and pore velocity v = U/porosity in the solution.

Usage:
  python3 thc_ogata_banks_check.py <tecplot_point_file> <time_in_days>
Defaults match regression_tests/thc/ogata_banks.in.

Author: Piyoosh Jaysaval, PNNL (piyoosh.jaysaval@pnnl.gov)
Date: August 20, 2026

"""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from analytical_solutions import AnalyticalSolution

# deck parameters (regression_tests/thc/ogata_banks.in)
C0 = 1.e-3            # inlet concentration [M]
U = 1.0 / 86400.      # Darcy flux [m/s]
POROSITY = 0.25
TORTUOSITY = 1.0
SATURATION = 1.0
DISPERSIVITY = 0.5    # [m]
D_MOL = 1.e-9         # [m^2/s]


def read_tecplot_point(path):
    x, c = [], []
    icol = None
    with open(path) as f:
        for line in f:
            if line.lstrip().startswith('VARIABLES'):
                names = [s.strip().strip('"') for s in
                         line.split('=', 1)[1].split(',')]
                icol = next(i for i, n in enumerate(names)
                            if 'Solute' in n)
            elif line.lstrip().startswith('ZONE') or icol is None:
                continue
            else:
                vals = line.split()
                if len(vals) > icol:
                    x.append(float(vals[0]))
                    c.append(float(vals[icol]))
    return x, c


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'ogata_banks-001.tec'
    t_days = float(sys.argv[2]) if len(sys.argv) > 2 else 1.25
    t = t_days * 86400.

    soln = AnalyticalSolution(0.0, C0, U, SATURATION, D_MOL, DISPERSIVITY,
                              TORTUOSITY, POROSITY, retardation=1.0,
                              half_life=1.e30) \
        if 'retardation' in AnalyticalSolution.__init__.__code__.co_varnames \
        else AnalyticalSolution(0.0, C0, U, SATURATION, D_MOL, DISPERSIVITY,
                                TORTUOSITY, POROSITY, 1.0, 1.e30)

    x, c_num = read_tecplot_point(path)
    max_abs = max_rel = l2 = 0.
    for xi, ci in zip(x, c_num):
        ca = soln.ogata_banks(xi, t)
        max_abs = max(max_abs, abs(ci - ca))
        max_rel = max(max_rel, abs(ci - ca) / C0)
        l2 += (ci - ca) ** 2
    l2 = math.sqrt(l2 / len(x)) / C0
    print(f'cells: {len(x)}  t = {t_days} d')
    print(f'max |C_num - C_ana|      : {max_abs:.6e} M')
    print(f'max |C_num - C_ana|/C0   : {max_rel:.6e}')
    print(f'rms error / C0           : {l2:.6e}')
    return 0 if max_rel < 0.05 else 1


if __name__ == '__main__':
    raise SystemExit(main())
