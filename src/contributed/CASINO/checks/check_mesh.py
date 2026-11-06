#!/usr/bin/env python3
# Copyright (c) 2026 Mike Towler. MIT licence; see LICENSE in src/contributed/CASINO.
"""The radial mesh's adequacy for a cubic-spline reader: for every valence
radial function of lapwfn.data, the cubic spline in ln r through the
even-indexed points evaluated at the odd ones (beyond the join radius),
max |error| / max |u|. Usage: check_mesh.py <deck dir>
"""
import os
import sys
import numpy as np
from scipy.interpolate import CubicSpline
W = os.environ.get('WIEN2KQMC', os.path.expanduser('~/PROJECTS/wien2kqmc'))
sys.path.insert(0, W)
from wien2kqmc.lapwfn import read_lapwfn
lapw = read_lapwfn(f'{sys.argv[1]}/lapwfn.data')
print("# type l j kind  n_mesh  ratio  max_interp_err/max|u| (odd points from even, beyond the join; the full mesh has about a sixteenth of it)")
worst = 0.0
for it, t in enumerate(lapw['types'], 1):
    r = t['r']; x = np.log(r)
    for l in range(t['lmax'] + 1):
        for j, fu in enumerate(t['functions'][l], 1):
            u = fu['u']
            i0 = int(np.searchsorted(r, fu['rjoin'])) + 2
            ev = np.arange(i0 % 2, len(r), 2); od = np.arange(1 - i0 % 2, len(r), 2)
            od = od[(od > i0) & (od < ev[-1])]        # no extrapolation beyond the last knot
            s = CubicSpline(x[ev], u[ev])
            err = np.abs(s(x[od]) - u[od]).max() / np.abs(u).max()
            worst = max(worst, err)
            print(f'{it} {l} {j} {fu["kind"]:4s} {len(r)} {r[1]/r[0]:.4f} {err:.3e}')
print(f'# worst {worst:.3e}')
