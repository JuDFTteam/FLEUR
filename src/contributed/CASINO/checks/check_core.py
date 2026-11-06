#!/usr/bin/env python3
# Copyright (c) 2026 Mike Towler. MIT licence; see LICENSE in src/contributed/CASINO.
"""The written core states: the norm over all space of what the file carries (the inner series
from the origin to the first point by its leading terms, the table by a cubic spline in ln r
re-integrated densely, the tabulated tail), the tail's norm beyond the sphere three ways (the
written mesh's trapezoid in ln r, the dense spline, the analytic integral of u(R) exp(-kappa (r-R))
r^2 from R to infinity), the last point's size against the maximum (the floor the writer
tabulates to), and the tail's agreement with the exponential at every written point.
Usage: check_core.py <deck dir>   (WIEN2KQMC names the reader's tree, default ~/PROJECTS/wien2kqmc)
"""
import os
import sys
import numpy as np
from scipy.interpolate import CubicSpline
W = os.environ.get('WIEN2KQMC', os.path.expanduser('~/PROJECTS/wien2kqmc'))
sys.path.insert(0, W)
from wien2kqmc.lapwfn import read_lapwfn

deck = sys.argv[1]
lapw = read_lapwfn(f'{deck}/lapwfn.data')
print('# type state: R | tail norm: trapezoid, dense spline, analytic | inside (dense) | origin piece | total | last/max | max dev of the tail from the exponential')
worst_total, worst_dev = 0.0, 0.0
for t in lapw['types']:
    R = t['rmt']
    Z = t['Z'] if 'Z' in t else t.get('zatom', 0)
    for c in t['core']:
        r = np.asarray(c['r']); u = np.asarray(c['u']); k = c['kappa']
        f = u*u*r*r
        x = np.log(r)
        iR = int(np.argmin(np.abs(r - R)))
        trap_out = np.trapz(f[iR:]*r[iR:], x[iR:])
        sp = CubicSpline(x, f*r)
        xd = np.linspace(x[0], x[-1], 200001)
        dense = sp(xd)
        dense_out = np.trapz(dense[xd >= x[iR]], xd[xd >= x[iR]])
        dense_in = np.trapz(dense[xd <= x[iR]], xd[xd <= x[iR]])
        uR = u[iR]
        ana = uR*uR*(R*R/(2*k) + 2*R/(2*k)**2 + 2/(2*k)**3)
        l = c['l']
        A, c2 = c['A'], c['c2']
        r1 = r[0]
        zz = float(Z)
        orig = A*A*r1**(2*l+3)*(1/(2*l+3) - 2*zz/(l+1)*r1/(2*l+4) + (zz*zz/(l+1)**2 + 2*c2)*r1**2/(2*l+5))
        total = orig + dense_in + dense_out
        tail_model = uR*np.exp(-k*(r[iR:] - R))
        dev = np.abs(u[iR:] - tail_model).max()/max(abs(uR), 1e-300)
        last = abs(u[-1])/np.abs(u).max()
        worst_total = max(worst_total, abs(total - 1))
        worst_dev = max(worst_dev, dev)
        print(f"{t['label']:6s} n={c['n']} l={l}: R {R:.3f} | {trap_out:.4e} {dense_out:.4e} {ana:.4e} | {dense_in:.6f} | {orig:.2e} | {total:.8f} | {last:.1e} | {dev:.1e}")
print(f'# worst |total norm - 1| {worst_total:.1e}; worst tail deviation from the exponential {worst_dev:.1e}')
bad = []
if worst_total > 1e-5:
    bad.append(f'|norm - 1| {worst_total:.1e} > 1e-5')
if worst_dev > 1e-10:
    bad.append(f'tail deviation {worst_dev:.1e} > 1e-10')
print('RESULT: ' + ('PASS' if not bad else 'FAIL ' + '; '.join(bad)))
sys.exit(1 if bad else 0)
