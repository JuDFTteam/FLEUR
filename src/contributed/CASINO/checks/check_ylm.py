#!/usr/bin/env python3
# Copyright (c) 2026 Mike Towler. MIT licence; see LICENSE in src/contributed/CASINO.
"""FLEUR's harmonics against the reference reader's: lapwfn_ylm.out (the writer writes
ylm4 at fifty directions, l <= 10) against wien2kqmc's Harmonics, the Y_lm the file's
coefficients are meant to multiply. Usage: check_ylm.py <lapwfn_ylm.out>
"""
import os
import sys
import numpy as np
W = os.environ.get('WIEN2KQMC', os.path.expanduser('~/PROJECTS/wien2kqmc'))
sys.path.insert(0, W)
from wien2kqmc.harmonics import Harmonics

lmax = 10
H = Harmonics(lmax)
worst, where, n = 0.0, None, 0
for line in open(sys.argv[1]):
    if line.startswith('#'):
        continue
    t = line.split()
    l, m = int(t[0]), int(t[1])
    v = np.array([float(x) for x in t[2:5]])
    yf = complex(float(t[5]), float(t[6]))
    Y = H.ylm(v[None, :])[:, 0]
    lm = l * (l + 1) + m
    d = abs(yf - Y[lm])
    n += 1
    if d > worst:
        worst, where = d, (l, m)
print(f'max |Y_fleur - Y_ref| over {n} values (fifty directions, l <= {lmax}): {worst:.3e} at (l, m) = {where}')
print('RESULT: ' + ('PASS' if worst < 1e-10 else f'FAIL harmonics differ by {worst:.2e} (bound 1e-10)'))
sys.exit(0 if worst < 1e-10 else 1)
