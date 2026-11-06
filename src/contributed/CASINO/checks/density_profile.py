#!/usr/bin/env python3
# Copyright (c) 2026 Mike Towler. MIT licence; see LICENSE in src/contributed/CASINO.
"""The valence (NOCC lowest bands at Gamma) density of two lapwfn.data decks along radial lines
from the first atom, and the difference profile: where the two codes' densities differ.
Usage: density_profile.py <deck A> <deck B> NOCC"""
import os, sys
import numpy as np
W = os.environ.get('WIEN2KQMC', os.path.expanduser('~/PROJECTS/wien2kqmc'))
sys.path.insert(0, W)
from wien2kqmc.lapwfn import read_lapwfn
from wien2kqmc.pwfn import read_pwfn
from wien2kqmc.evaluate import Evaluator
da, db, nocc = sys.argv[1], sys.argv[2], int(sys.argv[3])
la, pa = read_lapwfn(f'{da}/lapwfn.data'), read_pwfn(f'{da}/pwfn.data')
lb, pb = read_lapwfn(f'{db}/lapwfn.data'), read_pwfn(f'{db}/pwfn.data')
lat = np.asarray(la['lattice'])
# the translation A = B + s from the atom lists
def same(x, y):
    f = np.linalg.solve(lat.T, x - y); return np.allclose(f, np.round(f), atol=1e-6)
s = None
for atb in lb['atoms']:
    cand = np.asarray(la['atoms'][0]['pos']) - np.asarray(atb['pos'])
    if all(any(same(np.asarray(x['pos']), np.asarray(y['pos']) + cand) for y in lb['atoms']) for x in la['atoms']):
        s = cand; break
ea, eb = Evaluator(la, pa), Evaluator(lb, pb)
pos = np.asarray(la['atoms'][0]['pos'])
def rho(ev, X):
    v, _, _, _ = ev.evaluate(0, list(range(nocc)), X, nder=0)
    return 2 * (np.abs(v) ** 2).sum(axis=1)
for name, d in (('111', np.array([1, 1, 1]) / 3 ** 0.5), ('100', np.array([1.0, 0, 0])), ('110', np.array([1, 1, 0]) / 2 ** 0.5)):
    rs = np.array([0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.1, 1.25, 1.3, 1.35, 1.5, 1.8, 2.2, 2.6, 2.9])
    X = pos + rs[:, None] * d[None, :]
    ra, rb = rho(ea, X), rho(eb, X + s)
    print(f'# direction {name}: r, rho_A, rho_B, (A-B)/B, A-B')
    for r, a, b in zip(rs, ra, rb):
        print(f'{r:6.3f} {a:.6e} {b:.6e} {(a - b) / b:+.3e} {a - b:+.3e}')
