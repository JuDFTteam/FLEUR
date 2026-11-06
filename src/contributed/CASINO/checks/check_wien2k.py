#!/usr/bin/env python3
# Copyright (c) 2026 Mike Towler. MIT licence; see LICENSE in src/contributed/CASINO.
"""One lapwfn.data (FLEUR's, say) against another (WIEN2k's, Elk's) for the same crystal: electron counts, the
band energies at Gamma (lowest NB, after one common shift) and the occupied
density at Gamma (the NOCC lowest bands, per cell) at random points. The
two structures may differ by a translation s (A = B + s), found by
matching the atom lists modulo the lattice: the density of A is taken at
x + s where B's is taken at x.
Usage: check_wien2k.py <deck A> <deck B> NOCC NB [npts seed]
"""
import os
import sys
import numpy as np
W = os.environ.get('WIEN2KQMC', os.path.expanduser('~/PROJECTS/wien2kqmc'))
sys.path.insert(0, W)
from wien2kqmc.lapwfn import read_lapwfn
from wien2kqmc.pwfn import read_pwfn
from wien2kqmc.evaluate import Evaluator

de, dw = sys.argv[1], sys.argv[2]
nocc, nb = int(sys.argv[3]), int(sys.argv[4])
npts, seed = (int(sys.argv[5]), int(sys.argv[6])) if len(sys.argv) > 6 else (200, 2)
le, pe = read_lapwfn(f'{de}/lapwfn.data'), read_pwfn(f'{de}/pwfn.data')
lw, pw = read_lapwfn(f'{dw}/lapwfn.data'), read_pwfn(f'{dw}/pwfn.data')
lat = np.asarray(lw['lattice'])

def same(x, y):
    f = np.linalg.solve(lat.T, x - y)
    return np.allclose(f, np.round(f), atol=1e-6)

s = None
for ae in le['atoms']:
    cand = ae['pos'] - lw['atoms'][0]['pos']
    ok = all(any(same(aw['pos'] + cand, a2['pos']) and round(aw['z']) == round(a2['z']) for a2 in le['atoms'])
             for aw in lw['atoms'])
    if ok:
        s = cand; break
if s is None:
    raise SystemExit('no translation maps the atoms of B onto those of A')
print('translation A = B + s, s =', ' '.join(f'{v:.6f}' for v in s), 'bohr')
print(f'electrons: A {le["nelec"]} (valence {le["nvalence"]}), B {lw["nelec"]} (valence {lw["nvalence"]})')
# Gamma: the k point with |k| = 0 in each
ike = next(i for i, k in enumerate(le['kpoints']) if np.linalg.norm(k['k']) < 1e-9)
ikw = next(i for i, k in enumerate(lw['kpoints']) if np.linalg.norm(k['k']) < 1e-9)
Ee = np.array([b['E'] for b in le['kpoints'][ike]['bands']])[:nb]
Ew = np.array([b['E'] for b in lw['kpoints'][ikw]['bands']])[:nb]
shift = (Ee - Ew).mean()
print(f'Gamma band energies, lowest {nb}: mean shift A - B {shift:.6f} Ha; max |dE| after the shift {np.abs(Ee - Ew - shift).max():.3e} Ha')
for i in range(nb):
    print(f'  band {i + 1}: A {Ee[i]:.8f}  B {Ew[i]:.8f}  diff-shift {Ee[i] - Ew[i] - shift:+.3e}')
print(f'Fermi energies: A {le["efermi"]:.6f}, B {lw["efermi"]:.6f}; totals A {le["etot"]:.8f}, B {lw["etot"]:.8f} Ha')
rng = np.random.default_rng(seed)
X = rng.random((npts, 3)) @ lat
eve, evw = Evaluator(le, pe), Evaluator(lw, pw)
ve, _, _, we = eve.evaluate(ike, list(range(nocc)), X + s, nder=0)
vw, _, _, ww = evw.evaluate(ikw, list(range(nocc)), X, nder=0)
rhoe = (np.abs(ve) ** 2).sum(axis=1)
rhow = (np.abs(vw) ** 2).sum(axis=1)
d = np.abs(rhoe - rhow)
ins = we >= 0
print(f'occupied density at Gamma at {npts} random points: max|drho|/max rho {d.max() / rhow.max():.3e} '
      f'(in spheres {d[ins].max() / rhow.max() if ins.any() else 0:.3e} over {ins.sum()} points, '
      f'interstitial {d[~ins].max() / rhow.max() if (~ins).any() else 0:.3e}); regions agree {(ins == (ww >= 0)).sum()}/{npts}')
print(f'  mean rho A {rhoe.mean():.6e}, B {rhow.mean():.6e}; max rho {rhow.max():.6e}')
# the two codes' densities differ by what their linearisation energies make them (the README): only the
# electron counts are a pass or fail here
ok = le['nelec'] == lw['nelec'] and le['nvalence'] == lw['nvalence']
print('RESULT: ' + ('PASS' if ok else 'FAIL electron counts differ'))
sys.exit(0 if ok else 1)
