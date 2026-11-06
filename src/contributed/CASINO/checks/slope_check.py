#!/usr/bin/env python3
# Copyright (c) 2026 Mike Towler. MIT licence; see LICENSE in src/contributed/CASINO.
"""slope_check.py <deck dir>
The sphere surface channel by channel, with the derivative rule stated and varied (an independent check's
lesson: a continuity figure is checked with an accurate derivative as well as the code's own rule).
On the plane-wave side each (l, m) component of the sum over G of c_G exp(i(k+G).x) on the sphere of
radius R about an atom at tau is, by the Rayleigh expansion, 4 pi sum_G c_G exp(i(k+G).tau) i^l
j_l(qR) Y*_lm(q^) with q = k+G, and its radial derivative the same with q j_l'(qR); on the sphere
side they are sum_j c_lm,j u_j(R) and sum_j c_lm,j u_j'(R), with u_j'(R) from a polynomial in ln r
through the last 4, 6 and 8 mesh points (degree one less) and a degree-8 fit through the last 12
(wien2kqmc's rule, the one check_surface.py uses). Per k point and atom the largest channel mismatch
of the value and of the slope relative to the largest channel value (slope) over l <= L, m and the
bands, for channels holding at least 5 per cent of the surface weight. The pointwise difference on
the surface, which also holds the plane waves' content beyond L, is check_surface.py's. Fails when
the value mismatch exceeds 1e-9 or the 12-point slope mismatch 1e-5. Adapted from excitingqmc's
slope_check.py. Copyright (c) 2026 Mike Towler, MIT licence."""
import os
import sys
import numpy as np
from scipy.special import spherical_jn
W = os.environ.get('WIEN2KQMC', os.path.expanduser('~/PROJECTS/wien2kqmc'))
sys.path.insert(0, W)
from wien2kqmc.lapwfn import read_lapwfn
from wien2kqmc.pwfn import read_pwfn
from wien2kqmc.harmonics import Harmonics, ylm_wien

deck = sys.argv[1]
lapw = read_lapwfn(f'{deck}/lapwfn.data')
pw = read_pwfn(f'{deck}/pwfn.data')
G = np.asarray(pw['gvecs'], float)
lmax = max(t['lmax'] for t in lapw['types'])
harm = Harmonics(lmax)
RULES = (('4pt', 4, 3), ('6pt', 6, 5), ('8pt', 8, 7), ('12pt', 12, 8))


def end_deriv(r, u, npts, deg):
    """value and d/dr at r[-1] of the polynomial of the given degree in ln r through the last npts points"""
    x = np.log(r[-npts:] / r[-1])
    c = np.polyfit(x, u[-npts:], deg)
    p = np.poly1d(c)
    return p(0.0), p.deriv()(0.0) / r[-1]


print(f'# {deck}: the sphere surface channel by channel (l <= L, channels of 5 per cent or more of the surface weight); '
      'the slope on the sphere side by a polynomial in ln r through the last 4, 6 and 8 mesh points (degree one less) '
      'and a degree-8 fit through the last 12')
print('# ik  atom  rule   max|dv_lm|/max|v_lm|   max|dslope_lm|/max|slope_lm|   (over bands, l, m)')
worst = {}
for ik, kp in enumerate(lapw['kpoints']):
    k = kp['k']
    nb = len(kp['bands'])
    q = G + k[None, :]
    qn = np.linalg.norm(q, axis=1)
    Yq = ylm_wien(q, lmax, harm).conj()            # (nlm, nG), Y*_lm(q^)
    for ia, at in enumerate(lapw['atoms']):
        t = lapw['types'][at['type'] - 1]
        r = t['r']; R = t['rmt']
        L = t['lmax']; nlm = (L + 1) ** 2
        ph = np.exp(1j * (q @ at['pos']))         # (nG,)
        Av = np.zeros((nlm, len(G)), complex); As = np.zeros((nlm, len(G)), complex)
        for l in range(L + 1):
            jl = spherical_jn(l, qn * R); djl = spherical_jn(l, qn * R, derivative=True)
            f = 4 * np.pi * (1j ** l) * ph
            for m in range(-l, l + 1):
                i = l * (l + 1) + m
                Av[i] = f * jl * Yq[i]
                As[i] = f * qn * djl * Yq[i]
        ends = {rule: {(l, j): end_deriv(r, fu['u'], npts, deg)
                       for l in range(L + 1) for j, fu in enumerate(t['functions'][l])}
                for rule, npts, deg in RULES}
        vpw = np.zeros((nlm, nb), complex); spw = np.zeros((nlm, nb), complex)
        vin = {rule: np.zeros((nlm, nb), complex) for rule in ends}
        sin = {rule: np.zeros((nlm, nb), complex) for rule in ends}
        for ib in range(nb):
            c = np.asarray(pw['kpoints'][ik]['bands'][ib]['c'])
            vpw[:, ib] = Av @ c
            spw[:, ib] = As @ c
            coef = kp['bands'][ib]['coef'][ia]        # (nlm, nmax)
            for rule, uR in ends.items():
                for l in range(L + 1):
                    sl = slice(l * l, (l + 1) ** 2)
                    for j in range(len(t['functions'][l])):
                        u0, u1 = uR[(l, j)]
                        vin[rule][sl, ib] += coef[sl, j] * u0
                        sin[rule][sl, ib] += coef[sl, j] * u1
        # the channels of weight: (l, band) pairs holding 5 per cent or more of the band's surface weight
        wl = np.zeros((L + 1, nb))
        for l in range(L + 1):
            wl[l] = (np.abs(vpw[l * l:(l + 1) ** 2, :]) ** 2).sum(axis=0)
        wl /= np.maximum(wl.sum(axis=0), 1e-300)
        mask = np.zeros((nlm, nb), bool)
        for l in range(L + 1):
            mask[l * l:(l + 1) ** 2, :] = wl[l][None, :] >= 0.05
        for rule in ends:
            dv = np.abs((vin[rule] - vpw)[mask]).max() / np.abs(vpw).max()
            ds = np.abs((sin[rule] - spw)[mask]).max() / np.abs(spw).max()
            print(f'{ik + 1:3d}  {ia + 1:4d}  {rule:4s}   {dv:.3e}              {ds:.3e}')
            w = worst.get(rule, (0.0, 0.0))
            worst[rule] = (max(w[0], dv), max(w[1], ds))
for rule, (dv, ds) in worst.items():
    print(f'# worst over k and atoms, {rule}: value {dv:.3e}, slope {ds:.3e}')
bad = []
if max(v for v, s in worst.values()) > 1e-9:
    bad.append(f'value mismatch {max(v for v, s in worst.values()):.1e} > 1e-9')
if worst['12pt'][1] > 1e-5:
    bad.append(f'12-point slope mismatch {worst["12pt"][1]:.1e} > 1e-5')
print('RESULT: ' + ('PASS' if not bad else 'FAIL ' + '; '.join(bad)))
sys.exit(1 if bad else 0)
