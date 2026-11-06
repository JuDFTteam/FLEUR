#!/usr/bin/env python3
# Copyright (c) 2026 Mike Towler. MIT licence; see LICENSE in src/contributed/CASINO.
"""check_levels.py <deck dir> <norealise deck dir> <lapwfn_points.in> <NOCC>
The real combinations within degenerate levels at k = G/2, tested against a deck written from the same
stored density with the combinations switched off (output/lapwfn/@realise="F"), which carries FLEUR's own
vectors. Per k point and level (consecutive eigenvalues within 1e-6 Ha in the reference deck): the
eigenvalue labels' largest difference (the deck's are the expectation energies of its bands), the
difference of the level projectors V V^+ on the full coefficient vector (spheres and plane waves),
which no choice of basis within the level changes, the flags, and for each deck the closure residual
under the conjugation map of a real function (c_G to conj c_{-G-2k}, c_lm to (-1)^m conj c_{l,-m}) and
the singular values of the real span (n of order one and the rest at rounding when a real basis
exists). Then the occupied density (the NOCC lowest bands at every k point) at the points of
lapwfn_points.in from both decks by the reference evaluator. Fails when a projector differs by more
than 1e-10, a label by more than 1e-6 Ha, the density by more than 1e-10 of its maximum, or a level at
k = G/2 is not written real. The conjugation map and the real span follow an independent check's closure.py."""
import os
import sys
import numpy as np
W = os.environ.get('WIEN2KQMC', os.path.expanduser('~/PROJECTS/wien2kqmc'))
sys.path.insert(0, W)
from wien2kqmc.lapwfn import read_lapwfn
from wien2kqmc.pwfn import read_pwfn
from wien2kqmc.evaluate import Evaluator

DA, DB, fin, nocc = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
LA, PA = read_lapwfn(f'{DA}/lapwfn.data'), read_pwfn(f'{DA}/pwfn.data')
LB, PB = read_lapwfn(f'{DB}/lapwfn.data'), read_pwfn(f'{DB}/pwfn.data')
lat = np.asarray(LA['lattice'])
bad = []


def vectors(L, P, ik):
    kp, pk = L['kpoints'][ik], P['kpoints'][ik]
    V = np.array([np.concatenate([np.concatenate([c.ravel() for c in b['coef']]), pb['c']])
                  for b, pb in zip(kp['bands'], pk['bands'])]).T
    return V, np.array([b['E'] for b in kp['bands']]), ''.join('T' if b['real'] else 'F' for b in kp['bands'])


def conjmap_factory(L, P, ik):
    """the map K on the full coefficient vector at k point ik (None when k is not at G/2)"""
    kp = L['kpoints'][ik]
    kf = kp['k'] @ lat.T / (2 * np.pi)
    tk = np.rint(2 * kf)
    if np.abs(2 * kf - tk).max() > 1e-8:
        return None
    G = np.asarray(P['gvecs'])
    Gi = np.rint(G @ lat.T / (2 * np.pi)).astype(int)
    index = {tuple(g): i for i, g in enumerate(Gi)}
    partner = np.array([index.get(tuple(-g - tk.astype(int)), -1) for g in Gi])
    ty = [L['types'][a['type'] - 1] for a in L['atoms']]
    perm, sign, off = [], [], 0
    for ia, t in enumerate(ty):
        nlm = (t['lmax'] + 1) ** 2
        nf = kp['bands'][0]['coef'][ia].shape[1]
        for l in range(t['lmax'] + 1):
            for m in range(-l, l + 1):
                for j in range(nf):
                    perm.append(off + (l * l + l - m) * nf + j)
                    sign.append((-1) ** m)
        off += nlm * nf
    perm, sign, ns = np.array(perm), np.array(sign), off

    def K(v):
        s, p = v[:ns], v[ns:]
        return np.concatenate([sign * np.conj(s[perm]), np.where(partner >= 0, np.conj(p[np.maximum(partner, 0)]), 0)])
    return K


def span_numbers(V, K):
    """closure residual of span(V) under K, and the singular values (relative) n and n+1 of the real span"""
    n = V.shape[1]
    Q, _ = np.linalg.qr(V)
    CV = np.array([K(V[:, q]) for q in range(n)]).T
    res = np.abs(CV - Q @ (Q.conj().T @ CV)).max() / np.abs(V).max()
    A = np.concatenate([V + CV, 1j * (V - CV)], axis=1)
    S = np.linalg.svd(np.concatenate([A.real, A.imag], axis=0), compute_uv=False)
    return res, S[n - 1] / S[0], (S[n] / S[0] if len(S) > n else 0.0)


print('# ik bands(level)  flags A/B  |dE label|  |dP|/max|P|  closure A  closure B  s_n A  s_n+1 A  (A the deck, B realise="F")')
for ik in range(len(LA['kpoints'])):
    VA, EA, fA = vectors(LA, PA, ik)
    VB, EB, fB = vectors(LB, PB, ik)
    if VA.shape != VB.shape:
        bad.append(f'k {ik + 1}: the decks hold {VA.shape[1]} and {VB.shape[1]} bands')
        print(f'{ik + 1:3d} band counts differ: {VA.shape[1]} against {VB.shape[1]}')
        continue
    K = conjmap_factory(LA, PA, ik)
    nb = len(EB)
    i = 0
    while i < nb:
        j = i
        while j + 1 < nb and EB[j + 1] - EB[j] < 1e-6:
            j += 1
        sel = slice(i, j + 1)
        dE = np.abs(EA[sel] - EB[sel]).max()
        Pa, Pb = VA[:, sel] @ VA[:, sel].conj().T, VB[:, sel] @ VB[:, sel].conj().T
        dP = np.abs(Pa - Pb).max() / np.abs(Pb).max()
        if K is not None:
            ra, sna, sn1a = span_numbers(VA[:, sel], K)
            rb, _, _ = span_numbers(VB[:, sel], K)
            kstr = f'{ra:9.1e}  {rb:9.1e}  {sna:7.1e} {sn1a:8.1e}'
            if 'F' in fA[sel]:
                bad.append(f'k {ik + 1} bands {i + 1}-{j + 1}: not written real at k = G/2 (flags {fA[sel]})')
        else:
            kstr = '   (not at G/2)'
        print(f'{ik + 1:3d} {i + 1:3d}-{j + 1:<3d} ({j - i + 1})   {fA[sel]}/{fB[sel]}   {dE:9.1e}  {dP:9.1e}  {kstr}')
        if dP > 1e-10:
            bad.append(f'k {ik + 1} bands {i + 1}-{j + 1}: projectors differ by {dP:.1e}')
        if dE > 1e-6:
            bad.append(f'k {ik + 1} bands {i + 1}-{j + 1}: labels differ by {dE:.1e} Ha')
        i = j + 1
# the occupied density at the points
X = np.loadtxt(fin, skiprows=1)
eva, evb = Evaluator(LA, PA), Evaluator(LB, PB)
rhoa = np.zeros(len(X)); rhob = np.zeros(len(X))
for ik in range(len(LA['kpoints'])):
    va, _, _, _ = eva.evaluate(ik, list(range(nocc)), X, nder=0)
    vb, _, _, _ = evb.evaluate(ik, list(range(nocc)), X, nder=0)
    rhoa += (np.abs(va) ** 2).sum(axis=1)
    rhob += (np.abs(vb) ** 2).sum(axis=1)
drho = np.abs(rhoa - rhob).max() / rhob.max()
print(f'# occupied density ({nocc} bands at every k point) at {len(X)} points: max|drho|/max rho {drho:.2e}')
if drho > 1e-10:
    bad.append(f'the occupied density differs by {drho:.1e}')
print('RESULT: ' + ('PASS' if not bad else 'FAIL ' + '; '.join(bad)))
sys.exit(1 if bad else 0)
