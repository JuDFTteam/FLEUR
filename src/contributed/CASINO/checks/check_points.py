#!/usr/bin/env python3
# Copyright (c) 2026 Mike Towler. MIT licence; see LICENSE in src/contributed/CASINO.
"""The file against FLEUR's own orbitals at points: wien2kqmc's reference evaluator on
lapwfn.data and pwfn.data against lapwfn_points.out (the writer run with lapwfn_points.in),
whose first seven columns are ik ib ipt atom r re im, FLEUR's local-frame value in its own
frame. Reports, per k point, the largest |psi_file - psi_fleur| relative to max|psi_fleur|
over the points, separately for points on mesh radii (FLEUR's value exact there), random
sphere points (two interpolations), surface points and interstitial points, over all bands;
and beside it the same with each band measured against its own maximum (the stricter number,
the one PREDICTIONS registered).
Usage: check_points.py <deck dir> <lapwfn_points.out> <lapwfn_points.in> <NM> <NS> <NSF> <nat>
"""
import os
import sys
import numpy as np
W = os.environ.get('WIEN2KQMC', os.path.expanduser('~/PROJECTS/wien2kqmc'))
sys.path.insert(0, W)
from wien2kqmc.lapwfn import read_lapwfn
from wien2kqmc.pwfn import read_pwfn
from wien2kqmc.evaluate import Evaluator

deck, fpts, fin = sys.argv[1], sys.argv[2], sys.argv[3]
NM, NS, NSF, nat = (int(x) for x in sys.argv[4:8])
lapw = read_lapwfn(f'{deck}/lapwfn.data')
pw = read_pwfn(f'{deck}/pwfn.data')
ev = Evaluator(lapw, pw)
P = np.loadtxt(fin, skiprows=1)
npts = len(P)
kind = np.array(sum([['mesh'] * NM + ['sphere'] * NS + ['surface'] * NSF for _ in range(nat)], []) + ['inter'] * (npts - nat * (NM + NS + NSF)))
E = {}
for line in open(fpts):
    if line.startswith('#'):
        continue
    t = line.split()
    E[(int(t[0]), int(t[1]), int(t[2]))] = (int(t[3]), float(t[4]), complex(float(t[5]), float(t[6])))
print('# ik  kind      npts  max|dpsi|/max|psi|  per band    regions')
worst = {}
worstb = {}
agree_all = True
for ik, kp in enumerate(lapw['kpoints']):
    nb = len(kp['bands'])          # per k point: the list may end at a different level boundary at each
    val, _, _, where = ev.evaluate(ik, list(range(nb)), P, nder=0)
    for kd in ('mesh', 'sphere', 'surface', 'inter'):
        sel = np.where(kind == kd)[0]
        if len(sel) == 0:
            continue
        dmax, ref, dband = 0.0, 0.0, 0.0
        for ib in range(nb):
            ve = np.array([E[(ik + 1, ib + 1, i + 1)][2] for i in sel])
            vf = val[sel, ib]
            ref = max(ref, np.abs(ve).max())
            dd = np.abs(vf - ve).max()
            dmax = max(dmax, dd)
            dband = max(dband, dd / max(np.abs(ve).max(), 1e-300))
        reg = [E[(ik + 1, 1, i + 1)][0] for i in sel]
        regf = [where[i] + 1 if where[i] >= 0 else 0 for i in sel]
        agree = sum(a == b for a, b in zip(reg, regf))
        agree_all = agree_all and agree == len(sel)
        print(f'{ik + 1:3d}  {kd:7s} {len(sel):6d}  {dmax / ref:.3e}   {dband:.3e}   regions agree {agree}/{len(sel)}')
        worst[kd] = max(worst.get(kd, 0.0), dmax / ref)
        worstb[kd] = max(worstb.get(kd, 0.0), dband)
    # the real flag: a band flagged real must evaluate real (the evaluator keeps every band complex)
    flags = [b['real'] for b in kp['bands']]
    imx = 0.0
    for ib in range(nb):
        if flags[ib]:
            imx = max(imx, np.abs(val[:, ib].imag).max()/max(np.abs(val[:, ib]).max(), 1e-300))
    print(f'{ik + 1:3d}  real-flagged bands {sum(flags)}/{nb}: max |Im psi|/max|psi| over the points {imx:.3e}')
    worst['real-flag'] = max(worst.get('real-flag', 0.0), imx)
print('# worst over k: ' + ', '.join(f'{k} {v:.3e}' for k, v in worst.items()))
print('# worst over k, each band by its own maximum: ' + ', '.join(f'{k} {v:.3e}' for k, v in worstb.items()))
# the bounds a sound file meets by a wide margin (the registered figures are an order tighter): the
# script fails on a broken file, not on a case's own truncation
BOUND = {'mesh': 1e-10, 'sphere': 5e-6, 'surface': 1e-8, 'inter': 1e-10}
bad = [f'{k} {worstb[k]:.2e} > {b:.0e}' for k, b in BOUND.items() if k in worstb and worstb[k] > b]
if worst.get('real-flag', 0.0) > 1e-10:
    bad.append(f'real-flag {worst["real-flag"]:.2e} > 1e-10')
if not agree_all:
    bad.append('regions disagree')
print('RESULT: ' + ('PASS' if not bad else 'FAIL ' + '; '.join(bad)))
sys.exit(1 if bad else 0)
