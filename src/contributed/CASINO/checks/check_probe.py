#!/usr/bin/env python3
# Copyright (c) 2026 Mike Towler. MIT licence; see LICENSE in src/contributed/CASINO.
"""The frame test, from lapwfn_points.out alone: at the points just inside each sphere's
surface the writer gives FLEUR's local-frame value with the harmonics taken at R x (columns
6, 7), the same with R transposed (8, 9), and the plane-wave sum at the point (10, 11). The
right frame matches the plane-wave sum to the LAPW matching's accuracy; the wrong one does
not, unless R is symmetric (the identity, or an inversion partner's -1), where the two agree
and the test is empty, which is said. Per k point and atom: max |psi_R - psi_pw| and
max |psi_RT - psi_pw| over the surface points and bands, relative to max |psi_pw|, and
max |psi_R - psi_RT| over all sphere points.
Usage: check_probe.py <lapwfn_points.out> <lapwfn_points.in> <NM> <NS> <NSF> <nat>
"""
import sys
import numpy as np

fpts, fin = sys.argv[1], sys.argv[2]
NM, NS, NSF, nat = (int(x) for x in sys.argv[3:7])
npts = int(open(fin).readline())
kind = np.array(sum([['mesh'] * NM + ['sphere'] * NS + ['surface'] * NSF for _ in range(nat)], []) + ['inter'] * (npts - nat * (NM + NS + NSF)))
rows = {}
for line in open(fpts):
    if line.startswith('#'):
        continue
    t = line.split()
    ik, ib, ip, at = int(t[0]), int(t[1]), int(t[2]), int(t[3])
    v = [complex(float(t[5 + 2 * j]), float(t[6 + 2 * j])) for j in range(3)]
    rows.setdefault((ik, at), []).append((ip, ib, v))
print('# ik atom | surface points: max|R-pw|/max|pw|  max|RT-pw|/max|pw| | all sphere points: max|R-RT|/max|R|')
worst_r, worst_t, worst_d = 0.0, 0.0, 0.0
bad = []
for (ik, at) in sorted(rows):
    if at == 0:
        continue
    L = rows[(ik, at)]
    surf = [x for x in L if kind[x[0] - 1] == 'surface']
    if not surf:
        continue
    pw = np.array([x[2][2] for x in surf]); r = np.array([x[2][0] for x in surf]); rt = np.array([x[2][1] for x in surf])
    ref = np.abs(pw).max()
    dr, dt = np.abs(r - pw).max() / ref, np.abs(rt - pw).max() / ref
    ra = np.array([x[2][0] for x in L]); ta = np.array([x[2][1] for x in L])
    dd = np.abs(ra - ta).max() / np.abs(ra).max()
    note = '' if dd > 1e-10 else '  (R symmetric: the two frames coincide, no test)'
    print(f'{ik:3d} {at:4d} | {dr:.3e}  {dt:.3e} | {dd:.3e}{note}')
    if dr > 2e-2:
        bad.append(f'k {ik} atom {at}: R {dr:.2e} > 2e-2')
    if dd > 1e-10 and dt < 5 * dr:
        bad.append(f'k {ik} atom {at}: R transposed {dt:.2e} not above 5 x R {dr:.2e}')
    worst_r = max(worst_r, dr); worst_t = max(worst_t, dt); worst_d = max(worst_d, dd)
print(f'# worst over k and atoms: R {worst_r:.3e}, R transposed {worst_t:.3e}, |R - RT| {worst_d:.3e}')
# the right frame matches the plane-wave sum to the lmax truncation (below 2e-2 of the maximum in every
# case here); where the frames differ the transposed one must miss by at least five times as much
if bad:
    print('RESULT: FAIL ' + '; '.join(bad))
    sys.exit(1)
print('RESULT: PASS')
