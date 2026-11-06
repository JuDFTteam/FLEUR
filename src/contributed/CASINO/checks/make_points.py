#!/usr/bin/env python3
# Copyright (c) 2026 Mike Towler. MIT licence; see LICENSE in src/contributed/CASINO.
"""Writes lapwfn_points.in for the writer's point check from a first lapwfn.data of the same
case: per atom, NM points on exact mesh radii (indices beyond twice every join radius,
random directions), NS random points in the sphere between 0.05 bohr and the radius, NSF
random points just inside the sphere surface (R(1 - 1e-6), for the frame test against the
plane-wave sum), and NI random interstitial points of the cell (outside every sphere,
further than 1e-3 bohr from every surface).
Usage: make_points.py <lapwfn.data> <out file> [NM NS NSF NI seed]
"""
import os
import sys
import numpy as np
W = os.environ.get('WIEN2KQMC', os.path.expanduser('~/PROJECTS/wien2kqmc'))
sys.path.insert(0, W)
from wien2kqmc.lapwfn import read_lapwfn

fn, out = sys.argv[1], sys.argv[2]
NM, NS, NSF, NI, seed = (int(x) for x in sys.argv[3:8]) if len(sys.argv) > 7 else (30, 100, 20, 200, 1)
rng = np.random.default_rng(seed)
lapw = read_lapwfn(fn)
lat = np.asarray(lapw['lattice'])
pts = []
for at in lapw['atoms']:
    t = lapw['types'][at['type'] - 1]
    r = t['r']
    rjmax = max(fu['rjoin'] for l in range(t['lmax'] + 1) for fu in t['functions'][l])
    ok = np.where(r > 2 * rjmax)[0]
    idx = rng.choice(ok[:-1], NM, replace=False)
    for i in idx:
        d = rng.normal(size=3); d /= np.linalg.norm(d)
        pts.append(at['pos'] + r[i] * d)
    for _ in range(NS):
        d = rng.normal(size=3); d /= np.linalg.norm(d)
        rr = 0.05 + (t['rmt'] - 0.05) * rng.random() ** (1 / 3)
        pts.append(at['pos'] + rr * d)
    for _ in range(NSF):
        d = rng.normal(size=3); d /= np.linalg.norm(d)
        pts.append(at['pos'] + t['rmt'] * (1 - 1e-6) * d)
n = 0
while n < NI:
    f = rng.random(3)
    x = f @ lat
    inside = False
    for at in lapw['atoms']:
        t = lapw['types'][at['type'] - 1]
        d = x - at['pos']
        fr = np.linalg.solve(lat.T, d); fr -= np.round(fr)
        for sh in np.ndindex(3, 3, 3):
            p = (fr + np.array(sh) - 1) @ lat
            if np.linalg.norm(p) < t['rmt'] + 1e-3:
                inside = True
    if not inside:
        pts.append(x); n += 1
with open(out, 'w') as f:
    f.write(f'{len(pts)}\n')
    for p in pts:
        f.write(' %.15e %.15e %.15e\n' % tuple(p))
print(f'{len(pts)} points written to {out}: per atom {NM} on mesh radii, {NS} random in the sphere and {NSF} at the surface; {NI} interstitial')
