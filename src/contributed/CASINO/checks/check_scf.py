#!/usr/bin/env python3
# Copyright (c) 2026 Mike Towler. MIT licence; see LICENSE in src/contributed/CASINO.
"""The file against the SCF it came from (out.xml of the converged run): the eigenvalues at
Gamma (the first k point of a Gamma-centred mesh, the last iteration's), the Fermi energy, the
total energy line, the core energies and FLEUR's lost core charge against the writer's leakage
(lapwfn_measurements.out), and the electron counts.
Usage: check_scf.py <lapwfn.data> <out_scf.xml> <lapwfn_measurements.out>
"""
import re
import sys

L = open(sys.argv[1]).read().splitlines()
X = open(sys.argv[2]).read()
M = open(sys.argv[3]).read().splitlines()


def after(label):
    return L[L.index(label) + 1]


ef = float(after('Fermi energy (au)'))
et = float(after('Total energy (au per primitive cell)'))
nel = int(after('Number of electrons per primitive cell (core included)'))
nval = int(after('Number of valence electrons per primitive cell (in the WAVE FUNCTION section)'))
# the file's Gamma eigenvalues
eg = []
for i, s in enumerate(L):
    if s.startswith('k-point #'):
        t = L[i + 1].split()
        if all(abs(float(x)) < 1e-9 for x in t[3:6]):
            j = i + 2
            while j < len(L) and not L[j].startswith('k-point #'):
                if L[j].startswith('Band, spin'):
                    eg.append(float(L[j + 1].split()[2]))
                j += 1
            break
# the SCF's last iteration
its = X.split('<iteration ')
last = its[-1]
m = re.search(r'<eigenvaluesAt spin="1" ikpt="1"[^>]*>(.*?)</eigenvaluesAt>', last, re.S)
es = [float(x) for x in m.group(1).split()]
efs = float(re.search(r'<FermiEnergy value="\s*([-0-9.]+)"', last).group(1))
ets = float(re.search(r'<totalEnergy\s+value="\s*([-0-9.]+)"', last).group(1))
n = min(len(eg), len(es))
print("# Gamma eigenvalues: the file against the SCF's last iteration (Ha)")
for i in range(n):
    print(f'{i + 1:3d} {eg[i]:.10f} {es[i]:.10f} {eg[i] - es[i]:+.3e}')
bad = []
dmax = max(abs(eg[i] - es[i]) for i in range(n))
print(f'# max |diff| over {n} bands: {dmax:.3e}')
print(f'# Fermi energy: file {ef:.10f}, SCF {efs:.10f}, diff {ef - efs:+.3e}')
print(f'# total energy: file {et:.10f}, SCF {ets:.10f}, diff {et - ets:+.3e}')
if dmax > 1e-6:
    bad.append(f'Gamma eigenvalues differ from the SCF by {dmax:.1e} Ha (bound 1e-6)')
if abs(ef - efs) > 1e-8:
    bad.append(f'Fermi energy differs by {ef - efs:.1e} Ha (bound 1e-8)')
if abs(et - ets) > 1e-9:
    bad.append(f'total energy differs by {et - ets:.1e} Ha (bound 1e-9)')
# core states
cs = re.findall(r'<coreStates atomType="\s*(\d+)"[^>]*lostElectrons="\s*([0-9.]+)">(.*?)</coreStates>', last, re.S)
fc = {}
for i, s in enumerate(L):
    if s.startswith('Core state: type'):
        t = L[i + 1].split()
        fc[(int(t[0]), int(t[2]), int(t[3]))] = (float(t[4]), float(t[5]))
leak = {}
for s in M:
    if s.startswith('core '):
        t = s.split()
        kv = dict(zip(t[6::2], t[7::2]))
        leak[(int(t[1]), int(t[2]), int(t[3]))] = (float(kv['leak_tail']), float(kv['leak_fleur']))
print('# core states: type n l | file energy, SCF energy (the last iteration, kappa-averaged), diff | file occupation | leakage: tail, FLEUR\'s own (per electron), SCF lostElectrons/electrons')
for typ, lost, body in cs:
    typ = int(typ)
    states = {}
    for nn, ll, en, wt in re.findall(r'<state n="(\d+)" l="(\d+)" j="[0-9.]+" energy="\s*([-0-9.]+)" weight="\s*([0-9.]+)"', body):
        k = (int(nn), int(ll))
        e, w = states.get(k, (0.0, 0.0))
        states[k] = (e + float(en) * float(wt), w + float(wt))
    ntot = sum(2 * w for e, w in states.values())
    for (nn, ll), (ew, w) in sorted(states.items()):
        e = ew / w
        f = fc.get((typ, nn, ll))
        lk = leak.get((typ, nn, ll))
        if f is None:
            print(f'{typ:3d} {nn:2d} {ll:2d} | SCF {e:.10f} | not in the file')
            bad.append(f'core state type {typ} n {nn} l {ll} not in the file')
            continue
        if abs(f[0] - e) > 1e-6 or abs(f[1] - 2 * w) > 1e-6:
            bad.append(f'core state type {typ} n {nn} l {ll}: energy or occupation differs from the SCF')
        print(f'{typ:3d} {nn:2d} {ll:2d} | {f[0]:.10f} {e:.10f} {f[0] - e:+.3e} | occupation {f[1]:.1f} (SCF {2 * w:.1f}) | '
              + (f'leak tail {lk[0]:.4e} fleur {lk[1]:.4e}' if lk else 'no leak line') + f' lost/electrons {float(lost) / ntot:.4e}')
print(f'# electrons: file {nel} (valence {nval})')
print('RESULT: ' + ('PASS' if not bad else 'FAIL ' + '; '.join(bad)))
sys.exit(1 if bad else 0)
