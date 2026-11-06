#!/usr/bin/env python3
"""The k points CASINO's simulation cell needs, as an explicit FLEUR k-point list. Usage:
   fleurqmc_klist.py M11 M12 M13 M21 M22 M23 M31 M32 M33 [t1 t2 t3] [--inp inp.xml] [--kpts kpts.xml]
                     [--name qmc] [--numbands N] [--etot E] [--title TEXT] [--write FILE]
The rows of M are the supercell vectors in the primitive basis (integers), t the twist in
fractions of the supercell's reciprocal basis (default 0). Prints the det(M) points in
fractions of the primitive reciprocal basis (FLEUR's kPoint coordinates), writes them as a
kPointList element to FILE (default fleurqmc_kpts.xml), and with --inp and --kpts puts the
list into kpts.xml, selects it in inp.xml, sets numbands to N + 1 and adds (or replaces) the
element <lapwfn totalEnergy="E" title="TEXT"/> in inp.xml's output section: the inp.xml of a
writer run on the SCF's converged density (run fleur in a copy of the SCF directory). E is the
converged total energy (totalEnergy in the SCF's out.xml); without --etot the element is left
alone and must be added by hand. At most N bands are written at a k point: the writer never writes
the highest band FLEUR computed, whose degenerate level cannot be known complete, and ends the list
at the last level boundary below it (so a k point whose Nth band is degenerate with the (N+1)th
gets fewer). FLEUR raises numbands to its own floor, max(5, 0.75 zelec + 1), when asked for less.
Every point gets weight 1; the weights play no part in the writer.
Copyright (c) 2026 Mike Towler. MIT licence; see LICENSE in src/contributed/CASINO.
"""
import re
import sys
from fractions import Fraction

argv = sys.argv[1:]
args, opts, i = [], {}, 0
while i < len(argv):
    a = argv[i]
    if a.startswith('--'):
        if a not in ('--inp', '--kpts', '--name', '--numbands', '--etot', '--title', '--write'):
            raise SystemExit(f'unknown option {a}\n' + __doc__)
        if i + 1 >= len(argv):
            raise SystemExit(f'{a} needs a value')
        opts[a] = argv[i + 1]
        i += 2
    else:
        args.append(a)
        i += 1
out = opts.get('--write', 'fleurqmc_kpts.xml')
name = opts.get('--name', 'qmc')
if len(args) not in (9, 12):
    raise SystemExit(__doc__)
M = [[int(x) for x in args[3*i:3*i+3]] for i in range(3)]
tw = [Fraction(x) for x in args[9:12]] or [Fraction(0)] * 3
Mf = [[Fraction(M[i][j]) for j in range(3)] for i in range(3)]


def det3(A):
    return (A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]) - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
            + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]))


def inv3(A):
    d = det3(A)
    C = [[None]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            m = [[A[r][c] for c in range(3) if c != j] for r in range(3) if r != i]
            C[i][j] = (-1)**(i+j) * (m[0][0]*m[1][1]-m[0][1]*m[1][0])
    return [[C[j][i]/d for j in range(3)] for i in range(3)]


n = abs(int(det3(Mf)))
if n == 0:
    raise SystemExit('singular supercell matrix')
# exact: the rows of M are the supercell vectors S_i = sum_j M_ij a_j, so a k point of the
# supercell's Gamma (plus the twist t) satisfies k.S_i = 2 pi (m_i + t_i), that is M k = m + t
# with k in fractions of the primitive reciprocal basis: k = M^-1 (m + t), reduced into [0, 1)
Minv = inv3(Mf)
ks, seen = [], set()
rng = range(-2*n, 2*n+1)
for m1 in rng:
    for m2 in rng:
        for m3 in rng:
            mm = [tw[0]+m1, tw[1]+m2, tw[2]+m3]
            kf = tuple((sum(Minv[i][j]*mm[j] for j in range(3))) % 1 for i in range(3))
            if kf not in seen:
                seen.add(kf)
                ks.append(kf)
            if len(ks) == n:
                break
        if len(ks) == n:
            break
    if len(ks) == n:
        break
if len(ks) != n:
    raise SystemExit(f'found {len(ks)} points for a supercell of {n} cells')
for k in ks:     # the check of the condition itself: M k - t is a vector of integers
    for i in range(3):
        v = sum(Mf[i][j]*k[j] for j in range(3)) - tw[i]
        if v.denominator != 1:
            raise SystemExit(f'the point {k} fails M k = m + t (row {i + 1} gives {v})')
print(f'{n} k points (fractions of the primitive reciprocal basis, FLEUR kPoint coordinates):')
for k in ks:
    print('  ' + ' '.join(f'{float(x):.12f}' for x in k) + '   = ' + ' '.join(str(x) for x in k))
lines = [f'            <kPointList name="{name}" count="{n}">']
for k in ks:
    lines.append('               <kPoint weight="1.0">   ' + '    '.join(f'{float(x):.12f}' for x in k) + '</kPoint>')
lines.append('            </kPointList>')
block = '\n'.join(lines) + '\n'
with open(out, 'w') as f:
    f.write(block)
print(f'wrote {out}')
if '--kpts' in opts:
    p = opts['--kpts']
    s = open(p).read()
    s = re.sub(r'\n *<kPointList name="%s"[^>]*>.*?</kPointList>' % re.escape(name), '', s, flags=re.S)
    if '</kPointLists>' not in s:
        raise SystemExit(f'{p} has no </kPointLists>')
    s = s.replace('</kPointLists>', block + '         </kPointLists>', 1)
    open(p, 'w').write(s)
    print(f'{p}: the list "{name}" added')
if '--inp' in opts:
    p = opts['--inp']
    s = open(p).read()
    s, c1 = re.subn(r'<kPointListSelection listName="[^"]*"', f'<kPointListSelection listName="{name}"', s)
    if c1 != 1:
        raise SystemExit(f'{p}: kPointListSelection not found once ({c1})')
    c3 = 0
    if '--numbands' in opts:
        s, c3 = re.subn(r'numbands="[^"]*"', f'numbands="{int(opts["--numbands"]) + 1}"', s)
    msg = f'{p}: list "{name}" selected'
    if '--etot' in opts:
        float(opts['--etot'])
        el = f'<lapwfn totalEnergy="{opts["--etot"]}"' + (f' title="{opts["--title"]}"' if '--title' in opts else '') + '/>'
        s = re.sub(r'\n *<lapwfn\b[^>]*/>', '', s)
        s, c2 = re.subn(r'(<output\b[^>]*>)', lambda m: m.group(1) + '\n      ' + el, s, count=1)
        if c2 != 1:
            raise SystemExit(f'{p}: no <output> element found')
        msg += ', ' + el + ' added'
    open(p, 'w').write(s)
    print(msg + (f', numbands {int(opts["--numbands"]) + 1} (at most {opts["--numbands"]} bands written: the highest computed band is never written)' if c3 else ''))
