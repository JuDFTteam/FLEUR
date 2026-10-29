"""Readers for the files a wannierlib run writes.

The operator exports and the interpolation output are FLEUR's own formats, not XML, so the
testset has to parse them itself. They live here rather than in the test module so that
test_wannier.py reads as what it is -- the run plus the claims made about it.
"""
import re


def spin_sumrule_values(path):
    """The |<s>| column of every 'spin sum-rule check' block FLEUR prints."""
    txt = open(path, errors="ignore").read()
    out = []
    for blk in re.findall(r"wannierlib spin sum-rule check, k = \d+\n.*?\n"
                          r"((?:\s+\d+\s+[-\d.]+.*\n)+)", txt):
        for line in blk.strip().split("\n"):
            p = line.split()
            if len(p) == 6:
                out.append(float(p[5]))
    return out


def outxml_eigenvalues(path):
    """{(spin, ikpt): [eigenvalues, Htr]} from the <eigenvaluesAt> blocks of out.xml."""
    import re
    out = {}
    for m in re.finditer(
            r'<eigenvaluesAt spin="(\d+)" ikpt="(\d+)"[^>]*>(.*?)</eigenvaluesAt>',
            open(path).read(), re.S):
        out[(int(m.group(1)), int(m.group(2)))] = [float(x) for x in m.group(3).split()]
    return out


def dat_rows(path):
    """Numeric rows of a bands_wann_*.dat: kdist first, then the per-band payload."""
    rows = []
    for line in open(path):
        if line.lstrip().startswith("#"):
            continue
        f = line.split()
        if len(f) > 1:
            rows.append([float(x) for x in f])
    return rows


def rspauli_r0_diagonal_sums(path):
    """Per-component sum of Re O_nn over the R=0 diagonal of a 'generic'-format O(R) file."""
    tot = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 8:
                continue
            if f[0] == f[1] == f[2] == "0" and f[3] == f[4]:
                tot[f[5]] = tot.get(f[5], 0.0) + float(f[6])
    return tot


def anglmom_r0(path):
    """R=0 block of a 'generic'-format O(R) file as {comp: {(i,j): complex}}."""
    blocks = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 8 or not (f[0] == f[1] == f[2] == "0"):
                continue
            blocks.setdefault(f[5], {})[(int(f[3]), int(f[4]))] = complex(
                float(f[6]), float(f[7]))
    return blocks


def anglmom_r0_traces(path):
    """Per-component trace of the R=0 block: sum_n <w_0n|L|w_0n>."""
    return {c: sum(v.real for (i, j), v in b.items() if i == j)
            for c, b in anglmom_r0(path).items()}


def anglmom_r0_hermiticity(path):
    """Largest |L_ij - conj(L_ji)| over the R=0 block, worst component."""
    worst = 0.0
    for b in anglmom_r0(path).values():
        for (i, j), v in b.items():
            w = b.get((j, i))
            if w is not None:
                worst = max(worst, abs(v - w.conjugate()))
    return worst


def rspauli_r0_diagonal_max(path):
    """Largest |Re O_nn| over the R=0 diagonal of a 'generic'-format O(R) file.
    Layout (see m_matrixelement_io): R1 R2 R3  i j comp  Re Im."""
    worst = 0.0
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 8:
                continue
            if f[0] == f[1] == f[2] == "0" and f[3] == f[4]:
                worst = max(worst, abs(float(f[6])))
    return worst


def as_tuple(v):
    """A reference is either one number or one per wannierised spin channel."""
    return v if isinstance(v, tuple) else (v,)


def last_n(values, n):
    """The last n matches, so a per-iteration echo of the same label cannot shift them."""
    vals = values if isinstance(values, list) else [values]
    assert len(vals) >= n, f"expected {n} values in the output, found {len(vals)}"
    return vals[-n:]


def onsite_weight(path):
    """Share of |O(R)|^2 that sits at R = 0, in per cent.

    The export has one row per (R, i, j, component) ending in Re and Im, so the weight is
    read off the file without knowing which operator wrote it."""
    total = onsite = 0.0
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 7:
                continue
            a = float(f[-2])**2 + float(f[-1])**2
            total += a
            if f[0] == f[1] == f[2] == "0":
                onsite += a
    return 100.0 * onsite / total if total > 0.0 else 0.0


def nonzero_entries(path):
    """Number of entries of an O(R) file that are not exactly zero.

    Reads the real and imaginary parts as the last two fields rather than at a fixed
    column: the number of index columns before them is not the same in every file, since
    a spinor operator is indexed by two spin labels where a vector operator carries one
    component label."""
    n = 0
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 8:
                continue
            if float(f[-2]) != 0.0 or float(f[-1]) != 0.0:
                n += 1
    return n


def torque_perp_residual(path, axis):
    """max |T(R).m_hat| / max |T(R)| over every entry of a three-component O(R) file.

    T = -(s x m_hat), so this vanishes identically. The floor is set by the file, not
    by the physics: the operator is written with eight decimals, so an entry whose
    largest component is |T|max carries a rounding error of 1e-8 and the ratio cannot
    go below 1e-8 / |T|max.
    """
    comps = {}
    for line in open(path):
        c = line.split()
        if len(c) != 8:
            continue
        comps.setdefault(tuple(c[:5]), {})[int(c[5])] = complex(float(c[6]), float(c[7]))
    if not comps:
        raise AssertionError("%s has no data rows" % path)
    scale = max(abs(v) for entry in comps.values() for v in entry.values())
    if scale == 0.0:
        raise AssertionError("%s is identically zero" % path)
    worst = max(abs(sum(entry.get(i + 1, 0.0) * axis[i] for i in range(3)))
                for entry in comps.values())
    return worst / scale


def perpmag_census(path):
    """(nbnd, nkpts, nwfs, rows, complete, trace, peak) of an updown.perpmag file.

    `trace` is sum_k Tr(M M+). It is invariant under a unitary rotation within either
    spin's degenerate multiplets, but that is not the only freedom a run has: the reduction
    order reaches the density, so the value moves by about 1e-3 relative with the compiler,
    the rank count and the thread count alike. Its caller pins it loosely for that reason.
    """
    with open(path) as fh:
        fh.readline()
        nb, nk, nw = (int(x) for x in fh.readline().split())
        seen = set()
        trace = 0.0
        peak = 0.0
        rows = 0
        for line in fh:
            c = line.split()
            if len(c) != 5:
                continue
            i, j, k = int(c[0]), int(c[1]), int(c[2])
            v = complex(float(c[3]), float(c[4]))
            seen.add((i, j, k))
            trace += abs(v) ** 2
            peak = max(peak, abs(v))
            rows += 1
    return nb, nk, nw, rows, len(seen) == nb * nk * nw, trace, peak
