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
