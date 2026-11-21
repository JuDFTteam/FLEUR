"""
Differential regression test for <operatorsKs>, in both of its halves.

What it guards is a JOIN, not a number. <operatorsKs> hands over what the run holds in
k; <operatorsRs> writes O(R) after the Fourier sum. The two share their inputs and
nothing else, so if either export starts writing a different object -- a transposed
one, a differently distributed one, a stale slice -- the identities below stop holding
while every self-contained check on the new file keeps passing.

The block has two halves and this covers one join per half:

  spin, orbital       Bloch basis, no gauge. The gauge is applied here, from the
                      checkpoint, and must land on O(R = 0).
  bmn                 already contracted over neighbours, so it carries the gauge
                      and all that is left between it and the .dat is a phase sum.

They fail for different reasons -- the first for a gauge or a band-axis error, the
second for an R set, a component order or a unit -- which is why both are here.

WHY THE MESH IS 4x4x4, AND WHY IT MUST NOT BE SHRUNK. The second join is a phase sum,
and on a 2x2x2 mesh every k component is 0 or 1/2, so every k.R is a multiple of 1/2
and every phase is +-1: real. The sign of the exponent is then not observable at all --
flipping it leaves the result bit-identical, and so does sending R to -R. That is not a
weak assertion but a case with no information in it, and it was measured: on 2x2x2 the
inverted-sign mutation PASSES, on 4x4x4 it fails by 1.46 against a tolerance of 2e-6.
Fifty-six more k-points are what buys this test its sensitivity, and dropping back to
the coarser mesh would blind it without turning anything red.

That is the whole point: the file is Hermitian, has the right norm, the right
eigenvalues and the right shape whether or not its two band axes are the way round
the header claims, so no test that looks only at WF1_operators_ks.hdf can see the
error. This one closes the loop through a second, already-trusted writer:

    O(R = 0) = (1/Nk) sum_k V(k)^dag O(k) V(k),     V(k) = u_opt(k) . u_matrix(k)

with V read back out of the Wannier90 checkpoint the same run wrote. Measured
residual on this case is ~7e-9, which is the precision rspauli.1 is printed with.

The same case also asks for <export bmesh="T"/>. Nothing else writes that file, so
there is no second writer to join it to; what stands in for one is its layout, pinned
shape by shape, and the two identities the shells have to obey for a finite difference
over them to be a centred derivative.

NOT RUN BY CI, and it needs no marker of its own to stay out: the test venv installs
neither numpy nor h5py, so the importorskip below turns it into a skip there before
anything is executed. It is a cross-check of two writers against each other, not a
gate on the pipeline. Where both packages are present it runs like any other test:

    ./run_tests.sh -k operators_ks_reproduces
"""

import pytest

BASE = "./inputfiles/wannier/WannFeBccKs"

#> Both files are printed with eight decimals, so agreement cannot be asked for below
#> that. Anything that actually breaks the join is O(1), not O(1e-7).
TOL = 1e-7


def _records(path):
    """Fortran unformatted sequential: every record sits between two 4-byte markers."""
    blob = open(path, "rb").read()
    out, i = [], 0
    while i < len(blob):
        n = int.from_bytes(blob[i:i + 4], "little")
        out.append(blob[i + 4:i + 4 + n])
        i += 4 + n + 4
    return out


def _read_checkpoint(path, np):
    """The gauge Wannier90 stored: V(k) = u_opt(k) . u_matrix(k).

    Record order is Wannier90's own write_chkpt. Only the fields needed here are
    decoded; the rest are stepped over by position, which is why the two optional
    blocks have to be counted rather than assumed.
    """
    rec = _records(path)
    as_int = lambda k: int(np.frombuffer(rec[k], "<i4")[0])
    num_bands, num_kpts, num_wann = as_int(1), as_int(6), as_int(10)
    have_disentangled = bool(as_int(12))

    k = 13
    u_opt = None
    if have_disentangled:
        k += 1                                   # omega_invariant
        k += 1                                   # lwindow
        k += 1                                   # ndimwin
        #> Fortran writes u_matrix_opt(band, wann, k) with the band index fastest, so the
        #> flat buffer reshapes (k, wann, band) and the two matrix axes swap back here.
        u_opt = np.frombuffer(rec[k], "<c16").reshape(num_kpts, num_wann, num_bands)
        u_opt = u_opt.transpose(0, 2, 1)
        k += 1
    u_mat = np.frombuffer(rec[k], "<c16").reshape(num_kpts, num_wann, num_wann)
    u_mat = u_mat.transpose(0, 2, 1)

    gauge = u_mat if u_opt is None else np.einsum("kbw,kwv->kbv", u_opt, u_mat)
    return num_bands, num_wann, num_kpts, gauge


def _read_operator_r(path, num_wann, np):
    """O(R) at R = 0 from the 'generic' layout: R1 R2 R3  i j comp  Re Im, no header."""
    raw = np.loadtxt(path)
    at_origin = (raw[:, 0] == 0) & (raw[:, 1] == 0) & (raw[:, 2] == 0)
    assert at_origin.any(), f"{path} holds no R = 0 block"
    rows = raw[at_origin]
    ncomp = int(raw[:, 5].max())
    out = np.zeros((ncomp, num_wann, num_wann), complex)
    out[rows[:, 5].astype(int) - 1,
        rows[:, 3].astype(int) - 1,
        rows[:, 4].astype(int) - 1] = rows[:, 6] + 1j * rows[:, 7]
    return out


def _read_operator_k(group, np):
    """One operator out of WF<n>_operators_ks.hdf, as (nk, ncomp, row, column).

    HDF reverses the dimension order, so the Fortran (2, nb, nb, ncomp, nk) arrives as
    (nk, ncomp, nb, nb, 2) with the two band axes SWAPPED -- axis 2 is the column and
    axis 3 the row. Undoing that is the one step nothing downstream can check for you:
    the transpose of a Hermitian matrix is Hermitian, keeps its norm and keeps its
    eigenvalues, and answers a different question in silence.
    """
    block = np.asarray(group)
    return np.swapaxes(block[..., 0] + 1j * block[..., 1], -1, -2)


def _check_bmesh(path, np, h5py):
    """WF<n>_bmesh.hdf, on its own terms: the layout and the two shell identities.

    THE LAYOUT IS THE POINT. Nothing else writes this file, so there is no second
    writer to join it to and the shapes are the contract. HDF reverses the dimension
    order, so every array is read back to front relative to the Fortran that wrote it:
    nnlist goes out (nk, nntot) and arrives (nntot, nk), gkpb goes out (3, nk, nntot)
    and arrives (nntot, nk, 3), bk goes out (3, nntot, nk) and arrives (nk, nntot, 3).
    They are pinned one by one below, because a reader that guesses wrong gets arrays
    of a plausible shape and silently pairs the wrong b with the wrong k.

    THE IDENTITIES. sum_b w_b b_a b_c = delta_ac is what makes a finite difference over
    the shells a derivative, and sum_b w_b b = 0 is what makes it centred. Both hold in
    the units the file is written in and in no other, so between them they catch a unit,
    a weight or a convention moving. Both are checked at EVERY k and not just the first:
    Wannier90 chooses the shells per k-point and bk carries a k index for that reason.

    WHAT IS NOT CHECKED, and why. Not that neighbour index b names the same b vector for
    every k: which b a given slot holds is Wannier90's business, and bk carries a k index
    for exactly that reason. Tying nnlist and gkpb to bk needs the reciprocal lattice,
    which this file does not carry and a consumer has anyway.
    """
    with h5py.File(path, "r") as handle:
        nntot = int(np.asarray(handle["nntot"])[0])
        wb = np.asarray(handle["wb"])
        bk = np.asarray(handle["bk"])
        nnlist = np.asarray(handle["nnlist"])
        gkpb = np.asarray(handle["gkpb"])
        bkf = np.asarray(handle["kpts"])

    nk = bkf.shape[0]
    for name, got, want in (("wb", wb.shape, (nntot,)),
                            ("bk", bk.shape, (nk, nntot, 3)),
                            ("nnlist", nnlist.shape, (nntot, nk)),
                            ("gkpb", gkpb.shape, (nntot, nk, 3)),
                            ("kpts", bkf.shape, (nk, 3))):
        assert got == want, (
            f"{name} came back {got} and the layout says {want} for {nk} k-points and "
            f"{nntot} neighbours; the index order of the export moved")

    completeness = np.abs(np.einsum("b,kbi,kbj->kij", wb, bk, bk) - np.eye(3)[None]).max()
    assert completeness < 1e-8, (
        f"the exported b-shells do not satisfy sum_b w_b b_a b_c = delta_ac "
        f"(off by {completeness:.3e}); the weights or the units moved")

    centred = np.abs(np.einsum("b,kba->ka", wb, bk)).max()
    assert centred < 1e-8, (
        f"the shells are not centred, sum_b w_b b = {centred:.3e}; the difference they "
        f"define is first order in b and not a derivative")

    #> 1-based, Fortran's own. An out-of-range entry would index silently in python.
    assert nnlist.min() >= 1 and nnlist.max() <= nk, (
        f"nnlist runs {nnlist.min()}..{nnlist.max()} and must index the {nk} stored "
        f"k-points, 1-based")
    assert np.abs(gkpb).max() > 0, (
        "every gkpb is zero, so no neighbour leaves the zone. On a full-zone mesh some "
        "must, and a gkpb that is always zero is the signature of the field never being "
        "filled rather than of a mesh that does not need it")

    #> Handed back because the phase sum below needs the very k-list the file was built
    #> on, and reconstructing it from inp.xml would be a second convention to get wrong.
    return bkf


def _read_bmn_dat(path, num_wann, np):
    """WF<n>_bmn.dat in the 'bmn' layout: R1 R2 R3  i j  then the three components.

    One row per (R, i, j) with all three Cartesian components on it, so a row is
    5 + 2*3 = 11 numbers. The file opens with a header line that happens to split into
    eleven words too, which is why the first five are also required to be integers --
    counting the columns alone lets the header through and poisons the R set.

    Returns the R vectors as written and B(R) as (comp, row, column, R).
    """
    rvec, rows = [], []
    with open(path) as handle:
        for line in handle:
            tok = line.split()
            if len(tok) != 11:
                continue
            try:
                head = [int(x) for x in tok[:5]]
            except ValueError:
                continue
            rvec.append(head[:3])
            rows.append((head, [float(x) for x in tok[5:]]))
    assert rows, f"{path} holds no data row in the 'bmn' layout"

    uniq = sorted({tuple(r) for r in rvec})
    index = {r: i for i, r in enumerate(uniq)}
    out = np.zeros((3, num_wann, num_wann, len(uniq)), complex)
    for head, vals in rows:
        ir = index[tuple(head[:3])]
        for a in range(3):
            out[a, head[3] - 1, head[4] - 1, ir] = vals[2 * a] + 1j * vals[2 * a + 1]
    return np.array(uniq, float), out


@pytest.mark.wannierlib
@pytest.mark.bulk
@pytest.mark.soc
def test_operators_ks_reproduces_the_real_space_route(execute_fleur):
    """Both halves of <operatorsKs> must land on what <operatorsRs> wrote."""
    #> These two stay INSIDE the test and must not be hoisted to module level: they are
    #> what keeps it off CI. Collection imports the module, so a module-level "import
    #> numpy" raises there on a machine that lacks it -- which is every CI runner -- and
    #> fails the whole run. Inside the function the miss is a skip, which is the point.
    np = pytest.importorskip("numpy", reason="the Bloch cross-check needs numpy")
    h5py = pytest.importorskip("h5py", reason="the Bloch cross-check needs h5py")

    #> Two ranks rather than one on purpose: with a single rank the per-rank gather in
    #> the exporter is the identity and a mis-ordered k map would never show.
    res = execute_fleur(BASE, rm_files=["."], mpi_procs=2)

    for name in ("WF1_operators_ks.hdf", "WF1_modern_k.hdf", "WF1_bmesh.hdf", "WF1.chk",
                 "rspauli.1", "anglmomrs.1", "WF1_bmn.dat"):
        assert name in res, (
            f"{name} was not written. This case asks for <operatorsKs>, <operatorsRs>, "
            f"a b-mesh and a checkpoint, and the cross-check needs all of them")

    bkf = _check_bmesh(res["WF1_bmesh.hdf"], np, h5py)

    num_bands, num_wann, num_kpts, gauge = _read_checkpoint(res["WF1.chk"], np)

    #> If the checkpoint and the export disagree about the manifold there is nothing to
    #> compare, and saying so here beats a shape error twenty lines down.
    unitarity = np.abs(np.einsum("kbn,kbm->knm", np.conj(gauge), gauge)
                       - np.eye(num_wann)[None]).max()
    assert unitarity < 1e-10, (
        f"the checkpoint gauge is not unitary on the manifold (max |V^dag V - 1| = "
        f"{unitarity:.3e}); the cross-check below would be meaningless")

    with h5py.File(res["WF1_operators_ks.hdf"], "r") as handle:
        nb_file, nk_file = (int(v) for v in handle["dims"][:2])
        assert (nb_file, nk_file) == (num_bands, num_kpts), (
            f"the export says {nb_file} bands and {nk_file} k-points, the checkpoint "
            f"says {num_bands} and {num_kpts}")

        for dataset, real_space in (("spin", "rspauli.1"), ("orbital", "anglmomrs.1")):
            bloch = _read_operator_k(handle[dataset], np)
            mine = np.einsum("kbn,kcbd,kdm->cnm",
                             np.conj(gauge), bloch, gauge) / num_kpts
            theirs = _read_operator_r(res[real_space], num_wann, np)

            scale = np.abs(theirs).max()
            assert scale > 1e-3, f"{real_space} at R = 0 is all but empty (max {scale:.3e})"
            deviation = np.abs(theirs - mine).max() / scale
            assert deviation < TOL, (
                f"{dataset} in the Bloch export no longer reproduces {real_space} at "
                f"R = 0: relative deviation {deviation:.3e} against a tolerance of {TOL:.0e}. "
                f"The two writers share the coarse matrices, so either the gauge rotation, "
                f"the k ordering or the band-axis order of one of them has moved")

    #> The second half. B is already contracted over the neighbours and already gauged,
    #> so no rotation happens here: what separates the two files is the phase sum alone,
    #> over the R set the .dat itself carries and with no degeneracy weight -- the cart2
    #> family writes none and the sum takes none, unlike H(R).
    #>
    #> The tolerance is ABSOLUTE and does not scale: the .dat is printed with six
    #> decimals, so 5e-7 of disagreement is the file's own rounding whatever the size of
    #> B, and anything that breaks the join is O(1).
    irvec, theirs = _read_bmn_dat(res["WF1_bmn.dat"], num_wann, np)
    with h5py.File(res["WF1_modern_k.hdf"], "r") as handle:
        nw_file, nk_file = (int(v) for v in handle["dims"][:2])
        assert (nw_file, nk_file) == (num_wann, num_kpts), (
            f"the k-space tensors say {nw_file} Wannier functions and {nk_file} "
            f"k-points, the checkpoint says {num_wann} and {num_kpts}")
        bloch = _read_operator_k(handle["bmn"], np)          # (nk, 3, row, column)

    phase = np.exp(-2j * np.pi * (bkf @ irvec.T)) / num_kpts  # (nk, nR)
    mine = np.einsum("kcnm,kr->cnmr", bloch, phase)

    scale = np.abs(theirs).max()
    assert scale > 1e-3, f"WF1_bmn.dat is all but empty (max {scale:.3e})"
    deviation = np.abs(theirs - mine).max()
    assert deviation < 2e-6, (
        f"B(k) no longer Fourier-sums to WF1_bmn.dat: max deviation {deviation:.3e} "
        f"against a tolerance of 2e-6 on a scale of {scale:.3e}. Nothing but a phase sum "
        f"sits between the two, so the suspects are the R set, the order of the three "
        f"Cartesian components, the band-axis order of the HDF read, or the units")
