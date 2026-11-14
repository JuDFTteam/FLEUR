"""
Differential regression test for <operatorsKs>, the Bloch-basis operator export.

What it guards is a JOIN, not a number. <operatorsKs> writes O(k) before any gauge
exists; <operatorsRs> writes O(R) after the gauge has been applied and Fourier
transformed. The two share the coarse matrices and nothing else, so if the Bloch
export ever starts writing a different object -- a transposed one, a differently
distributed one, a stale slice -- the identity below stops holding while every
self-contained check on the new file keeps passing.

That is the whole point: the file is Hermitian, has the right norm, the right
eigenvalues and the right shape whether or not its two band axes are the way round
the header claims, so no test that looks only at WF1_operators_ks.hdf can see the
error. This one closes the loop through a second, already-trusted writer:

    O(R = 0) = (1/Nk) sum_k V(k)^dag O(k) V(k),     V(k) = u_opt(k) . u_matrix(k)

with V read back out of the Wannier90 checkpoint the same run wrote. Measured
residual on this case is ~7e-9, which is the precision rspauli.1 is printed with.

NOT RUN BY CI, and it needs no marker of its own to stay out: the test venv installs
neither numpy nor h5py, so the importorskip below turns it into a skip there before
anything is executed. It is a cross-check of two writers against each other, not a
gate on the pipeline. Where both packages are present it runs like any other test:

    ./run_tests.sh testing/tests/invariants/test_wannier_bloch_invariants.py
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


@pytest.mark.wannierlib
@pytest.mark.bulk
@pytest.mark.soc
def test_bloch_operators_reproduce_the_real_space_route(execute_fleur):
    """Rotating <operatorsKs> with the checkpoint gauge must give <operatorsRs> at R = 0."""
    #> These two stay INSIDE the test and must not be hoisted to module level: they are
    #> what keeps it off CI. Collection imports the module, so a module-level "import
    #> numpy" raises there on a machine that lacks it -- which is every CI runner -- and
    #> fails the whole run. Inside the function the miss is a skip, which is the point.
    np = pytest.importorskip("numpy", reason="the Bloch cross-check needs numpy")
    h5py = pytest.importorskip("h5py", reason="the Bloch cross-check needs h5py")

    #> Two ranks rather than one on purpose: with a single rank the per-rank gather in
    #> the exporter is the identity and a mis-ordered k map would never show.
    res = execute_fleur(BASE, rm_files=["."], mpi_procs=2)

    for name in ("WF1_operators_ks.hdf", "fleur_wlib_internal.chk",
                 "rspauli.1", "anglmomrs.1"):
        assert name in res, (
            f"{name} was not written. This case asks for <operatorsKs>, <operatorsRs> "
            f"and a checkpoint, and the cross-check needs all three")

    num_bands, num_wann, num_kpts, gauge = _read_checkpoint(res["fleur_wlib_internal.chk"], np)

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

        #> The b-shell completeness relation, sum_b w_b b_a b_b = delta_ab. It is the one
        #> statement about the shells that is true in the right units and convention and
        #> false in any other, and the shells travel in this file precisely because no
        #> reader can rebuild them.
        wb = np.asarray(handle["bmesh_wb"])
        bk = np.asarray(handle["bmesh_bk"])
        completeness = np.abs(np.einsum("b,bi,bj->ij", wb, bk[0], bk[0]) - np.eye(3)).max()
        assert completeness < 1e-8, (
            f"the exported b-shells do not satisfy sum_b w_b b_a b_b = delta_ab "
            f"(off by {completeness:.3e}); the weights or the units moved")

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
