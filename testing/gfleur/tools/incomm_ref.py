"""NumPy reference for the incommensurate embedding transformer (gf_embtransform).

The reference evaluates the substrate amplitudes sigma_u(p) = <p+u|Sigma|p> from an
analytic function on the same k-mesh as the synthetic substrate file and interpolates
them with the same 4x4 Lagrange stencil as gf_sigtab. It builds

    Sigma0(i,j) = exp(-i u.d) sigma_u((kappa+G_j) M),   u = (G_i-G_j) M integer,

and the Umklapp correction

    Sigma_eff = Sigma0 - Delta_0Q [D_Q + Delta_QQ]^-1 Delta_Q0

independently of the Fortran code. It also writes the synthetic files in the layout
of gf_io2dmat. HDF5 datasets are stored in Fortran order, so h5py sees reversed
shapes: 2Dmat(ng,ng,2,spins,nE,nk) is arr[k,e,s,ri,g2,g1] for element (g1,g2).
"""
import numpy as np

try:
    import h5py
except ImportError:  # the tests skip without h5py
    h5py = None

TWO_PI = 2.0 * np.pi


# --------------------------------------------------------------------------- lattices
def hex_amat(a, angle_deg=0.0):
    """In-plane lattice vectors as columns, a1 along x rotated by angle."""
    t = np.radians(angle_deg)
    rot = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    return rot @ np.array([[a, 0.5 * a], [0.0, 0.5 * np.sqrt(3.0) * a]])


def recip(amat):
    """Reciprocal vectors as rows: B = 2 pi A^-1."""
    return TWO_PI * np.linalg.inv(amat)


def rotation(theta_deg):
    t = np.radians(theta_deg)
    return np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])


def point_group_ops(amat, cart_ops):
    """Integer real-space operations x -> mrot x for Cartesian 2x2 operations."""
    ops = []
    ainv = np.linalg.inv(amat)
    for r in cart_ops:
        m = ainv @ r @ amat
        assert np.allclose(m, np.rint(m), atol=1e-8)
        ops.append(np.rint(m).astype(int))
    return ops


def c6v_cartesian():
    ops = []
    for n in range(6):
        r = rotation(60.0 * n)
        ops.append(r)
        ops.append(r @ np.diag([1.0, -1.0]))
    return ops


def mesh(n1, n2):
    """Gamma-centered mesh as gfprep: k in (-1/2, 1/2]."""
    k = []
    for j in range(n2):
        for i in range(n1):
            q = np.array([i / n1, j / n2])
            k.append(q - np.rint(q))
    return np.array(k)


def irreducible(kpts, bmat, cart_ops, timerev=True):
    """Representatives of the orbits of the k-points (Cartesian operations + TR)."""
    binv = np.linalg.inv(bmat)
    left = list(range(len(kpts)))
    reps = []
    while left:
        i = left.pop(0)
        reps.append(kpts[i])
        kc = kpts[i] @ bmat
        images = []
        for r in cart_ops:
            for s in ([1, -1] if timerev else [1]):
                images.append((s * (r @ kc)) @ binv)
        keep = []
        for j in left:
            if any(np.allclose(kpts[j] - im, np.rint(kpts[j] - im), atol=1e-8) for im in images):
                continue
            keep.append(j)
        left = keep
    return np.array(reps)


def glist(bmat, rk):
    """G-list of gfleur (gf_apws_globalmap): |G+s| <= rk for one of the shifts
    s in {-1/2,0,1/2}^2, i.e. the G of |k+G| <= rk for k in the first zone"""
    n = int(rk / min(np.linalg.norm(bmat, axis=1))) + 4
    shifts = [np.array([a, b]) for a in (-0.5, 0.0, 0.5) for b in (-0.5, 0.0, 0.5)]
    g = [(i, j) for i in range(-n, n + 1) for j in range(-n, n + 1)
         if min(np.linalg.norm((np.array([i, j]) + s) @ bmat) for s in shifts) <= rk]
    g.sort(key=lambda x: (np.linalg.norm(np.array(x) @ bmat), x))
    return np.array(g, dtype=int)


# --------------------------------------------------------------------------- amplitudes
def sigma_model(e, spin):
    """Analytic substrate amplitudes: rotation and time-reversal invariant."""
    v0 = 0.3 + 0.05 * (spin - 1)
    amp = (0.06 + 0.03j) * (1.0 + 0.1 * e) * (1.0 - 0.1 * (spin - 1))

    def f(p2, pu2, u2, is_zero):
        diag = -0.5 * np.sqrt(p2 + 2.0 * (v0 - e) + 0j)
        off = amp * np.exp(-0.12 * (p2 + pu2) - 0.05 * u2)
        return np.where(is_zero, diag, off)
    return f


def sigma_model_noco(e, js):
    """noco substrate: file spin index js = 11, 22, 12, 21; spin off-diagonal amplitudes"""
    if js <= 2:
        return sigma_model(e, js)
    amp = (0.02 - 0.015j) if js == 3 else (0.01 + 0.02j)

    def f(p2, pu2, u2, is_zero):
        return amp * (1.0 + 0.2 * e) * np.exp(-0.15 * (p2 + pu2) - 0.05 * u2)
    return f


def soc_matrix(e, pc, puc):
    """Rashba-like 2x2 amplitudes from cartesian p and p' = p+u (arrays [..., 2]):
    f [a0 + i l1 (p x p')_z s_z + l2 ((p+p') x z).s], invariant under rotations about z,
    vertical mirrors (spin as axial vector) and s_y sigma^T s_y time reversal"""
    p2 = np.sum(pc ** 2, -1)
    pu2 = np.sum(puc ** 2, -1)
    f = np.exp(-0.12 * (p2 + pu2)) * (1.0 + 0.1 * e)
    zero = np.all(np.abs(puc - pc) < 1e-12, axis=-1)
    a0 = np.where(zero, -0.5 * np.sqrt(p2 + 2.0 * (0.3 - e) + 0j), (0.06 + 0.03j) * f)
    c = pc[..., 0] * puc[..., 1] - pc[..., 1] * puc[..., 0]
    big_p = pc + puc
    rx, ry = big_p[..., 1], -big_p[..., 0]
    l1, l2 = (0.02 + 0.01j) * f, (0.03 - 0.02j) * f
    m = np.empty(np.shape(c) + (2, 2), dtype=complex)
    m[..., 0, 0] = a0 + 1j * l1 * c
    m[..., 1, 1] = a0 - 1j * l1 * c
    m[..., 0, 1] = l2 * (rx - 1j * ry)
    m[..., 1, 0] = l2 * (rx + 1j * ry)
    return m


def sigma_model_soc(e, js):
    """file spin index js = 11, 22, 12, 21 of soc_matrix (kind 'subvec')"""
    i, j = {1: (0, 0), 2: (1, 1), 3: (0, 1), 4: (1, 0)}[js]

    def f(pc, puc, uc, is_zero):
        return soc_matrix(e, pc, puc)[..., i, j]
    return f


def c3v_cartesian():
    ops = []
    for n in range(3):
        r = rotation(120.0 * n)
        ops.append(r)
        ops.append(r @ np.diag([1.0, -1.0]))
    return ops


def sigma_vacuum(e, spin, vvac=0.2):
    def f(p2, pu2, u2, is_zero):
        return np.where(is_zero, -0.5 * np.sqrt(p2 + 2.0 * (vvac - e) + 0j), 0.0)
    return f


def up_model(e, spin):
    """Synthetic Sigma_up(q)(m1,m2), m = q+G Cartesian."""
    def f(m1, m2):
        d2 = np.sum((m1 - m2) ** 2, axis=-1)
        a1 = np.sum(m1 ** 2, axis=-1)
        a2 = np.sum(m2 ** 2, axis=-1)
        val = (0.03 - 0.01j) * np.exp(-0.2 * (a1 + a2) - 0.3 * d2) * (1.0 + 0.05 * spin)
        diag = -0.5 * np.sqrt(a1 + 2.0 * (0.1 - e) + 0j)
        return np.where(d2 < 1e-20, diag, val)
    return f


def up_model_noco(e, js):
    """noco Sigma_up: spin index js = 11, 22, 12, 21"""
    if js <= 2:
        return up_model(e, js)
    amp = (0.004 + 0.002j) if js == 3 else (0.003 - 0.001j)

    def f(m1, m2):
        a1 = np.sum(m1 ** 2, axis=-1)
        a2 = np.sum(m2 ** 2, axis=-1)
        d2 = np.sum((m1 - m2) ** 2, axis=-1)
        return amp * np.exp(-0.2 * (a1 + a2) - 0.3 * d2)
    return f


def up_model_diag(model):
    """noco Sigma_up from a collinear model: 11 = spin 1, 22 = spin 2, no spin coupling"""
    def m(e, js):
        if js <= 2:
            return model(e, js)
        return lambda m1, m2: np.zeros(np.broadcast(m1[..., 0], m2[..., 0]).shape)
    return m


# --------------------------------------------------------------------------- files
def _att(obj, name, value, dtype):
    obj.attrs.create(name, np.atleast_1d(np.asarray(value, dtype=dtype)))


def write_embfile(name, amat, kpts, gvecs, energies, spins, model, rkmax, ops=None,
                  timerev=False, model_kind="sub", soc=False):
    """Embedding-potential file in the gf_io2dmat layout.

    model(e, spin) returns f(p2, pu2, u2, is_zero) (kind 'sub', p=k+g2, u=g1-g2)
    or f(m1, m2) with Cartesian momenta (kind 'up').
    """
    bmat = recip(amat)
    ng, nk, ne = len(gvecs), len(kpts), len(energies)
    data = np.zeros((nk, ne, spins, 2, ng, ng))
    for ik, k in enumerate(kpts):
        m = (k[None, :] + gvecs) @ bmat
        inside = np.sum(m ** 2, axis=1) <= rkmax ** 2 * (1.0 + 1e-10)
        g1 = gvecs[:, None, :]
        g2 = gvecs[None, :, :]
        u = g1 - g2
        for ie, e in enumerate(energies):
            for s in range(spins):
                if model_kind == "subvec":
                    p = k + g2
                    mat = model(e, s + 1)(np.broadcast_to(p @ bmat, u.shape), (p + u) @ bmat, u @ bmat,
                                          np.all(u == 0, axis=-1))
                elif model_kind == "sub":
                    p = k + g2
                    pc = p @ bmat
                    puc = (p + u) @ bmat
                    uc = u @ bmat
                    mat = model(e, s + 1)(np.sum(pc ** 2, -1), np.sum(puc ** 2, -1),
                                          np.sum(uc ** 2, -1), np.all(u == 0, axis=-1))
                else:
                    m1 = (k + g1) @ bmat
                    m2 = (k + g2) @ bmat
                    m1, m2 = np.broadcast_arrays(m1, m2)
                    mat = model(e, s + 1)(m1, m2)
                mat = np.where(inside[:, None] & inside[None, :], mat, 0.0)
                data[ik, ie, s, 0] = mat.real.T
                data[ik, ie, s, 1] = mat.imag.T
    with h5py.File(name + ".hdf", "w") as h:
        h.create_dataset("2Dmat", data=data)
        h.create_dataset("energies", data=np.array([[e.real, e.imag] for e in energies]))
        h.create_dataset("gvecs", data=np.asarray(gvecs, dtype=np.int32))
        kd = h.create_dataset("kpts", data=np.asarray(kpts, dtype=float))
        ops = ops if ops is not None else [np.eye(2, dtype=int)]
        _att(kd, "nsym", len(ops), np.int32)
        kd.attrs.create("symop", np.array([o.T for o in ops], dtype=np.int32))
        kd.attrs.create("symtau", np.zeros((len(ops), 2)))
        kd.attrs.create("amat2d", np.asarray(amat, dtype=float).T)
        _att(kd, "rkmax", rkmax, float)
        _att(kd, "l_timerev", int(timerev), np.int32)
        _att(kd, "l_noco", int(spins == 4), np.int32)
        _att(kd, "l_soc", int(soc), np.int32)


def read_embfile(name):
    with h5py.File(name + ".hdf", "r") as h:
        d = h["2Dmat"][...]
        out = {
            "mat": d[:, :, :, 0] + 1j * d[:, :, :, 1],  # [k,e,s,g2,g1]
            "energies": h["energies"][...] @ np.array([1.0, 1j]),
            "kpts": h["kpts"][...],
            "gvecs": h["gvecs"][...],
            "attrs": {k: v for k, v in h["kpts"].attrs.items()},
        }
    out["mat"] = np.swapaxes(out["mat"], -1, -2)  # [k,e,s,g1,g2]
    return out


def read_klist(path):
    lines = open(path).read().split("\n")
    n = int(lines[0])
    return np.array([[float(x) for x in lines[i + 1].split()] for i in range(n)])


# --------------------------------------------------------------------------- reference table
class SigRef:
    """sigma_u(p) as gf_sigtab: exact on the mesh, else 4x4 Lagrange (fallback 2x2)."""

    def __init__(self, amat, nmesh, gvecs, rkmax, model, e, spin, interp="cubic"):
        self.b = recip(amat)
        self.n = np.asarray(nmesh)
        self.gset = {tuple(g) for g in gvecs}
        self.rk2 = rkmax ** 2 * (1.0 + 1e-10)
        self.f = model(e, spin)
        self.interp = interp

    def node(self, ip, u):
        c = np.mod(ip + self.n // 2, self.n) - self.n // 2  # k = c/N in [-1/2,1/2)
        g = (ip - c) // self.n
        if tuple(g) not in self.gset:
            return None
        p = ip / self.n
        pc, puc, uc = p @ self.b, (p + u) @ self.b, np.asarray(u) @ self.b
        p2, pu2 = pc @ pc, puc @ puc
        if p2 > self.rk2 or pu2 > self.rk2:
            return None
        return complex(self.f(p2, pu2, uc @ uc, np.all(np.asarray(u) == 0)))

    def __call__(self, p, u):
        x = np.asarray(p) * self.n
        if np.all(np.abs(x - np.rint(x)) < 1e-6):
            return self.node(np.rint(x).astype(int), u)
        if self.interp == "exact":
            return None
        i0 = np.floor(x).astype(int) - 1
        c = [[self.node(i0 + np.array([i, j]), u) for j in range(4)] for i in range(4)]
        t = x - i0
        if self.interp == "cubic" and all(v is not None for row in c for v in row):
            wx, wy = lagrange4(t[0]), lagrange4(t[1])
            return sum(wx[i] * wy[j] * c[i][j] for i in range(4) for j in range(4))
        if all(c[i][j] is not None for i in (1, 2) for j in (1, 2)):
            lx = [2.0 - t[0], t[0] - 1.0]
            ly = [2.0 - t[1], t[1] - 1.0]
            return sum(lx[i] * ly[j] * c[i + 1][j + 1] for i in range(2) for j in range(2))
        return None


class SpinTab:
    """2x2 amplitudes: from a noco table set (spins 11,22,12,21) with the spin frame rotated
    by theta, or the diagonal of a collinear set"""
    nsp = 2

    def __init__(self, tabs, theta=0.0):
        self.tabs = tabs
        self.th = np.radians(theta)

    def __call__(self, p, u):
        v = [t(p, u) for t in self.tabs]
        if any(x is None for x in v):
            return None
        if len(v) == 2:
            return np.diag(v).astype(complex)
        m = np.array([[v[0], v[2]], [v[3], v[1]]], dtype=complex)
        m[1, 0] *= np.exp(1j * self.th)
        m[0, 1] *= np.exp(-1j * self.th)
        return m


def _blk(tab, p, u):
    v = tab(p, u)
    return None if v is None else np.atleast_2d(np.asarray(v, dtype=complex))


def spinor(mat4, nsp):
    """full matrix (index (s-1)*n+G) from the file spins 11,22,12,21 ([..., spin, G, G'])"""
    if nsp == 1:
        return mat4[0]
    return np.block([[mat4[0], mat4[2]], [mat4[3], mat4[1]]])


def lagrange4(t):
    return [-(t - 1) * (t - 2) * (t - 3) / 6.0, t * (t - 2) * (t - 3) / 2.0,
            -t * (t - 1) * (t - 3) / 2.0, t * (t - 1) * (t - 2) / 6.0]


# --------------------------------------------------------------------------- transformer
class Geometry:
    def __init__(self, ao, asub, theta=0.0, rk=0.0):
        self.ao = np.asarray(ao)
        self.bo = recip(self.ao)
        self.bbo = self.bo @ self.bo.T
        self.ab = rotation(theta) @ np.asarray(asub)
        self.bb = recip(self.ab)
        self.m = self.bo @ np.linalg.inv(self.bb)
        self.minv = np.linalg.inv(self.m)
        self.rk = rk

    def home(self, kap, glist):
        return [n for n, g in enumerate(glist)
                if (kap + g) @ self.bbo @ (kap + g) <= self.rk * self.rk]

    def as_sub(self, w, eps=1e-6):
        """substrate vector of the ad-layer fractional difference w, or None"""
        x = np.asarray(w) @ self.m
        return np.rint(x).astype(int) if np.all(np.abs(x - np.rint(x)) < eps) else None


def sigma0(geo, tab, kap, glist, shift=(0.0, 0.0)):
    """Sigma0 in the global basis (zero outside the kappa sphere); spinor tables give the
    full matrix with the index (s-1)*ng+G"""
    nsp = getattr(tab, "nsp", 1)
    ng = len(glist)
    out = np.zeros((nsp * ng, nsp * ng), dtype=complex)
    hl = geo.home(kap, glist)
    for b in hl:
        p = (kap + glist[b]) @ geo.m
        for a in hl:
            u = geo.as_sub(glist[a] - glist[b])
            if u is None:
                continue
            v = _blk(tab, p, u)
            assert v is not None, ("missing amplitude", p, u)
            ph = np.exp(-1j * (u @ geo.bb) @ np.asarray(shift))
            for s1 in range(nsp):
                for s2 in range(nsp):
                    out[s1 * ng + a, s2 * ng + b] = ph * v[s1, s2]
    return out


def shells(geo, nshell):
    cand = [(i, j) for i in range(-12, 13) for j in range(-12, 13) if (i, j) != (0, 0)]
    lens = {c: np.linalg.norm(np.array(c) @ geo.bb) for c in cand}
    radii = []
    for _ in range(nshell):
        rest = [l for l in lens.values() if not radii or l > radii[-1] * (1 + 1e-6)]
        radii.append(min(rest))
    return [np.array(c) for c in cand if lens[c] <= radii[-1] * (1 + 1e-6)]


def find_k(kfile, q, eps=1e-5):
    for n, k in enumerate(kfile):
        d = k - q
        if np.all(np.abs(d - np.rint(d)) < eps):
            return n, np.rint(d).astype(int)
    raise KeyError(q)


def sigma_eff(geo, tab, kap, glist, q, up, e_up, spin, order):
    """Sigma0 + Umklapp correction for one kappa (home block, index (s-1)*nh+G).
    Returns the home indices and the matrix."""
    nsp = getattr(tab, "nsp", 1)
    hl = geo.home(kap, glist)
    nh = len(hl)
    ng = len(glist)
    idx = [s1 * ng + a for s1 in range(nsp) for a in hl]
    s0 = sigma0(geo, tab, kap, glist)[np.ix_(idx, idx)]
    gpos = {tuple(glist[a]): i for i, a in enumerate(hl)}

    def put(mat, ir, nr, ic, nc, v):
        for s1 in range(nsp):
            for s2 in range(nsp):
                mat[s1 * nr + ir, s2 * nc + ic] += v[s1, s2]

    def upmat(ik, sel):
        if nsp == 1:
            m = up["mat"][ik, e_up, min(spin, up["mat"].shape[2] - 1)]
            return m[np.ix_(sel, sel)]
        m4 = up["mat"][ik, e_up]
        if m4.shape[0] == 4:
            return spinor([x[np.ix_(sel, sel)] for x in m4], 2)
        return np.block([[m4[0][np.ix_(sel, sel)], np.zeros((len(sel), len(sel)))],
                         [np.zeros((len(sel), len(sel))), m4[min(1, m4.shape[0] - 1)][np.ix_(sel, sel)]]])

    blocks = []
    for u in q:
        qq = kap + u @ geo.minv
        ik, nshift = find_k(up["kpts"], qq)
        kf = up["kpts"][ik]
        sel = [m for m, g in enumerate(up["gvecs"]) if (kf + g) @ geo.bbo @ (kf + g) <= geo.rk ** 2]
        grel = [up["gvecs"][m] + nshift for m in sel]
        nb = len(sel)
        d = upmat(ik, sel).astype(complex)
        du0 = np.zeros((nsp * nb, nsp * nh), dtype=complex)
        d0u = np.zeros((nsp * nh, nsp * nb), dtype=complex)
        for n, g in enumerate(grel):
            put(d, n, nb, n, nb, _blk(tab, (qq + g) @ geo.m, np.zeros(2, int)))
            h = gpos.get(tuple(g))
            if h is None:
                continue
            p = (kap + g) @ geo.m
            put(du0, n, nb, h, nh, _blk(tab, p, u))
            put(d0u, h, nh, n, nb, _blk(tab, p + u, -u))
        blocks.append(dict(u=u, q=qq, g=grel, nb=nb, d=d, du0=du0, d0u=d0u))

    def duv(bi, bj):
        out = np.zeros((nsp * bi["nb"], nsp * bj["nb"]), dtype=complex)
        for n2, g2 in enumerate(bj["g"]):
            for n1, g1 in enumerate(bi["g"]):
                if np.all(g1 == g2):
                    put(out, n1, bi["nb"], n2, bj["nb"], _blk(tab, (bj["q"] + g2) @ geo.m, bi["u"] - bj["u"]))
        return out

    s2 = sum(-b["d0u"] @ np.linalg.solve(b["d"], b["du0"]) for b in blocks)
    if order == 2:
        return hl, s0 + s2
    if order == 3:
        s3 = np.zeros_like(s0)
        for i, bi in enumerate(blocks):
            for j, bj in enumerate(blocks):
                if i != j:
                    s3 += bi["d0u"] @ np.linalg.solve(bi["d"], duv(bi, bj)) @ np.linalg.solve(bj["d"], bj["du0"])
        return hl, s0 + s2 + s3
    sizes = [nsp * b["nb"] for b in blocks]
    off = np.concatenate([[0], np.cumsum(sizes)])
    big = np.zeros((off[-1], off[-1]), dtype=complex)
    rhs = np.zeros((off[-1], nsp * nh), dtype=complex)
    d0q = np.zeros((nsp * nh, off[-1]), dtype=complex)
    for i, bi in enumerate(blocks):
        sl = slice(off[i], off[i + 1])
        big[sl, sl] = bi["d"]
        rhs[sl] = bi["du0"]
        d0q[:, sl] = bi["d0u"]
        for j, bj in enumerate(blocks):
            if i != j:
                big[sl, off[j]:off[j + 1]] = duv(bi, bj)
    return hl, s0 - d0q @ np.linalg.solve(big, rhs)
