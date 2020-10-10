"""
Code for calculating a PFSS extrapolation.
"""
import numpy as np

import pfsspy


HAS_NUMBA = False
try:
    import numba
    HAS_NUMBA = True
except Exception:
    pass


def _eigh(A):
    return np.linalg.eigh(A)


def _compute_r_term(m, k, ns, Q, brt, lam, ffm, nr, ffp, psi, psir):
    for l in range(ns):
        # Ignore the l=0 and m=0 term; for a globally divergence free field
        # this term is zero anyway, but numerically it may be small which
        # causes numerical issues when solving for c, d
        if l == 0 and m == 0:
            continue
        # - sum (c_{lm} + d_{lm}) * lam_{l}
        # lam[l] is small so this blows up
        cdlm = np.dot(Q[:, l], brt[:, m]) / lam[l]
        # - ratio c_{lm}/d_{lm} [numerically safer this way up]
        ratio = (ffm[l]**(nr - 1) - ffm[l]**nr) / (ffp[l]**nr - ffp[l]**(nr - 1))
        dlm = cdlm / (1.0 + ratio)
        clm = ratio * dlm
        psir[:, l] = clm * ffp[l]**k + dlm * ffm[l]**k

    # - compute entry for this m in psit = Sum_l c_{lm}Q_{lm}**j
    psi[:, :, m] = np.dot(psir, Q.T)
    return psi, psir


def _als_alp(nr, nphi, Fs, psi, Fp, als, alp):
    for j in range(nr + 1):
        for i in range(nphi + 1):
            als[i, :, j] = Fs * (psi[j, :, ((i - 1) % nphi)] - psi[j, :, ((i) % nphi)])
        for i in range(nphi):
            alp[i, 1:-1, j] = Fp[1:-1] * (psi[j, 1:, i] - psi[j, :-1, i])
    return als, alp


def _A_diag(A, ns, Vg, Uc, mu, m):
    for j in range(ns):
        A[j, j] = Vg[j] + Vg[j + 1] + Uc[j] * mu[m]
    return A


if HAS_NUMBA:
    _eigh = numba.jit(nopython=True)(_eigh)
    _compute_r_term = numba.jit(nopython=True)(_compute_r_term)
    _als_alp = numba.jit(nopython=True)(_als_alp)
    _A_diag = numba.jit(nopython=True)(_A_diag)


def pfss(input):
    r"""
    Compute PFSS model.

    Extrapolates a 3D PFSS using an eigenfunction method in :math:`r,s,p`
    coordinates, on the dumfric grid
    (equally spaced in :math:`\rho = \ln(r/r_{sun})`,
    :math:`s= \cos(\theta)`, and :math:`p=\phi`).

    Parameters
    ----------
    input : :class:`Input`
        Input parameters.

    Returns
    -------
    out : :class:`Output`

    Notes
    -----
    In order to avoid numerical issues, the monopole term (which should be zero
    for a physical magnetic field anyway) is explicitly excluded from the
    solution.

    The output should have zero current to machine precision,
    when computed with the DuMFriC staggered discretization.
    """
    br0 = input.br
    nr = input.grid.nr
    ns = input.grid.ns
    nphi = input.grid.nphi

    # Coordinates:
    ds = input.grid.ds
    dp = input.grid.dp
    dr = input.grid.dr

    rg = input.grid.rg
    rc = input.grid.rc

    sg = input.grid.sg
    sc = input.grid.sc

    k = np.linspace(0, nr, nr + 1)

    Fp = sg * 0  # Lp/Ls on p-ribs
    Fp[1:-1] = np.sqrt(1 - sg[1:-1]**2) / (np.arcsin(sc[1:]) - np.arcsin(sc[:-1])) * dp
    Vg = Fp / ds / dp
    Fs = (np.arcsin(sg[1:]) - np.arcsin(sg[:-1])) / np.sqrt(1 - sc**2) / dp  # Ls/Lp on s-ribs
    Uc = Fs / ds / dp

    # FFT in phi of photospheric distribution at each latitude:
    brt = np.fft.rfft(br0, axis=1)
    brt = brt.astype(np.complex128)

    # Prepare tridiagonal matrix:
    # - create off-diagonal part of the matrix:
    A = np.zeros((ns, ns))
    for j in range(ns - 1):
        A[j, j + 1] = -Vg[j + 1]
        A[j + 1, j] = A[j, j + 1]
    # - term required for m-dependent part of matrix:
    mu = np.fft.fftfreq(nphi)
    mu = 4 * np.sin(np.pi * mu)**2
    # - initialise:
    psir = np.zeros((nr + 1, ns), dtype='complex')
    psi = np.zeros((nr + 1, ns, nphi), dtype='complex')
    e1 = np.exp(dr)
    fact = np.sinh(dr) * (e1 - 1)

    # Loop over azimuthal modes (positive m):
    for m in range(nphi // 2 + 1):
        # - set diagonal terms of matrix:
        A = _A_diag(A, ns, Vg, Uc, mu, m)

        # - compute eigenvectors Q_{lm} and eigenvalues lam_{lm}:
        #   (note that A is symmetric so use special solver)
        lam, Q = _eigh(A)
        Q = Q.astype(np.complex128)
        # - solve quadratic:
        Flm = 0.5 * (1 + e1 + lam * fact)
        ffp = Flm + np.sqrt(Flm**2 - e1)
        ffm = e1 / ffp

        # - compute radial term for each l (for this m):
        psi, psir = _compute_r_term(m, k, ns, Q, brt, lam, ffm, nr, ffp, psi, psir)

        if (m > 0):
            psi[:, :, nphi - m] = np.conj(psi[:, :, m])

    # Past this point only psi, Fs, Fp are needed
    # Compute psi by inverse fft:
    psi = np.real(np.fft.ifft(psi, axis=2))

    # Hence compute vector potential [note index order, for netcdf]:
    # (note that alr is zero by definition)
    alr = np.zeros((nphi + 1, ns + 1, nr))
    als = np.zeros((nphi + 1, ns, nr + 1))
    alp = np.zeros((nphi, ns + 1, nr + 1))

    als, alp = _als_alp(nr, nphi, Fs, psi, Fp, als, alp)

    r = np.exp(rg)
    th = np.arccos(sg)
    ph = np.linspace(0, 2 * np.pi, nphi + 1)

    return pfsspy.Output(alr, als, alp, input.grid, input.map)
