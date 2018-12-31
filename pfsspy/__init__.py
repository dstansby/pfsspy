import numpy as n
import scipy.linalg as la
import pfsspy.output


def pfss(br0, nr, ns, np, rss, filename='', output='a', testQ=False):
    """
    Compute PFSS model.

    Extrapolates a 3D PFSS using an eigenfunction method in r,s,p coordinates,
    on the dumfric grid
    (equally spaced in rho=ln(r/rsun), s=cos(theta0), and p=phi).

    The output should have zero current to machine precision,
    when computed with the DuMFriC staggered discretization.

    Output depends on the flag 'output':

    - output='none': as it says
    - output='a': ar*Lr, as*Ls, ap*Lp on cell edges.
    - output='bc': br, bs, bp on the centres of the cell faces.
    - output='bg': br, bs, bp (weighted) averaged to grid points.

    Paramters
    ---------
    br0

    nr

    ns

    np

    rss

    filename : str
        Output filename.

    ouput : str
        String from ``('none', 'a', 'bc', 'bg')``.

    testQ : bool
        If ``True``, compare the discrete eigenfunctions Qj_{lm} to
        Plm(cos(th)).
    """

    # Coordinates:
    ds = 2.0 / ns
    dp = 2 * n.pi / np
    dr = n.log(rss) / nr
    rg = n.linspace(0, n.log(rss), nr + 1)
    rc = n.linspace(0.5 * dr, n.log(rss) - 0.5 * dr, nr)
    sg = n.linspace(-1, 1, ns + 1)
    sc = n.linspace(-1 + 0.5 * ds, 1 - 0.5 * ds, ns)

    k = n.linspace(0, nr, nr + 1)

    Fp = sg * 0  # Lp/Ls on p-ribs
    Fp[1:-1] = n.sqrt(1 - sg[1:-1]**2) / (n.arcsin(sc[1:]) - n.arcsin(sc[:-1])) * dp
    Vg = Fp / ds / dp
    Fs = (n.arcsin(sg[1:]) - n.arcsin(sg[:-1]))/n.sqrt(1 - sc**2) / dp  # Ls/Lp on s-ribs
    Uc = Fs / ds / dp

    # FFT in phi of photospheric distribution at each latitude:
    brt = n.fft.rfft(br0, axis=1)

    # Prepare tridiagonal matrix:
    # - create off-diagonal part of the matrix:
    A = n.zeros((ns, ns))
    for j in range(ns-1):
        A[j, j+1] = -Vg[j+1]
        A[j+1, j] = A[j, j+1]
    # - term required for m-dependent part of matrix:
    mu = n.fft.fftfreq(np)
    mu = 4 * n.sin(n.pi * mu)**2
    # - initialise:
    psir = n.zeros((nr + 1, ns), dtype='complex')
    psi = n.zeros((nr + 1, ns, np), dtype='complex')
    e1 = n.exp(dr)
    fact = n.sinh(dr) * (e1 - 1)

    if (testQ):
        import scipy.special as sp
        import matplotlib.pyplot as plt
        plt.figure()

    # Loop over azimuthal modes (positive m):
    for m in range(np // 2 + 1):
        # - set diagonal terms of matrix:
        for j in range(ns):
            A[j,j] = Vg[j] + Vg[j+1] + Uc[j]*mu[m]
        # - compute eigenvectors Q_{lm} and eigenvalues lam_{lm}:
        #   (note that A is symmetric so use special solver)
        lam, Q = la.eigh(A)
        # - solve quadratic:
        Flm = 0.5*(1 + e1 + lam*fact)
        ffp = Flm + n.sqrt(Flm**2 - e1)
        ffm = e1/ffp
        # - compute radial term for each l (for this m):
        for l in range(ns):
            # - sum c_{lm} + d_{lm}
            cdlm = n.dot(Q[:,l], brt[:,m])/lam[l]
            # - ratio c_{lm}/d_{lm} [numerically safer this way up]
            ratio = (ffm[l]**(nr-1) - ffm[l]**nr)/(ffp[l]**nr - ffp[l]**(nr-1))
            dlm = cdlm/(1.0 + ratio)
            clm = ratio*dlm
            psir[:,l] = clm*ffp[l]**k + dlm*ffm[l]**k
        # - compute entry for this m in psit = Sum_l c_{lm}Q_{lm}**j
        psi[:,:,m] = n.dot(psir, Q.T)
        if (m > 0):
            psi[:,:,np-m] = n.conj(psi[:,:,m])

        if (testQ & (m==6)):
            isrt = n.argsort(lam, axis=0) # sort eigenvalues
            lam = lam[isrt]
            istat = n.indices((ns, ns))
            Q = Q[istat[0],isrt]
            plt.clf()
            for l in range(5):
                plm = sp.lpmv(m, m+l, sc)
                Ql = Q[:,l]*Q[1,l]/n.abs(Q[1,l]) # normalise and match sign
                plt.plot(sc, Ql/n.max(n.abs(Ql)), 'ko')
                plm = plm*plm[1]/n.abs(plm[1])
                plt.plot(sc, plm/n.max(n.abs(plm)), label='l=%i'%l)
                plt.xlabel(r'$\cos(\theta)$')
            plt.title('m = %i' % m)
            plt.legend()
            plt.savefig('Q.png',bbox_inches='tight')
            plt.show()

    del(psir, mu, A)

    # Compute psi by inverse fft:
    psi = n.real(n.fft.ifft(psi, axis=2))

    # Hence compute vector potential [note index order, for netcdf]:
    alr = n.zeros((np+1, ns+1, nr))
    als = n.zeros((np+1, ns, nr+1))
    alp = n.zeros((np, ns+1, nr+1))

    for j in range(nr+1):
        for i in range(np+1):
            als[i,:,j] = Fs*(psi[j,:,((i-1)%np)] - psi[j,:,((i)%np)])
        for i in range(np):
            alp[i,1:-1,j] = Fp[1:-1]*(psi[j,1:,i] - psi[j,:-1,i])

    # Output to netcdf file:
    r = n.exp(rg)
    th = n.arccos(sg)
    ph = n.linspace(0, 2*n.pi, np+1)

    if (output=='none'):
         return alr, als, alp

    if (output=='a'):
        pfsspy.output.a(filename, r, th, ph, alr, als, alp)

    if ((output=='bc') | (output=='bg')):
        rc = n.linspace(-0.5*dr, n.log(rss)+0.5*dr, nr+2)
        rrc = n.exp(rc)
        thc = n.zeros(ns+2) - 1
        thc[1:-1] = n.arccos(sc)
        pc = n.linspace(-0.5*dp, 2*n.pi+0.5*dp, np+2)

        # Required face normals:
        dnp = n.zeros((ns+2,2))
        dns = n.zeros((ns+1,2))
        dnr = n.zeros(ns+2)
        for k in range(2):
            for j in range(1,ns+1):
                dnp[j,k] = rrc[k]*n.sqrt(1 - sc[j-1]**2)*dp
            dnp[0,k] = dnp[1,k]
            dnp[-1,k] = dnp[-2,k]
            for j in range(1,ns):
                dns[j,k] = rrc[k]*(n.arcsin(sc[j]) - n.arcsin(sc[j-1]))
            dns[0,k] = dns[1,k]
            dns[-1,k] = dns[-2,k]
        for j in range(ns+2):
            dnr[j] = rrc[0]*(n.exp(dr) - 1)
        dnr[0] = -dnr[0]
        dnr[-1] = -dnr[-1]

        # Required area factors:
        Sbr = n.zeros((ns+2,nr+1))
        for k in range(nr+1):
            Sbr[1:-1,k] = n.exp(2*rg[k])*ds*dp
            Sbr[0,k] = Sbr[1,k]
            Sbr[-1,k] = Sbr[-2,k]
        Sbs = n.zeros((ns+1,nr+2))
        for k in range(nr+2):
            for j in range(1,ns):
                Sbs[j,k] = 0.5*n.exp(2*rc[k] - dr)*dp*(n.exp(2*dr)-1)*n.sqrt(1 - sg[j]**2)
            Sbs[0,k] = Sbs[1,k]
            Sbs[-1,k] = Sbs[-2,k]
        Sbp = n.zeros((ns+2,nr+2))
        for k in range(nr+2):
            for j in range(1,ns+1):
                Sbp[j,k] = 0.5*n.exp(2*rc[k] - dr)*(n.exp(2*dr) - 1)*(n.arcsin(sg[j]) - n.arcsin(sg[j-1]))
            Sbp[0,k] = Sbp[1,k]
            Sbp[-1,k] = Sbp[-2,k]

        # Compute br*Sbr, bs*Sbs, bp*Sbp at cell centres by Stokes theorem:
        br = n.zeros((np+2,ns+2,nr+1))
        bs = n.zeros((np+2,ns+1,nr+2))
        bp = n.zeros((np+1,ns+2,nr+2))
        br[1:-1,1:-1,:] = als[1:,:,:] - als[:-1,:,:] + alp[:,:-1,:] - alp[:,1:,:]
        bs[1:-1,:,1:-1] = alp[:,:,1:] - alp[:,:,:-1]
        bp[:,1:-1,1:-1] = als[:,:,:-1] - als[:,:,1:]

        del(alr, als, alp)

        # Fill ghost values with boundary conditions:
        # - zero-gradient at outer boundary:
        bs[1:-1,:,-1] = 2*bs[1:-1,:,-2] - bs[1:-1,:,-3]
        bp[:,1:-1,-1] = 2*bp[:,1:-1,-2] - bp[:,1:-1,-3]
        # - periodic in phi:
        bs[0,:,:] = bs[-2,:,:]
        bs[-1,:,:] = bs[1,:,:]
        br[0,:,:] = br[-2,:,:]
        br[-1,:,:] = br[1,:,:]
        # js = jp = 0 at photosphere:
        for i in range(np+1):
            bp[i,:,0] = Sbp[:,0]/dnp[:,0]*(bp[i,:,1]*dnp[:,1]/Sbp[:,1] + br[i,:,0]*dnr[:]/Sbr[:,0] - br[i+1,:,0]*dnr[:]/Sbr[:,0])
        for i in range(np+2):
            bs[i,:,0] = Sbs[:,0]/dns[:,0]*(bs[i,:,1]*dns[:,1]/Sbs[:,1] + br[i,:-1,0]*dnr[:-1]/Sbr[:-1,0] - br[i,1:,0]*dnr[1:]/Sbr[1:,0])
        # - polar boundaries as in dumfric:
        for i in range(np+2):
            i1 = (i + np//2) % np
            br[i,-1,:] = br[i1,-2,:]
            br[i,0,:] = br[i1,1,:]
            bs[i,-1,:] = 0.5*(bs[i,-2,:] - bs[i1,-2,:])
            bs[i,0,:] = 0.5*(bs[i,1,:] - bs[i1,1,:])
        for i in range(np+1):
            i1 = (i + np//2) % np
            bp[i,-1,:] = -bp[i1,-2,:]
            bp[i,0,:] = -bp[i1,1,:]

        if (output == 'bc'):
            # Remove area factors:
            for i in range(np + 2):
                br[i, :, :] = br[i, :, :] / Sbr
                bs[i, :, :] = bs[i, :, :] / Sbs
            for i in range(np + 1):
                bp[i, :, :] = bp[i, :, :] / Sbp

            pfsspy.output.bc(filename, r, th, ph, rrc, thc, pc, br, bs, bp)

        if (output == 'bg'):
            # Weighted average to grid points:
            brg = br[:-1, :-1, :] + br[1:,:-1,:] + br[1:,1:,:] + br[:-1,1:,:]
            bsg = bs[:-1, :, :-1] + bs[1:,:,:-1] + bs[1:,:,1:] + bs[:-1,:,1:]
            bpg = bp[:, :-1, :-1] + bp[:,1:,:-1] + bp[:,1:,1:] + bp[:,:-1,1:]
            for i in range(np+1):
                brg[i,:,:] /= 2*(Sbr[:-1,:] + Sbr[1:,:])
                bsg[i,:,:] /= 2*(Sbs[:,:-1] + Sbs[:,1:])
            for i in range(np+1):
                bpg[i,:,:] /= Sbp[:-1,:-1] + Sbp[1:,:-1] + Sbp[1:,1:] + Sbp[:-1,1:]

            pfsspy.output.bg(filename, r, th, ph, brg, bsg, bpg)
