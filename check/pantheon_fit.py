"""Fit the directional pieces of the local Hubble law to Pantheon+.

Usage:  python3 check/pantheon_fit.py DATA_DIR
DATA_DIR must hold the two public Pantheon+ files
    Pantheon+SH0ES.dat
    Pantheon+SH0ES_STAT+SYS.cov
from https://github.com/PantheonPlusSH0ES/DataRelease (Pantheon+_Data/4_DISTANCES_AND_COVAR).

Model for the distance-modulus residual at fixed redshift, with n the unit
vector to the supernova:
    Delta mu(n) = m0 - (5/ln 10) * delta(n),   delta = dH/H,
    (a) cosine toward a fixed direction s:  delta = A_s (n.s)
    (b) free dipole:                        delta = d.n
    (c) dipole + traceless quadrupole:      delta = d.n + n.Q.n
Generalised least squares with the full STAT+SYS covariance restricted to
the selected supernovae.  The background is flat LCDM with H0 = 73.04,
Omega_m = 0.334 (the offset m0 absorbs any H0 change).

What the post needs from the fit (section "What the data fix"):
    static observers : delta(n) = -(c/H0) C_r (n.rhat)                 -> dipole only
    ZAMO observers   : delta(n) = (c/H0)[-C_r (n.rhat) + D_r (n.rhat)(n.phihat)]
                       -> dipole along -rhat plus Q = (c/H0) D_r (rhat phihat + phihat rhat)/2,
                          eigenvalues (+q, 0, -q), null eigenvector = theta-hat,
                          dipole orthogonal to the null eigenvector.
"""
import sys, os
import numpy as np
from scipy.integrate import quad

C_KMS = 299792.458
H0, OM = 73.04, 0.334
SHAPLEY = (201.98, -31.50)   # RA, Dec of A3558, core of the Shapley supercluster


def unit(ra, dec):
    ra, dec = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)], -1)


def mu_lcdm(z):
    f = lambda x: 1/np.sqrt(OM*(1+x)**3 + 1 - OM)
    dc = np.array([quad(f, 0, zz)[0] for zz in z])*C_KMS/H0
    return 5*np.log10((1+z)*dc) + 25


def load(ddir):
    d = np.genfromtxt(os.path.join(ddir, 'Pantheon+SH0ES.dat'), names=True, dtype=None, encoding=None)
    with open(os.path.join(ddir, 'Pantheon+SH0ES_STAT+SYS.cov')) as f:
        n = int(f.readline())
        cov = np.loadtxt(f).reshape(n, n)
    return d, cov


def gls(X, y, C):
    Ci = np.linalg.inv(C)
    F = X.T @ Ci @ X
    Fi = np.linalg.inv(F)
    b = Fi @ X.T @ Ci @ y
    r = y - X @ b
    return b, Fi, float(r @ Ci @ r)


def quad_basis(n):
    x, y, z = n.T
    # traceless symmetric Q = sum_k q_k B_k
    return np.stack([x*x - z*z, y*y - z*z, 2*x*y, 2*x*z, 2*y*z], 1)


def Qmatrix(q):
    a, b, c, d, e = q
    return np.array([[a, c, d], [c, b, e], [d, e, -a - b]])


def fit(d, cov, zcol, lo, hi):
    z = d[zcol]
    sel = (z > lo) & (z < hi) & (d['IS_CALIBRATOR'] == 0)
    idx = np.where(sel)[0]
    C = cov[np.ix_(idx, idx)]
    n = unit(d['RA'][idx], d['DEC'][idx])
    y = d['MU_SH0ES'][idx] - mu_lcdm(z[idx])
    k = -5/np.log(10)
    out = {'N': len(idx), 'unique': len(set(d['CID'][idx]))}
    s = unit(*SHAPLEY)
    _, _, chi20 = gls(np.ones((len(idx), 1)), y, C)
    b, F, chi2 = gls(np.stack([np.ones(len(idx)), n @ s], 1), y, C)
    out['shapley_mag'] = (b[1], np.sqrt(F[1, 1]))
    out['shapley_dHH'] = (b[1]/k, np.sqrt(F[1, 1])/abs(k))
    out['shapley_dchi2'] = chi20 - chi2
    b, F, chi2d = gls(np.column_stack([np.ones(len(idx)), n]), y, C)
    dip = b[1:4]/k
    out['dipole'] = (dip, np.sqrt(np.diag(F)[1:4])/abs(k), chi20 - chi2d)
    X = np.column_stack([np.ones(len(idx)), n, quad_basis(n)])
    b, F, chi2q = gls(X, y, C)
    out['dip_q'] = b[1:4]/k
    out['Q'] = Qmatrix(b[4:]/k)
    out['Q_err'] = np.sqrt(np.diag(F)[4:])/abs(k)
    out['dchi2_quad'] = chi2d - chi2q
    return out


def fibonacci_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2*i/n)
    th = np.pi*(1 + 5**0.5)*i
    return np.stack([np.cos(th)*np.sin(phi), np.sin(th)*np.sin(phi), np.cos(phi)], 1)


def radec(v):
    v = v/np.linalg.norm(v)
    return np.degrees(np.arctan2(v[1], v[0])) % 360, np.degrees(np.arcsin(v[2]))


def main(ddir):
    d, cov = load(ddir)
    for zcol in ('zHD', 'zCMB'):
        for lo, hi in ((0.015, 0.06), (0.0233, 0.15)):
            o = fit(d, cov, zcol, lo, hi)
            print(f'--- {zcol} in ({lo}, {hi}): N = {o["N"]} light curves, {o["unique"]} distinct SNe')
            a, e = o['shapley_mag']
            print(f'  cosine toward Shapley: {a:+.4f} +- {e:.4f} mag  ->  dH/H = {o["shapley_dHH"][0]:+.4f} +- {o["shapley_dHH"][1]:.4f}; Delta chi2 = {o["shapley_dchi2"]:.1f} (1 dof)')
            dip, de, dchi = o['dipole']
            ra, dec = radec(dip)
            print(f'  free dipole |d| = {np.linalg.norm(dip):.4f} toward RA {ra:.0f}, Dec {dec:+.0f};'
                  f' component errors {np.round(de, 4)}; Delta chi2 vs none = {dchi:.1f} (3 dof)')
            w, V = np.linalg.eigh(o['Q'])
            print(f'  quadrupole eigenvalues {np.round(w, 4)} (errors on the 5 components ~ {np.round(o["Q_err"], 4)});'
                  f' Delta chi2 = {o["dchi2_quad"]:.1f} (5 dof)')
            mid = V[:, 1]
            dq = o['dip_q']
            print(f'  |cos(dipole, middle eigenvector)| = {abs(dq @ mid)/np.linalg.norm(dq):.2f}'
                  f'  (model: 0 with a middle eigenvalue near 0)')


def comoving(z):
    f = lambda x: 1/np.sqrt(OM*(1+x)**3 + 1 - OM)
    return np.array([quad(f, 0, zz)[0] for zz in z])*C_KMS/H0


BASES = ['b1 = -d n_r', 'b2 = -d^2(1-n_r^2)/2', 'b3 = -d^2 n_th^2/2', 'b4 = -d n_th',
         'b5 = d n_r n_ph', 'b6 = d^2(1-n_r^2) n_ph/2', 'b7 = d n_th n_ph', 'b8 = d^2 n_th^2 n_ph/2']
MODELS = {
    # leading order in the horizontal regime (Corollary 7 of the post):
    # C_r, D_r and, on the equator, the (dx)^2 terms are dropped
    'equator, static': [3],              # z = (C_theta vt/l) dx sin(gamma)
    'equator, ZAMO': [3, 6],             # z = (vt/l) dx sin(gamma) (C_theta + D_theta cos(gamma))
    'pole': [3, 1],                      # z = (C_theta/l) dx (vt sin(gamma) - dx/(2l))
    # for comparison, Proposition 6 in full (nothing dropped)
    'full, equator, static': [3, 0, 1, 2],
    'full, equator, ZAMO': [3, 0, 1, 2, 4, 5, 6, 7],
}


def fit_local_model(d, cov, zcol, lo, hi, ngrid=2000, npsi=36):
    """Fit Proposition 6, z(dx, gamma; dr), written for a supernova at comoving
    distance d in direction n with components (n_r, n_th, n_ph) in the observer's
    local frame (rhat, thetahat, phihat):
        dx = d sqrt(n_th^2 + n_ph^2),  cos(gamma) = -n_ph/sqrt(..), sin(gamma) = -n_th/sqrt(..)
        dr = r_O - r_E = -d n_r - dx^2/(2r)
    The photon runs from E to O, along -n.  z_patch is a sum of the bases
        z = C_r b1 + K2 b2 + K3 b3 + K4 b4 + D_r b5 + K6 b6 + K7 b7 + K8 b8
    with K2 = C_r/r (+ C_theta/l^2 at the pole), K3 = C_theta/l^2 (equator),
    K4 = C_theta vartheta/l, K6 = D_r/r, K7 = D_theta vartheta/l, K8 = D_theta/l^2.
    The supernova's residual moves by Delta mu = -(5/ln10)[(1+z) D'/D - 1] z_patch (redshift shift and dimming),
    which is -(5/ln10)(c/H0) z_patch/d at low z.
    The frame orientation (3 angles) is scanned; the coefficients are linear (GLS).
    An isotropic d^2 term (the deceleration of the background, unknown here) is
    fitted as a nuisance, so b2 and b3 measure only their anisotropic parts."""
    z = d[zcol]
    sel = (z > lo) & (z < hi) & (d['IS_CALIBRATOR'] == 0)
    idx = np.where(sel)[0]
    C = cov[np.ix_(idx, idx)]
    L = np.linalg.cholesky(C)
    wh = lambda X: np.linalg.solve(L, X)
    nv = unit(d['RA'][idx], d['DEC'][idx])
    dist = comoving(z[idx])
    N = len(idx)
    k = -5/np.log(10)
    E = np.sqrt(OM*(1 + z[idx])**3 + 1 - OM)
    # Delta mu per unit z_patch: the shift to a larger observed redshift and the (1 + z_patch)^2 dimming
    # (reciprocity); at low z this is (5/ln 10)(c/H0)/d
    w = k*((1 + z[idx])*(C_KMS/H0/E)/dist - 1)
    y = wh(d['MU_SH0ES'][idx] - mu_lcdm(z[idx]))
    one = wh(np.ones(N))
    T1 = wh((w*dist)[:, None]*nv)                                             # d n_i
    T2 = wh(w*dist**2)                                                        # d^2
    T2b = wh(((w*dist**2)[:, None, None]*nv[:, :, None]*nv[:, None, :]).reshape(N, 9))
    T5 = wh(((w*dist)[:, None, None]*nv[:, :, None]*nv[:, None, :]).reshape(N, 9))
    T6 = wh((w*dist**2)[:, None]*nv)
    T6b = wh(((w*dist**2)[:, None, None, None]*nv[:, :, None, None]*nv[:, None, :, None]
              * nv[:, None, None, :]).reshape(N, 27))
    X0 = np.column_stack([one, T2])
    chi2_iso = float(np.sum((y - X0 @ np.linalg.lstsq(X0, y, rcond=None)[0])**2))
    best = {m: None for m in MODELS}
    for e in fibonacci_sphere(ngrid):
        t0 = np.array([0., 0, 1]) if abs(e[2]) < 0.9 else np.array([1., 0, 0])
        a1 = np.cross(e, t0); a1 /= np.linalg.norm(a1)
        a2 = np.cross(e, a1)
        for psi in np.linspace(0, 2*np.pi, npsi, endpoint=False):
            f = np.cos(psi)*a1 + np.sin(psi)*a2            # phihat
            t = np.cross(f, e)                               # thetahat, (r, th, ph) right handed
            ee, tt, ef, tf = [np.outer(p, q).ravel() for p, q in ((e, e), (t, t), (e, f), (t, f))]
            eef = np.einsum('i,j,k->ijk', e, e, f).ravel()
            ttf = np.einsum('i,j,k->ijk', t, t, f).ravel()
            B = np.stack([-T1 @ e, -(T2 - T2b @ ee)/2, -(T2b @ tt)/2, -T1 @ t,
                          T5 @ ef, (T6 @ f - T6b @ eef)/2, T5 @ tf, (T6b @ ttf)/2], 1)
            for m, cols in MODELS.items():
                X = np.column_stack([one, T2, B[:, cols]])
                b = np.linalg.lstsq(X, y, rcond=None)[0]
                chi2 = float(np.sum((y - X @ b)**2))
                if best[m] is None or chi2 < best[m][0]:
                    err = np.sqrt(np.diag(np.linalg.inv(X.T @ X)))
                    best[m] = (chi2, e, t, f, b[2:], err[2:], cols)
    return dict(N=N, chi2_iso=chi2_iso, best=best)


def physical(ddir):
    d, cov = load(ddir)
    MPC = 3.0857e22
    for zcol in ('zHD', 'zCMB'):
        o = fit_local_model(d, cov, zcol, 0.015, 0.06)
        print(f'=== z(dx, gamma) fit, {zcol}, 0.015 < z < 0.06, N = {o["N"]}, isotropic chi2 = {o["chi2_iso"]:.1f}')
        for m, (chi2, e, t, f, b, err, cols) in o['best'].items():
            if b[0] < 0:          # thetahat -> -thetahat (with phihat -> -phihat): b4 changes sign, b7 does not
                t, f = -t, -f
                b = np.array([-v if c == 3 else v for c, v in zip(cols, b)])
            ra_r, de_r = radec(-e)
            ra_f, de_f = radec(f)
            ra_t, de_t = radec(t)
            print(f'  [{m}] chi2 = {chi2:.1f} (Delta chi2 = {o["chi2_iso"] - chi2:.1f} vs isotropic with free d^2 term);'
                  f' thetahat toward RA {ra_t:.0f} Dec {de_t:+.0f}, phihat toward RA {ra_f:.0f} Dec {de_f:+.0f}')
            for c, v, s_ in zip(cols, b, err):
                unit_ = '/Mpc' if c in (0, 3, 4, 6) else '/Mpc^2'
                print(f'      {BASES[c]:26s} coefficient {v:+.3e} +- {s_:.1e} {unit_}  ({abs(v)/s_:.1f} sigma)')
            if 6 in cols:
                k7 = b[cols.index(6)]; k4 = b[0]
                e7 = err[cols.index(6)]
                print(f'      D_theta/C_theta = b7/b4 = {k7/k4:+.2f} +- {e7/abs(k4):.2f}   (Kerr, Myers-Perry: -2a/r)')


if __name__ == '__main__':
    ddir = sys.argv[1] if len(sys.argv) > 1 else 'data'
    main(ddir)
    physical(ddir)
