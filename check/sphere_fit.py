"""Fit the three-sphere rotating along phi (Proposition 14 of the post) to every supernova sample at hand.

Usage:  python3 check/sphere_fit.py DATA_DIR
DATA_DIR holds
    Pantheon+SH0ES.dat, Pantheon+SH0ES_STAT+SYS.cov
        from https://github.com/PantheonPlusSH0ES/DataRelease (Pantheon+_Data/4_DISTANCES_AND_COVAR)
    des/repo/4_DISTANCES_COVMAT/{DES-Dovekie_HD.csv, DES-Dovekie_Metadata.csv, STAT+SYS.npz}
        from https://github.com/des-science/DES-SN5YR (the DES-Dovekie recalibration, Popovic et al. 2026)

On the sky the leading redshift of a patch of the three-sphere, for a source at distance d in direction n, is
    z_patch = A1 (-d n.e1) + A2 (-d^2 (1 - (n.e_a)^2)/2),     A1 = C_theta vartheta/r,  A2 = C_theta/r^2,
with e1 perpendicular to e_a (Proposition 14: e1 the observer's offset from the patch circle, e_a along it).
A redshift z_patch at fixed distance moves the source to a larger observed redshift, (1 + z) -> (1 + z)(1 + z_patch),
and, being a real frequency shift, also dims it: by reciprocity d_L = (1 + z)^2 d_A, so d_L grows by (1 + z_patch)^2.
The Hubble residual at the observed redshift therefore moves by
    Delta mu = -(5/ln 10) [(1 + z) D'(z)/D(z) - 1] z_patch,       D the comoving distance,
which reduces to -(5/ln 10)(c/H0) z_patch/d at low z.  Isotropic nuisances: an offset and an isotropic d^2
term (so A2 measures only the anisotropic part), and for the wide windows the derivatives of mu with
respect to Omega_m and w (the background is not modelled here).  Full STAT+SYS covariances.
Models:  iso; dipole (A1 alone, free direction: the equator at leading order);
         sphere (A1 with e1 perpendicular to e_a, and A2; e_a scanned over the sphere).
The script also turns the CMB dipole and quadrupole and the spatial-curvature bound into limits on A1, A2 and r.
"""
import sys
import os
import numpy as np
import pandas as pd
from scipy.integrate import quad
import pantheon_fit as PF

C_KMS = 299792.458
H0, OM = 73.04, 0.334


def comoving(z, om=OM, w=-1.0):
    f = lambda x: 1/np.sqrt(om*(1 + x)**3 + (1 - om)*(1 + x)**(3*(1 + w)))
    return np.array([quad(f, 0, zz)[0] for zz in np.atleast_1d(z)])*C_KMS/H0


def mu_of(z, om=OM, w=-1.0):
    return 5*np.log10((1 + z)*comoving(z, om, w)) + 25


def dmu_dz(z):
    h = 1e-5*np.maximum(z, 1e-3)
    return (mu_of(z + h) - mu_of(z - h))/(2*h)


def load_pantheon(ddir, zcol):
    d, cov = PF.load(ddir)
    keep = d['IS_CALIBRATOR'] == 0
    idx = np.where(keep)[0]
    return dict(name=f'Pantheon+ {zcol}', z=d[zcol][idx], mu=d['MU_SH0ES'][idx], ra=d['RA'][idx], dec=d['DEC'][idx],
                cov=cov[np.ix_(idx, idx)], cid=np.array([str(c).strip() for c in d['CID'][idx]]))


def load_des(ddir):
    base = os.path.join(ddir, 'des', 'repo', '4_DISTANCES_COVMAT')
    rows = [l.split() for l in open(os.path.join(base, 'DES-Dovekie_HD.csv')) if l.startswith('SN:')]
    hd = pd.DataFrame([r[1:] for r in rows], columns='CID IDSURVEY zHD zHEL MU MUERR MUERR_VPEC MUERR_SYS PROBIA'.split())
    for c in hd.columns[1:]:
        hd[c] = hd[c].astype(float)
    meta = pd.read_csv(os.path.join(base, 'DES-Dovekie_Metadata.csv'), sep=r'\s+', comment='#')
    meta['CID'] = meta['CID'].astype(str)
    meta = meta.drop_duplicates('CID').set_index('CID')
    dd, _ = PF.load(ddir)
    pos = {}
    for c, ra, de in zip(dd['CID'], dd['RA'], dd['DEC']):
        pos.setdefault(str(c).strip(), (ra, de))
    ra = np.full(len(hd), np.nan); dec = np.full(len(hd), np.nan)
    for i, (c, s) in enumerate(zip(hd.CID, hd.IDSURVEY)):
        if s == 10 and c in meta.index and meta.loc[c, 'HOST_RA'] > -900:
            ra[i], dec[i] = meta.loc[c, 'HOST_RA'], meta.loc[c, 'HOST_DEC']
        elif c in pos:
            ra[i], dec[i] = pos[c]
    npz = np.load(os.path.join(base, 'STAT+SYS.npz'))
    n = int(npz[npz.files[0]][0])
    icov = np.zeros((n, n))
    icov[np.triu_indices(n)] = npz[npz.files[1]]
    il = np.tril_indices(n, -1)
    icov[il] = icov.T[il]
    cov = np.linalg.inv(icov)
    good = np.isfinite(ra)
    idx = np.where(good)[0]
    return dict(name='DES-Dovekie zHD', z=hd.zHD.values[idx], mu=hd.MU.values[idx], ra=ra[idx], dec=dec[idx],
                cov=cov[np.ix_(idx, idx)], cid=hd.CID.values[idx], dropped=int((~good).sum()))


def fit(S, lo, hi, wide=False, naxis=3000):
    sel = (S['z'] > lo) & (S['z'] < hi)
    idx = np.where(sel)[0]
    z = S['z'][idx]
    C = S['cov'][np.ix_(idx, idx)]
    L = np.linalg.cholesky(C)
    wh = lambda X: np.linalg.solve(L, X)
    n = PF.unit(S['ra'][idx], S['dec'][idx])
    d = comoving(z)
    E = np.sqrt(OM*(1 + z)**3 + 1 - OM)
    w = -(5/np.log(10))*((1 + z)*(C_KMS/H0/E)/d - 1)     # Delta mu per unit z_patch (shift and dimming)
    y = wh(S['mu'][idx] - mu_of(z))
    nuis = [np.ones(len(idx)), w*d**2]
    if wide:
        e = 1e-3
        nuis += [(mu_of(z, OM + e) - mu_of(z, OM - e))/(2*e), (mu_of(z, OM, -1 + e) - mu_of(z, OM, -1 - e))/(2*e)]
    Xn = wh(np.column_stack(nuis))
    def solve(X):
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = y - X @ b
        return b, float(r @ r), np.linalg.inv(X.T @ X)
    _, chi_iso, _ = solve(Xn)
    D = wh((w*d)[:, None]*(-n))                 # columns: -d n_i  (times w)
    bd, chi_dip, Fd = solve(np.column_stack([Xn, D]))
    k = Xn.shape[1]
    dip = bd[k:k + 3]
    J = dip/np.linalg.norm(dip)
    sig_dip = float(np.sqrt(J @ Fd[k:k + 3, k:k + 3] @ J))
    T2 = (w*d**2)[:, None, None]*n[:, :, None]*n[:, None, :]
    T2 = wh(T2.reshape(len(idx), 9))
    best = None
    axes = PF.fibonacci_sphere(2*naxis)
    axes = axes[axes[:, 2] >= 0]
    chis = []
    for ea in axes:
        t0 = np.array([0., 0, 1]) if abs(ea[2]) < 0.9 else np.array([1., 0, 0])
        a1 = np.cross(ea, t0); a1 /= np.linalg.norm(a1)
        a2 = np.cross(ea, a1)
        Bq = 0.5*(T2 @ np.outer(ea, ea).ravel())                   # A2 basis: -d^2(1 - n_a^2)/2 = (+d^2 n_a^2 - d^2)/2; the
        X = np.column_stack([Xn, D @ a1, D @ a2, Bq])              # isotropic -d^2/2 part is absorbed by the nuisance
        b, chi, F = solve(X)
        chis.append(chi)
        if best is None or chi < best[0]:
            v = b[k]*a1 + b[k + 1]*a2
            A1 = np.linalg.norm(v)
            Jv = np.array([b[k], b[k + 1]])/max(A1, 1e-30)
            sA1 = float(np.sqrt(Jv @ F[k:k + 2, k:k + 2] @ Jv))
            best = (chi, ea, v/max(A1, 1e-30), A1, sA1, b[k + 2], float(np.sqrt(F[k + 2, k + 2])))
    chis = np.array(chis)
    # axes within Delta chi2 < 2.3 of the best (two angles): the 68% region of e_a
    near = axes[chis - chis.min() < 2.3]
    spread = np.degrees(np.arccos(np.clip(np.abs(near @ best[1]), -1, 1))).max()
    return dict(N=len(idx), zmed=float(np.median(z)), dmax=float(d.max()), chi_iso=chi_iso, chi_dip=chi_dip,
                dip=np.linalg.norm(dip), sdip=sig_dip, dipdir=dip/np.linalg.norm(dip), sphere=best, spread=spread)


def fit_shapley(S, lo, hi, wide=False, npsi=720):
    """the dipole fixed toward Shapley: z_patch = A_s d (n.s) + A2(-d^2 (1 - (n.e_a)^2)/2), e_a perpendicular to s.
    A_s > 0 means more redshift toward Shapley (-e1 = s).  On the equator patch e_a is the rotation direction."""
    sel = (S['z'] > lo) & (S['z'] < hi)
    idx = np.where(sel)[0]
    z = S['z'][idx]
    C = S['cov'][np.ix_(idx, idx)]
    L = np.linalg.cholesky(C)
    wh = lambda X: np.linalg.solve(L, X)
    n = PF.unit(S['ra'][idx], S['dec'][idx])
    d = comoving(z)
    E = np.sqrt(OM*(1 + z)**3 + 1 - OM)
    w = -(5/np.log(10))*((1 + z)*(C_KMS/H0/E)/d - 1)
    y = wh(S['mu'][idx] - mu_of(z))
    nuis = [np.ones(len(idx)), w*d**2]
    if wide:
        e = 1e-3
        nuis += [(mu_of(z, OM + e) - mu_of(z, OM - e))/(2*e), (mu_of(z, OM, -1 + e) - mu_of(z, OM, -1 - e))/(2*e)]
    Xn = wh(np.column_stack(nuis))
    s = PF.unit(*PF.SHAPLEY)
    Ds = wh(w*d*(n @ s))
    def solve(X):
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = y - X @ b
        return b, float(r @ r), np.linalg.inv(X.T @ X)
    _, chi_iso, _ = solve(Xn)
    bs, chi_s, Fs = solve(np.column_stack([Xn, Ds]))
    k = Xn.shape[1]
    t0 = np.array([0., 0, 1])
    a1 = np.cross(s, t0); a1 /= np.linalg.norm(a1)
    a2 = np.cross(s, a1)
    T2 = wh(((w*d**2)[:, None, None]*n[:, :, None]*n[:, None, :]).reshape(len(idx), 9))
    rows = []
    for psi in np.linspace(0, np.pi, npsi, endpoint=False):           # an axis: psi and psi + pi are the same
        ea = np.cos(psi)*a1 + np.sin(psi)*a2
        b, chi, F = solve(np.column_stack([Xn, Ds, 0.5*(T2 @ np.outer(ea, ea).ravel())]))
        rows.append((chi, psi, ea, b[k], np.sqrt(F[k, k]), b[k + 1], np.sqrt(F[k + 1, k + 1])))
    rows.sort(key=lambda t: t[0])
    best = rows[0]
    good = [t for t in rows if t[0] - best[0] < 1.0]                  # one angle: Delta chi2 < 1
    span = max(np.degrees(np.arccos(np.clip(abs(t[2] @ best[2]), -1, 1))) for t in good)
    return dict(N=len(idx), chi_iso=chi_iso, As=(bs[k], np.sqrt(Fs[k, k])), chi_s=chi_s, best=best, span=span)


def report_shapley(name, lo, hi, o):
    chi, psi, ea, As, sAs, A2, sA2 = o['best']
    ra, de = PF.radec(ea)
    if de < 0:
        ra, de = (ra + 180) % 360, -de
    print(f'  {name}, {lo} < z < {hi}: N = {o["N"]}; dipole toward Shapley alone A_s = ({o["As"][0]*1e6:+.2f} +- {o["As"][1]*1e6:.2f})e-6 /Mpc'
          f' (Delta chi2 = {o["chi_iso"] - o["chi_s"]:.1f}, 1 dof);'
          f' with the d^2 term: A_s = ({As*1e6:+.2f} +- {sAs*1e6:.2f})e-6, A2 = ({A2*1e9:+.1f} +- {sA2*1e9:.1f})e-9 /Mpc^2,'
          f' axis e_a toward RA {ra:.0f} Dec {de:+.0f} (+- {o["span"]:.0f} deg), Delta chi2 over the Shapley dipole {o["chi_s"] - chi:.1f} (2 dof)')
    return dict(name=name, lo=lo, hi=hi, As=As, sAs=sAs, A2=A2, sA2=sA2, ra=ra, de=de, span=o['span'],
                dchi_s=o['chi_iso'] - o['chi_s'], dchi_ax=o['chi_s'] - chi, As0=o['As'])


def report(name, lo, hi, o):
    chi, ea, e1, A1, sA1, A2, sA2 = o['sphere']
    print(f'  {name}, {lo} < z < {hi}: N = {o["N"]}, d_max = {o["dmax"]:.0f} Mpc')
    ra, de = PF.radec(o['dipdir'])
    print(f'    dipole (equator, leading order): A1 = ({o["dip"]*1e6:.2f} +- {o["sdip"]*1e6:.2f})e-6 /Mpc'
          f' along -n toward RA {ra:.0f} Dec {de:+.0f}; Delta chi2 = {o["chi_iso"] - o["chi_dip"]:.1f} (3 dof)')
    ra1, de1 = PF.radec(e1); raa, dea = PF.radec(ea)
    print(f'    three-sphere: A1 = ({A1*1e6:.2f} +- {sA1*1e6:.2f})e-6 /Mpc (e1 toward RA {ra1:.0f} Dec {de1:+.0f}),'
          f' A2 = ({A2*1e9:+.3f} +- {sA2*1e9:.3f})e-9 /Mpc^2 ({A2/sA2:+.1f} sigma);'
          f' e_a along RA {raa:.0f} Dec {dea:+.0f} (68% within {o["spread"]:.0f} deg);'
          f' Delta chi2 = {o["chi_iso"] - chi:.1f} (5 dof)')
    return dict(name=name, lo=lo, hi=hi, A1=A1, sA1=sA1, A2=A2, sA2=sA2, dip=o['dip'], sdip=o['sdip'],
                dchi_dip=o['chi_iso'] - o['chi_dip'], dchi_sph=o['chi_iso'] - chi, N=o['N'], dmax=o['dmax'],
                dipdir=o['dipdir'], ea=ea)


def bounds(results):
    """CMB and curvature, turned into limits on A1, A2 and r (see the post)."""
    d_lss = comoving(np.array([1089.9]))[0]
    d_sn = comoving(np.array([2.26]))[0]
    dT_dip = 1.234e-3            # CMB dipole Delta T/T (3362 muK / 2.7255 K)
    v_kin = 450./C_KMS           # allowance for a kinematic dipole of up to 450 km/s with any direction
    # rms Delta T/T of the CMB quadrupole: D_2 = 225.9 (+533.1) muK^2 (Planck 2018); the upper value, 759 muK^2,
    # gives C_2 = 2 pi D_2/6 and an rms 5 C_2/(4 pi) over the sky
    q_rms = np.sqrt(5*(2*np.pi*759.0/6)/(4*np.pi))/2.7255e6
    pat = np.sqrt(4/45)          # rms over the sphere of n_a^2 - 1/3
    print(f'  CMB: last scattering at d = {d_lss:.0f} Mpc; the most distant Pantheon+ supernova at d = {d_sn:.0f} Mpc')
    for lab, dd in (('d_LSS', d_lss), ('d(z = 2.26), conservative', d_sn)):
        A1max = (dT_dip + v_kin)/dd
        A2max = 2*q_rms/(pat*dd**2)
        print(f'    if the CMB comes from the patch beyond {lab}: |A1| < {A1max:.1e} /Mpc, |A2| < {A2max:.1e} /Mpc^2')
    cH = C_KMS/67.4
    for lab, om_min in (('Planck 2018 + BAO, Omega_K = 0.001 +- 0.002, 2 sigma', -0.003),
                        ('DESI DR2 + CMB, Omega_K = 0.0023 +- 0.0011, 3 sigma', 0.0023 - 3*0.0011)):
        rmin = cH/np.sqrt(-om_min)
        print(f'    a closed space of radius r has Omega_K = -(c/H0 r)^2: {lab} -> r > {rmin/1e3:.0f} Gpc')
    return d_lss, d_sn


def gal2eq(l, b):
    """galactic (l, b) -> unit vector in ICRS"""
    from astropy.coordinates import SkyCoord
    c = SkyCoord(l=l, b=b, unit='deg', frame='galactic').icrs
    return PF.unit(c.ra.deg, c.dec.deg)


def compare(res):
    """the dipole of the patch is linear in d: a bulk-flow survey of a sphere of radius R sees v = (3/4) c A1 R"""
    cmb = gal2eq(264.021, 48.253)        # Planck 2018: 369.82 km/s toward (l, b) = (264.021, 48.253)
    bf = gal2eq(297., -6.)               # Watkins et al. 2023 (CF4): 387 +- 28 km/s within 150/h Mpc toward (297, -6)
    print('  directions (RA, Dec): CMB dipole %.0f %+.0f, CF4 bulk flow %.0f %+.0f' % (*PF.radec(cmb), *PF.radec(bf)))
    R = 150/(H0/100)
    for r in res:
        if r['hi'] <= 0.15:
            v = 0.75*C_KMS*r['dip']*R
            sv = 0.75*C_KMS*r['sdip']*R
            up = -r['dipdir']                  # z_patch is largest toward -e1
            ang = lambda a, b: np.degrees(np.arccos(np.clip(a @ b, -1, 1)))
            print(f'    {r["name"]}, {r["lo"]}-{r["hi"]}: bulk flow within {R:.0f} Mpc from the fitted dipole,'
                  f' {v:.0f} +- {sv:.0f} km/s (CF4: 387 +- 28); its direction is {ang(up, bf):.0f} deg from the CF4 flow'
                  f' and {ang(up, cmb):.0f} deg from the CMB dipole')
    print(f'    with the CMB limit |A1| < 2.2e-7 /Mpc: at most {0.75*C_KMS*2.2e-7*R:.0f} km/s')


def d_lss_planck():
    """comoving distance to last scattering, Planck 2018 (H0 = 67.4, Omega_m = 0.315, radiation 9.1e-5)"""
    h, om, orad = 67.4, 0.315, 9.1e-5
    f = lambda x: 1/np.sqrt(orad*(1 + x)**4 + om*(1 + x)**3 + 1 - om - orad)
    return quad(f, 0, 1089.9, limit=200)[0]*C_KMS/h


def example_C():
    """C_theta of the ZAMO on the stable orbit of the worked example of the post (units of r_s)"""
    J, K, B, Q, c0, c1, c2, c3 = 0.5, -0.2, 0.005, -1e-4, -1.0, 0.0, -0.005, 0.02
    F = lambda r: J/r**2 + K + B*r**2 + Q*r**4
    def h0(r):
        L = np.log(r)
        return (J*J/(12*r**6) + J*K/(6*r**4) - K*K*L/(4*r**2) - J*Q*L + B*K*L + B*B*r**2/4
                + K*Q*r**2*L/2 - K*Q*r**2/8 + B*Q*r**4/2 + 19*Q*Q*r**6/96)
    def h2(r):
        L = np.log(r)
        return (J*K*L/(3*r**4) + K*K*L/(4*r**2) + K*K/(16*r**2) + B*J/(2*r**2) + B*K/2
                + K*Q*r**2*L/6 + B*Q*r**4/8 + 27*Q*Q*r**6/160)
    Hr = lambda r: c0/r**2 + c1 + c2*r**2 + c3*np.log(r)
    We = lambda r: -2*(Hr(r) + h0(r) + h2(r)) + F(r)**2/r**2
    from scipy.optimize import minimize_scalar
    rs = minimize_scalar(We, bounds=(2, 6), method='bounded', options=dict(xatol=1e-10)).x
    N2 = 1 + We(rs)
    return rs, (4*h2(rs) - F(rs)**2/rs**2)/N2, N2


def unknowns():
    dl = d_lss_planck()
    cH = C_KMS/67.4
    dT_dip, v_kin = 1.234e-3, 450./C_KMS
    q_rms = np.sqrt(5*(2*np.pi*759.0/6)/(4*np.pi))/2.7255e6
    kq = 2*q_rms/np.sqrt(4/45)
    print(f'  last scattering at d_LSS = {dl/1e3:.1f} Gpc (Planck 2018 background)')
    print(f'  CMB quadrupole: |C_theta| < {kq:.2e} (r/d_LSS)^2;  CMB dipole: |C_theta| vartheta < {dT_dip + v_kin:.2e} r/d_LSS')
    for lab, om in (('Planck 2018 + BAO (2 sigma)', -0.003), ('DESI DR2 + CMB (3 sigma)', 0.0023 - 3*0.0011)):
        r = cH/np.sqrt(-om)
        print(f'    {lab}: r > {r/1e3:.0f} Gpc -> at that r, |C_theta| < {kq*(r/dl)**2:.1e} and |C_theta| vartheta < {(dT_dip + v_kin)*r/dl:.3f}')
    rs, Ct, N2 = example_C()
    rmin = dl*np.sqrt(Ct/kq)
    print(f'  the worked example: stable orbit at r_s = {rs:.4f} in its own length unit, ZAMO C_theta = (4h2 - F^2/r^2)/N^2 = {Ct:.3e}'
          f' (N^2 = {N2:.4f}); it passes the CMB for r_s > {rmin/1e3:.0f} Gpc')


def main(ddir):
    res = []
    print('=== the three-sphere pattern, fitted sample by sample')
    for zcol in ('zHD', 'zCMB'):
        S = load_pantheon(ddir, zcol)
        for lo, hi, wide in ((0.015, 0.06, False), (0.0233, 0.15, False), (0.0233, 2.3, True)):
            res.append(report(S['name'], lo, hi, fit(S, lo, hi, wide)))
    S = load_des(ddir)
    print(f'  (DES-Dovekie: {S["dropped"]} low-z supernovae without a position dropped)')
    for lo, hi, wide in ((0.02, 0.1, False), (0.02, 1.2, True)):
        res.append(report(S['name'], lo, hi, fit(S, lo, hi, wide)))
    print('=== the dipole fixed toward the Shapley supercluster; the axis e_a scanned perpendicular to it')
    for zcol in ('zHD', 'zCMB'):
        S = load_pantheon(ddir, zcol)
        for lo, hi in ((0.015, 0.06), (0.0233, 0.15)):
            report_shapley(S['name'], lo, hi, fit_shapley(S, lo, hi))
    S = load_des(ddir)
    report_shapley(S['name'], 0.02, 0.1, fit_shapley(S, 0.02, 0.1))
    print('=== bounds from the CMB and from spatial curvature')
    bounds(res)
    print('=== other data: bulk flow and dipole directions')
    compare(res)
    print('=== what the data leave of the unknowns')
    unknowns()
    return res


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'data')
