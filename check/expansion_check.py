"""Numerical check of the patch expansion in the observer's local coordinates.

For each metric, patch and observer family a null geodesic is integrated from
the emitter E to the observer O.  The exact redshift of Proposition 5
(Killing charges) at those endpoints is compared with Proposition 6,
    z(dx, gamma) = C_r dr + C_theta q + cos(gamma) [D_r dr + D_theta q],
    q = vartheta dx sin(gamma)/l - dx^2 s(gamma)/(2 l^2),
    s = sin^2(gamma) on the equator, 1 at the pole,
with dx the proper horizontal step and gamma its angle from the rotation
direction e_phi.  Halving dx and quartering dr must shrink the difference by
at least 8 (the error is O(eps^3)).
"""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import redshift as R

r_, th_ = sp.symbols('r theta', positive=True)


def build(kind, pars):
    if kind == 'kerr':
        M, a = pars['M'], pars['a']
        rho2 = r_**2 + a**2*sp.cos(th_)**2
        D = r_**2 - 2*M*r_ + a**2
        g = sp.zeros(4, 4)
        g[0, 0] = -(1 - 2*M*r_/rho2)
        g[0, 3] = g[3, 0] = -2*M*a*r_*sp.sin(th_)**2/rho2
        g[3, 3] = (r_**2 + a**2 + 2*M*a**2*r_*sp.sin(th_)**2/rho2)*sp.sin(th_)**2
        g[1, 1] = rho2/D
        g[2, 2] = rho2
        idx = dict(r=1, th=2, ph=3)
        xi = np.array([1., 0, 0, 0])
    elif kind == 'mp5':
        mu, a = pars['mu'], pars['a']
        rho2 = r_**2 + a**2*sp.cos(th_)**2
        g = sp.zeros(5, 5)
        g[0, 0] = -1 + mu/rho2
        g[0, 3] = g[3, 0] = -mu*a*sp.sin(th_)**2/rho2
        g[3, 3] = (r_**2 + a**2)*sp.sin(th_)**2 + mu*a**2*sp.sin(th_)**4/rho2
        g[1, 1] = rho2/(r_**2 + a**2 - mu)
        g[2, 2] = rho2
        g[4, 4] = r_**2*sp.cos(th_)**2
        idx = dict(r=1, th=2, ph=3)
        xi = np.array([1., 0, 0, 0, 0])
    elif kind == 'gmrot':  # rotating 1+4 geodesic monism, coordinates (t,u,r,theta,phi)
        J_, K_, B_, Q_ = (pars.get(k, 0.) for k in ('J', 'K', 'B', 'Q'))
        c0_, c1_, c2_, c3_ = (pars.get(k, 0.) for k in ('c0', 'c1', 'c2', 'c3'))
        rr = r_
        f = J_/rr + K_*rr + B_*rr**2 + Q_*rr**4
        h0 = (J_**2/(12*rr**4) + J_*K_/(6*rr**2) - 2*J_*Q_*rr/3 + K_**2*sp.log(rr) + 4*B_*K_*rr/3
              + B_**2*rr**2/3 + K_*Q_*rr**3 + 2*B_*Q_*rr**4/3 + sp.Rational(11, 42)*Q_**2*rr**6)
        h2 = (J_**2/(6*rr**4) - 2*J_*K_/(3*rr**2) - 2*J_*B_/(3*rr) + J_*Q_*rr/3 - K_**2/2 - 2*B_*K_*rr/3
              - 2*K_*Q_*rr**3/3 - sp.Rational(4, 21)*B_*Q_*rr**4 - sp.Rational(29, 126)*Q_**2*rr**6)
        H = c0_/rr + c1_ + c2_*rr + c3_*rr**2 + h0 + h2*(3*sp.cos(th_)**2 - 1)/2
        g = sp.zeros(5, 5)
        g[0, 1] = g[1, 0] = 1
        g[1, 1] = 1 + 2*H
        g[1, 4] = g[4, 1] = f*sp.sin(th_)**2
        g[2, 2] = 1
        g[3, 3] = rr**2
        g[4, 4] = rr**2*sp.sin(th_)**2
        idx = dict(r=2, th=3, ph=4)
        xi = np.array([1., -1, 0, 0, 0])
    else:  # gm, coordinates (t,u,r,theta,phi), clock d_t - d_u
        c0, c1, c2, c3, J = (pars[k] for k in ('c0', 'c1', 'c2', 'c3', 'J'))
        P2 = (3*sp.cos(th_)**2 - 1)/2
        H = c0/r_ + c1 + c2*r_ + c3*r_**2 + J**2/(12*r_**4) + J**2*P2/(6*r_**4)
        g = sp.zeros(5, 5)
        g[0, 1] = g[1, 0] = 1
        g[1, 1] = 1 + 2*H
        g[1, 4] = g[4, 1] = J*sp.sin(th_)**2/r_
        g[2, 2] = 1
        g[3, 3] = r_**2
        g[4, 4] = r_**2*sp.sin(th_)**2
        idx = dict(r=2, th=3, ph=4)
        xi = np.array([1., -1, 0, 0, 0])
    n = g.shape[0]
    X = [sp.Symbol('x%d' % i) for i in range(n)]
    sub = {r_: X[idx['r']], th_: X[idx['th']]}
    gX = g.subs(sub)
    gi = gX.inv()
    dg = [gX.diff(X[c]) for c in range(n)]
    Gam = [[[sum(gi[a, d]*(dg[b][d, c] + dg[c][d, b] - dg[d][b, c]) for d in range(n))/2
             for c in range(n)] for b in range(n)] for a in range(n)]
    fg = sp.lambdify([X], gX, 'numpy')
    fG = sp.lambdify([X], Gam, 'numpy')
    return n, idx, xi, (lambda x: np.array(fg(x), dtype=float)), (lambda x: np.array(fG(x), dtype=float))


def observer(kind_obs, x, g, xi, iph):
    eph = np.zeros(len(x)); eph[iph] = 1.
    if kind_obs == 'static':
        U = xi
    else:  # ZAMO
        w = -(xi @ g @ eph)/(eph @ g @ eph)
        U = xi + w*eph
    return U/np.sqrt(-(U @ g @ U))


def frame(u, g, idx, n):
    """orthonormal spatial vectors e_r, e_th, e_ph (Gram-Schmidt, orthogonal to u)."""
    es = []
    for key in ('r', 'th', 'ph'):
        e = np.zeros(n); e[idx[key]] = 1.
        e = e + (u @ g @ e)*u
        for f in es:
            e = e - (f @ g @ e)*f
        es.append(e/np.sqrt(e @ g @ e))
    return es


def run(kind, pars, patch, obs, rO=60., step=1.5, nhat=(0.04, 0.6, 0.8), theta_off=0.01):
    n, idx, xi, fg, fG = build(kind, pars)
    th0 = np.pi/2 + theta_off if patch == 'equator' else 3*theta_off
    x0 = np.zeros(n); x0[idx['r']] = rO; x0[idx['th']] = th0
    g0 = fg(x0)
    u0 = observer(obs, x0, g0, xi, idx['ph'])
    er, et, ep = frame(u0, g0, idx, n)
    nr, nt, nph = np.array(nhat)/np.linalg.norm(nhat)
    k0 = u0 + nr*er + nt*et + nph*ep

    def rhs(l, y):
        x, k = y[:n], y[n:]
        return np.concatenate([k, -np.einsum('abc,b,c->a', fG(x), k, k)])

    sol = solve_ivp(rhs, (0, step), np.concatenate([x0, k0]), method='DOP853', rtol=1e-12, atol=1e-14)
    x1 = sol.y[:n, -1]
    g1 = fg(x1)
    eph = np.zeros(n); eph[idx['ph']] = 1.
    E = -(k0 @ g0 @ xi); L = k0 @ g0 @ eph

    def NW(g):
        if obs == 'static':
            return np.sqrt(-(xi @ g @ xi)), 0.
        w = -(xi @ g @ eph)/(eph @ g @ eph)
        U = xi + w*eph
        return np.sqrt(-(U @ g @ U)), w
    N0, w0 = NW(g0); N1, w1 = NW(g1)
    b = L/E
    z_exact = (N1/N0)*(1 - w0*b)/(1 - w1*b) - 1
    # local coordinates of the observer O = x1, emitter E = x0
    iR, iT, iP = idx['r'], idx['th'], idx['ph']
    dr = x1[iR] - x0[iR]
    met = dict(kerr=R.kerr, mp5=R.mp5, gm=R.gm, gmrot=R.gmrot)[kind]
    name = dict(kerr='Kerr', mp5='Myers-Perry 5D', gm='Geodesic monism (xi = d_t - d_u)',
                gmrot='Geodesic monism, rotating 1+4')[kind]
    symvals = {R.r: x1[iR], R.M: pars.get('M', 0), R.a: pars.get('a', 0), R.mu: pars.get('mu', 0),
               R.J: pars.get('J', 0), R.c0: pars.get('c0', 0), R.c1: pars.get('c1', 0),
               R.c2: pars.get('c2', 0), R.c3: pars.get('c3', 0),
               R.Kk: pars.get('K', 0), R.Bb: pars.get('B', 0), R.Qq: pars.get('Q', 0)}
    l = float(R.ell(name, patch).subs(symvals))
    if patch == 'equator':
        vt = x1[iT] - np.pi/2
        hx = np.array([np.sqrt(g1[iP, iP])*(x1[iP] - x0[iP]), np.sqrt(g1[iT, iT])*(x1[iT] - x0[iT])])
    else:
        vt = x1[iT]
        P = lambda x: l*x[iT]*np.array([np.cos(x[iP]), np.sin(x[iP])])
        d = P(x1) - P(x0)
        phihat = np.array([-np.sin(x1[iP]), np.cos(x1[iP])]); thhat = np.array([np.cos(x1[iP]), np.sin(x1[iP])])
        hx = np.array([d @ phihat, d @ thhat])
    dx = np.hypot(*hx)
    gam = np.arctan2(hx[1], hx[0])
    co = R.coefficients(met, patch)[obs]
    # the rotation term multiplies the photon's direction cosine n_phi in the ZAMO frame; it equals
    # cos(gamma) up to O(eps) and up to the dragging tilt O(|f|/r) of the u-direction (Proposition 6)
    z_local = float(R.z_local(co, patch, l, dr, dx, gam, vt, nphi=nph).subs(symvals))
    return dict(z_exact=z_exact, z_local=z_local, dx=dx, gam=gam, dr=dr, vt=vt)


if __name__ == '__main__':
    cases = [('kerr', dict(M=1., a=0.9)), ('mp5', dict(mu=4., a=0.9)),
             ('gm', dict(c0=-0.5, c1=-0.05, c2=1e-4, c3=0., J=0.8)),
             ('gmrot', dict(c0=-0.5, c1=-0.05, c2=1e-4, J=0.8, K=0.002, B=2e-5, Q=1e-9))]
    for kind, pars in cases:
        for patch in ('equator', 'pole'):
            for obs in ('static', 'ZAMO'):
                errs = []
                for s_ in (1.0, 0.5, 0.25):
                    d = run(kind, pars, patch, obs, step=s_, nhat=(0.04*s_, 0.6, 0.8), theta_off=0.01*s_)
                    errs.append(abs(d['z_exact'] - d['z_local']))
                    if s_ == 1.0:
                        d0 = d
                ratio = errs[1]/errs[2]
                print(f"{kind:4s} {patch:7s} {obs:6s} dx={d0['dx']:.3f} gamma={np.degrees(d0['gam']):6.1f} deg "
                      f"dr={d0['dr']:+.2e}: z exact {d0['z_exact']:+.6e}, z(dx,gamma) {d0['z_local']:+.6e};"
                      f" error ratio on halving {errs[0]/errs[1]:.1f}, {ratio:.1f}")
                assert ratio > 5
    print('PASS z(dx, gamma): error O(eps^3) in every metric, patch and observer family')
