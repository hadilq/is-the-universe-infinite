"""Orbits of the three-sphere rotating along phi (geodesic monism, 1 + 1 + 4).

The transverse space is flat R^4, so the geodesics are integrated in Cartesian coordinates
(t, u, x1, x2, x3, x4), where the rotation axis is not a coordinate singularity:
    ds^2 = 2 dt du + (1 + 2H) du^2 + 2 (F(r)/r^2)(x1 dx2 - x2 dx1) du + dx.dx,
    mu = (x1^2 + x2^2)/r^2,   H = H_rad + h0(r) + h2(r)(2 mu - 1).
Along every geodesic p = u' and E_u = k_u are constant, and g(k, k) = -kappa gives
    |x'|^2 - 2 p^2 H(x) = p^2 - kappa - 2 p E_u,
so the motion in R^4 is that of a charge in the potential -2p^2 H and the magnetic field p dA
(which does no work).  For zero angular momentum along phi the reduced potential is p^2 W,
    W(r, mu) = -2H + F^2 mu/r^2 = N^2 - 1,   linear in mu.
Checks:
  1. the ZAMO on the rotation plane (mu = 1) at the minimum r_s of W_e(r) = W(r, 1) is a
     geodesic, and the orbit is stable: W_e''(r_s) > 0 and dW/dmu < 0 there;
  2. the conserved |x'|^2 - 2p^2 H holds along integrated geodesics;
  3. geodesics launched from the orbit in random directions, null and timelike with
     increasing speed, either stay in a band of r or leak; the leak goes through the
     rotation axis, where W -> -infinity at both ends:
         W ~ J^2 (mu - 1/6)/r^6 near the centre,   W ~ Q^2 r^6 (13 mu/40 - 7/120) far out.
"""
import sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

PARS = dict(J=0.5, K=-0.2, B=0.005, Q=-1e-4, c0=-1.0, c1=0.0, c2=-0.005, c3=0.02)
xs = sp.symbols('x1 x2 x3 x4', real=True)


def radial(P, r):
    J, K, B, Q = P['J'], P['K'], P['B'], P['Q']
    L = sp.log(r)
    F = J/r**2 + K + B*r**2 + Q*r**4
    h0 = (J**2/(12*r**6) + J*K/(6*r**4) - K**2*L/(4*r**2) - J*Q*L + B*K*L + B**2*r**2/4
          + K*Q*r**2*L/2 - K*Q*r**2/8 + B*Q*r**4/2 + sp.Rational(19, 96)*Q**2*r**6)
    h2 = (J*K*L/(3*r**4) + K**2*L/(4*r**2) + K**2/(16*r**2) + B*J/(2*r**2) + B*K/2
          + K*Q*r**2*L/6 + B*Q*r**4/8 + sp.Rational(27, 160)*Q**2*r**6)
    Hr = P['c0']/r**2 + P['c1'] + P['c2']*r**2 + P['c3']*L
    return F, h0, h2, Hr


def build(P):
    R = sp.sqrt(sum(x**2 for x in xs))
    F, h0, h2, Hr = radial(P, R)
    mu = (xs[0]**2 + xs[1]**2)/R**2
    H = Hr + h0 + h2*(2*mu - 1)
    G = F/R**2
    A = [-G*xs[1], G*xs[0], 0, 0]
    X = [sp.Symbol('t'), sp.Symbol('u')] + list(xs)
    g = sp.zeros(6, 6)
    g[0, 1] = g[1, 0] = 1
    g[1, 1] = 1 + 2*H
    for i in range(4):
        g[1, 2 + i] = g[2 + i, 1] = A[i]
        g[2 + i, 2 + i] = 1
    # inverse: g^tu = 1, g^tt = -(1+2H) + |A|^2, g^ti = -A_i, g^ij = delta
    gi = sp.zeros(6, 6)
    gi[0, 1] = gi[1, 0] = 1
    gi[0, 0] = -(1 + 2*H) + sum(a**2 for a in A)
    for i in range(4):
        gi[0, 2 + i] = gi[2 + i, 0] = -A[i]
        gi[2 + i, 2 + i] = 1
    assert sp.simplify((g*gi - sp.eye(6)).subs({xs[0]: 0.3, xs[1]: -1.1, xs[2]: 0.7, xs[3]: 2.0})) == sp.zeros(6, 6)
    dg = [g.diff(X[c]) for c in range(6)]
    Gam = [[[sum(gi[a, d]*(dg[b][d, c] + dg[c][d, b] - dg[d][b, c]) for d in range(6))/2
             for c in range(6)] for b in range(6)] for a in range(6)]
    fg = sp.lambdify([X], g, 'numpy')
    fG = sp.lambdify([X], Gam, 'numpy')
    fH = sp.lambdify([X], H, 'numpy')
    r = sp.Symbol('r', positive=True)
    Fr, h0r, h2r, Hrr = radial(P, r)
    m = sp.Symbol('m')
    W = -2*(Hrr + h0r + h2r*(2*m - 1)) + Fr**2*m/r**2
    fW = sp.lambdify((r, m), W, 'numpy')
    fWr = sp.lambdify((r, m), sp.diff(W, r), 'numpy')
    fWrr = sp.lambdify((r, m), sp.diff(W, r, 2), 'numpy')
    fWm = sp.lambdify(r, sp.diff(W, m), 'numpy')
    fF = sp.lambdify(r, Fr, 'numpy')
    fHrm = sp.lambdify((r, m), Hrr + h0r + h2r*(2*m - 1), 'numpy')
    return (lambda x: np.array(fg(x), dtype=float)), (lambda x: np.array(fG(x), dtype=float)), fH, fW, fWr, fWrr, fWm, fF, fHrm


MU = np.concatenate([[0.], np.geomspace(1e-8, 0.5, 400), 1 - np.geomspace(1e-8, 0.5, 400)[::-1][1:], [1.]])


def vstar(fF, fHrm, r, p, kphi, kc):
    """V_*(r) = min over mu of V(r, mu) = -2p^2 H + (k_phi - p F mu)^2/(r^2 mu) + k_c^2/(r^2 (1 - mu))"""
    r = np.atleast_1d(r)[:, None]
    m = MU[None, :]
    with np.errstate(divide='ignore', invalid='ignore'):
        t1 = np.where(m > 0, (kphi - p*fF(r)*m)**2/(r**2*m), np.where(kphi == 0, 0., np.inf))
        t2 = np.where(m < 1, kc**2/(r**2*(1 - m)), np.where(kc == 0, 0., np.inf))
    return np.min(-2*p**2*fHrm(r, m) + t1 + t2, axis=1)


def band(fF, fHrm, r0, eps, p, kphi, kc, rlo=1e-3, rhi=1e3):
    """the component of {r : V_*(r) <= eps} that contains r0: (r_min, r_max); 0 or inf if it is open"""
    f = lambda x: float(vstar(fF, fHrm, x, p, kphi, kc)[0] - eps)
    down = np.geomspace(r0, rlo, 3000); up = np.geomspace(r0, rhi, 3000)
    tol = 1e-9*max(1., abs(eps))          # a launch tangent to the orbit starts at a turning point, V_*(r0) = eps
    vd = vstar(fF, fHrm, down, p, kphi, kc) - eps - tol
    vu = vstar(fF, fHrm, up, p, kphi, kc) - eps - tol
    vd[0] = vu[0] = min(vd[0], 0.)
    i = int(np.argmax(vd > 0)) if np.any(vd > 0) else None
    j = int(np.argmax(vu > 0)) if np.any(vu > 0) else None
    g = lambda x: f(x) - tol
    r1 = brentq(g, down[i], down[i - 1], xtol=1e-10) if i is not None else 0.
    r2 = brentq(g, up[j - 1], up[j], xtol=1e-10) if j is not None else np.inf
    return r1, r2


def integrate(fG, x0, k0, lmax):
    def rhs(l, y):
        x, k = y[:6], y[6:]
        return np.concatenate([k, -np.einsum('abc,b,c->a', fG(x), k, k)])

    def centre(l, y):
        return np.linalg.norm(y[2:6]) - 0.05
    centre.terminal = True

    def far(l, y):
        return np.linalg.norm(y[2:6]) - 150.
    far.terminal = True
    return solve_ivp(rhs, (0, lmax), np.concatenate([x0, k0]), method='DOP853', rtol=1e-11, atol=1e-13,
                     events=(centre, far), max_step=0.25)


def zamo(fg, x):
    g = fg(x)
    xi = np.array([1., -1, 0, 0, 0, 0])
    eph = np.array([0, 0, -x[3], x[2], 0, 0])          # d_phi = x1 d_2 - x2 d_1
    w = -(xi @ g @ eph)/(eph @ g @ eph)
    U = xi + w*eph
    return g, U/np.sqrt(-(U @ g @ U))


def main():
    ok = True
    fg, fG, fH, fW, fWr, fWrr, fWm, fF, fHrm = build(PARS)
    rr = np.geomspace(0.5, 30, 4000)
    We = fW(rr, 1.0)
    i = [j for j in range(1, len(rr) - 1) if We[j] < We[j - 1] and We[j] < We[j + 1]][0]
    rs = brentq(lambda x: fWr(x, 1.0), rr[i - 1], rr[i + 1], xtol=1e-14)
    print(f'parameters {PARS}')
    print(f'rotation plane: W_e = W(r, 1) has a minimum at r_s = {rs:.6f}; N^2 = 1 + W = {1 + fW(rs, 1.0):.4f},'
          f' W_e\'\'(r_s) = {fWrr(rs, 1.0):.3e}, dW/dmu = F^2/r^2 - 4h2 = {fWm(rs):.3e} (< 0: stable off the plane)')
    ok &= bool(fWrr(rs, 1.0) > 0 and fWm(rs) < 0)
    # the barrier: lowest W on any path from the orbit to r -> 0 or r -> infinity (mu in [0, 1])
    rg = np.geomspace(2e-3, 2e3, 2500)
    Wmin_mu = np.minimum(fW(rg, 0.0), fW(rg, 1.0))     # W is linear in mu: its minimum over the sphere is at mu = 0 or 1
    j0 = int(np.argmin(abs(rg - rs)))
    inner = Wmin_mu[:j0].max() if j0 > 0 else np.inf    # pass needed to get in: the lowest crest on the way to r -> 0
    outer = Wmin_mu[j0:].max()
    # r -> 0 and r -> infinity along mu = 0
    print(f'W(r, 0) on the axis at the ends: W(2e-3, 0) = {fW(2e-3, 0.0):.3g}, W(2e3, 0) = {fW(2e3, 0.0):.3g};'
          f' in the plane: W(2e-3, 1) = {fW(2e-3, 1.0):.3g}, W(2e3, 1) = {fW(2e3, 1.0):.3g}')
    print(f'barrier from the orbit: inward {inner - fW(rs, 1.0):.3f}, outward {outer - fW(rs, 1.0):.3f} (in units of W)')
    # 1. the ZAMO stays
    x0 = np.array([0, 0, rs, 0, 0, 0.])
    g0, U0 = zamo(fg, x0)
    sol = integrate(fG, x0, U0, 300.)
    drift = np.max(np.abs(np.linalg.norm(sol.y[2:6], axis=0) - rs))
    mudrift = np.max(1 - (sol.y[2]**2 + sol.y[3]**2)/np.sum(sol.y[2:6]**2, axis=0))
    print(f'ZAMO released at r_s in the rotation plane: max |r - r_s| = {drift:.1e}, max (1 - mu) = {mudrift:.1e} over proper time 300')
    ok &= bool(drift < 1e-6 and mudrift < 1e-6)
    # 2-3. random geodesics from the orbit
    rng = np.random.default_rng(0)
    worst = 0.
    for label, speeds in (('null', [1.0]), ('timelike', [0.3, 0.6, 0.9])):
        for v in speeds:
            leaks, rmin, rmax, mumin, n = 0, np.inf, 0., 1., 16
            inside, fill, plo, phi_ = True, [], np.inf, 0.
            for k in range(n):
                ang = rng.uniform(0, 2*np.pi)
                x0 = np.array([0, 0, rs*np.cos(ang), rs*np.sin(ang), 0, 0.])
                g0, U0 = zamo(fg, x0)
                es = []
                for j in (2, 3, 4, 5, 1):
                    e = np.zeros(6); e[j] = 1.
                    e = e + (U0 @ g0 @ e)*U0
                    for f in es:
                        e = e - (f @ g0 @ e)*f
                    es.append(e/np.sqrt(e @ g0 @ e))
                nv = rng.normal(size=5); nv /= np.linalg.norm(nv)
                if label == 'null':
                    k0 = U0 + sum(a*e for a, e in zip(nv, es))
                else:
                    gam = 1/np.sqrt(1 - v**2)
                    k0 = gam*(U0 + v*sum(a*e for a, e in zip(nv, es)))
                p = k0[1]
                E0 = k0[2:] @ k0[2:] - 2*p**2*fH(x0)
                sol = integrate(fG, x0, k0, 100.)
                R = np.linalg.norm(sol.y[2:6], axis=0)
                mu = (sol.y[2]**2 + sol.y[3]**2)/R**2
                Ered = np.sum(sol.y[8:12]**2, axis=0) - 2*sol.y[7]**2*fH(sol.y[:6])
                if sol.status == 0:
                    worst = max(worst, np.max(np.abs(Ered - E0))/max(1., abs(E0)))
                leaks += sol.status == 1
                rmin, rmax, mumin = min(rmin, R.min()), max(rmax, R.max()), min(mumin, mu.min())
                # the band predicted from the constants and the conserved p, k_phi, k_c and eps
                kphi = (g0 @ k0) @ np.array([0, 0, -x0[3], x0[2], 0, 0])
                kc = x0[4]*k0[5] - x0[5]*k0[4]
                r1, r2 = band(fF, fHrm, rs, E0, p, kphi, kc)
                plo, phi_ = min(plo, r1), max(phi_, r2)
                inside &= bool(R.min() >= r1*(1 - 1e-6) and R.max() <= r2*(1 + 1e-6))
                if np.isfinite(r2) and r1 > 0:
                    fill.append((R.max() - R.min())/(r2 - r1))
            print(f'  {label:8s} v = {v:.1f}: {n} geodesics from the orbit, {leaks} leak (reach r < 0.05 or r > 150);'
                  f' r in [{rmin:.3f}, {rmax:.3f}], smallest mu {mumin:.3f};'
                  f' predicted bands within [{plo:.3f}, {phi_:.3f}], every orbit inside its band: {inside},'
                  f' median fill {np.median(fill):.2f}')
            ok &= inside
            if label == 'null' or v <= 0.6:
                ok &= bool(leaks == 0)
    print(f'|x\'|^2 - 2p^2 H conserved to {worst:.1e} along the bounded geodesics')
    ok &= bool(worst < 1e-6)
    # 3b. radial kicks in the rotation plane: W_e(r_pm) = W_s + v^2 N_s^2 exactly, r_pm ~ r_s +- v N_s sqrt(2/W_e'')
    Ws, Ns2, Wpp = fW(rs, 1.0), 1 + fW(rs, 1.0), fWrr(rs, 1.0)
    for v in (0.1, 0.3, 0.6, 0.9, 1.0):
        x0 = np.array([0, 0, rs, 0, 0, 0.])
        g0, U0 = zamo(fg, x0)
        er = np.array([0, 0, 1., 0, 0, 0]); er = er + (U0 @ g0 @ er)*U0; er /= np.sqrt(er @ g0 @ er)
        k0 = U0 + er if v == 1.0 else (U0 + v*er)/np.sqrt(1 - v**2)
        sol = integrate(fG, x0, k0, 120.)
        R = np.linalg.norm(sol.y[2:6], axis=0)
        target = Ws + v**2*Ns2
        rm = brentq(lambda x: fW(x, 1.0) - target, 0.3, rs)
        rp = brentq(lambda x: fW(x, 1.0) - target, rs, 200.)
        small = v*np.sqrt(Ns2)*np.sqrt(2/Wpp)
        err = max(abs(R.min() - rm)/rm, abs(R.max() - rp)/rp)
        print(f'  radial kick v = {v:.1f}: integrated r in [{R.min():.5f}, {R.max():.5f}], W_e(r) = W_s + v^2 N_s^2 gives'
              f' [{rm:.5f}, {rp:.5f}] (relative difference {err:.1e}); small-v formula r_s -+ {small:.4f}')
        ok &= bool(err < 1e-4)
    # 3c. any direction in the ZAMO frame (no kick along u), speed v, direction cosine n_phi:
    #     eps/p^2 = W_s + v^2 N_s^2 + 2 v n_phi N_s F_s/r_s,  k_phi/p = -v n_phi N_s r_s,  k_c = 0
    Fs = fF(rs); Ns = np.sqrt(Ns2)
    for v in (0.3, 0.6, 0.9, 1.0):
        lo, hi = np.inf, 0.
        for nphi in np.linspace(-1, 1, 41):
            eps = Ws + v**2*Ns2 + 2*v*nphi*Ns*Fs/rs
            r1, r2 = band(fF, fHrm, rs, eps, 1.0, -v*nphi*Ns*rs, 0.)
            lo, hi = min(lo, r1), max(hi, r2)
        print(f'  every geodesic launched from the orbit at speed v = {v:.1f} (any direction, no kick along u) stays in'
              f' [{lo:.3f}, {hi:.3f}]')
        ok &= bool(lo > 0 and np.isfinite(hi))
    # 4. through the axis: released at rest on the rotation axis just inside the crest of W(r, 0), it falls in
    rg2 = np.geomspace(0.05, rs, 4000)
    Wax = fW(rg2, 0.0)
    ic = int(np.argmax(Wax))
    r_in = brentq(lambda x: fW(x, 0.0), rg2[0], rg2[ic])          # W = 0, so N^2 = 1 there
    x0 = np.array([0, 0, 0, 0, r_in, 0.])
    g0 = fg(x0)
    U0 = np.array([1., -1, 0, 0, 0, 0]); U0 = U0/np.sqrt(-(U0 @ g0 @ U0))
    sol = integrate(fG, x0, U0, 50.)
    fell = sol.status == 1 and np.linalg.norm(sol.y[2:6, -1]) < 0.06
    print(f'  the crest of W on the rotation axis is at r = {rg2[ic]:.3f} (W = {Wax[ic]:.2f}); an observer released at rest on'
          f' the axis at r = {r_in:.3f}, inside it, reaches the centre: {fell}')
    ok &= bool(fell)
    print('PASS three-sphere orbits' if ok else 'FAIL three-sphere orbits')
    return ok


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
