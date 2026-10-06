"""No leak into r: the three-sphere of geodesic monism rotating equally in both planes, 1 + 1 + 4
(the closing subsection of the 1 + 1 + 4 section of the post).

Metric, Hopf coordinates (t, u, r, theta, phi, psi), theta in [0, pi/2]:
    ds^2 = 2 dt du + (1+2H) du^2 + 2 F (sin^2 dphi + cos^2 dpsi) du
           + dr^2 + r^2 (dtheta^2 + sin^2 dphi^2 + cos^2 dpsi^2),
    F = J/r^2 + K + B r^2 + Q r^4,
    H = c0/r^2 + c1 + c2 r^2 + c3 ln r + J^2/(6r^6) + JK/(3r^4) - K^2 ln r/(2r^2)
        + KQ r^2 ln r + BQ r^4 + 19Q^2 r^6/48.
The sphere rotates along chi = d_phi + d_psi (both planes at the same rate), |chi| = r.
Along every geodesic, with p = u', P = k(chi) = k_phi + k_psi
and the Casimir l^2 of the three-sphere,
    r'^2 + V(r) = const,   V = -2 p^2 H + (1/r^2) [ l^2 - P^2 + (P - p F)^2 ].
This script checks that V is conserved along integrated geodesics (which checks
the radial reduction), that the ZAMO at the minimum r_s of U = -2H + F^2/r^2 is a
geodesic that stays there, and that random null and timelike geodesics launched
in any direction from r_s stay in the band V(r) <= const: no geodesic reaches the
centre or escapes, for J != 0 and Q != 0.
"""
import sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar, brentq

r_, th_ = sp.symbols('r theta', positive=True)
PARS = dict(J=2.0, K=0.2, B=0.004, Q=2e-4, c0=-0.3, c1=0.0, c2=0.0, c3=0.01)


def funcs(P):
    J, K, B, Q = P['J'], P['K'], P['B'], P['Q']
    f = J/r_**2 + K + B*r_**2 + Q*r_**4          # this is F
    H = (P['c0']/r_**2 + P['c1'] + P['c2']*r_**2 + P['c3']*sp.log(r_)
         + J**2/(6*r_**6) + J*K/(3*r_**4) - K**2*sp.log(r_)/(2*r_**2)
         + K*Q*r_**2*sp.log(r_) + B*Q*r_**4 + sp.Rational(19, 48)*Q**2*r_**6)
    return f, H


def build(P):
    f, H = funcs(P)
    X = [sp.Symbol('x%d' % i) for i in range(6)]
    g = sp.zeros(6, 6)
    s2, c2 = sp.sin(th_)**2, sp.cos(th_)**2
    g[0, 1] = g[1, 0] = 1
    g[1, 1] = 1 + 2*H
    g[1, 4] = g[4, 1] = f*s2
    g[1, 5] = g[5, 1] = f*c2
    g[2, 2] = 1
    g[3, 3] = r_**2
    g[4, 4] = r_**2*s2
    g[5, 5] = r_**2*c2
    gX = g.subs({r_: X[2], th_: X[3]})
    gi = gX.inv()
    dg = [gX.diff(X[k]) for k in range(6)]
    Gam = [[[sum(gi[a, d]*(dg[b][d, cc] + dg[cc][d, b] - dg[d][b, cc]) for d in range(6))/2
             for cc in range(6)] for b in range(6)] for a in range(6)]
    fg = sp.lambdify([X], gX, 'numpy')
    fG = sp.lambdify([X], Gam, 'numpy')
    fU = sp.lambdify(r_, -2*H + f**2/r_**2, 'numpy')
    fdU = sp.lambdify(r_, sp.diff(-2*H + f**2/r_**2, r_), 'numpy')
    fH = sp.lambdify(r_, H, 'numpy')
    ff = sp.lambdify(r_, f, 'numpy')
    return (lambda x: np.array(fg(x), dtype=float)), (lambda x: np.array(fG(x), dtype=float)), fU, fH, ff, fdU


def integrate(fG, x0, k0, lmax):
    def rhs(l, y):
        x, k = y[:6], y[6:]
        return np.concatenate([k, -np.einsum('abc,b,c->a', fG(x), k, k)])
    def centre(l, y):
        return y[2] - 1e-2
    centre.terminal = True
    return solve_ivp(rhs, (0, lmax), np.concatenate([x0, k0]), method='DOP853', rtol=1e-11, atol=1e-13,
                     events=centre, max_step=0.2)


def V_of(fH, ff, p, P3, l2):
    return lambda rr: -2*p**2*fH(rr) + 1/rr**2*(l2 - P3**2 + (P3 - p*ff(rr))**2)


def invariants(g, x, k, ff):
    """p = k^u, P = k(chi) = k_phi + k_psi, l^2 from the angular kinetic term"""
    p = k[1]
    kl = g @ k
    P3 = kl[4] + kl[5]
    ang = k[3:] @ g[3:, 3:] @ k[3:]           # h_ab xdot^a xdot^b on the sphere
    rr = x[2]
    l2 = rr**2*ang + P3**2 - (P3 - p*ff(rr))**2
    return p, P3, l2


def main():
    ok = True
    fg, fG, fU, fH, ff, fdU = build(PARS)
    res = minimize_scalar(fU, bounds=(0.5, 50), method='bounded')
    rs = brentq(fdU, 0.9*res.x, 1.1*res.x, xtol=1e-14)
    dU = fdU
    d2U = (fU(rs + 1e-3) - 2*fU(rs) + fU(rs - 1e-3))/1e-6
    print(f'parameters {PARS}')
    print(f'U = -2H + F^2/r^2 has its minimum at r_s = {rs:.6f} (ZAMO lapse N^2 = 1 + U = {1 + fU(rs):.4f}), U\'\'(r_s) = {d2U:.3e};'
          f' U(0.3) = {fU(0.3):.3g}, U(40) = {fU(40.):.3g}')
    ok &= bool(d2U > 0 and abs(dU(rs)) < 1e-6)
    print("  [check]", ok) if not ok else None
    # ZAMO at r_s
    x0 = np.array([0., 0., rs, 0.5, 0., 0.])
    g0 = fg(x0)
    xi = np.array([1., -1, 0, 0, 0, 0])
    e5 = np.zeros(6); e5[4] = e5[5] = 1.      # chi = d_phi + d_psi
    w = -(xi @ g0 @ e5)/(e5 @ g0 @ e5)
    U0 = xi + w*e5
    u0 = U0/np.sqrt(-(U0 @ g0 @ U0))
    sol = integrate(fG, x0, u0, 300.)
    drift = np.max(np.abs(sol.y[2] - rs))
    print(f'ZAMO released at r_s (any point of the sphere): max |r - r_s| over proper time 300 = {drift:.2e}')
    ok &= bool(drift < 1e-6)
    print("  [check]", ok) if not ok else None
    # random geodesics from the orbit, any direction on the sphere and in u, r
    rng = np.random.default_rng(0)
    worst_cons, rmin, rmax, fell = 0., np.inf, 0., 0
    for i in range(40):
        x0 = np.array([0., 0., rs, rng.uniform(0.15, np.pi/2 - 0.15), rng.uniform(0, 2*np.pi), rng.uniform(0, 2*np.pi)])
        g0 = fg(x0)
        w = -(xi @ g0 @ e5)/(e5 @ g0 @ e5)
        U0 = xi + w*e5
        u0 = U0/np.sqrt(-(U0 @ g0 @ U0))
        es = []
        for j in (2, 3, 4, 5, 1):
            e = np.zeros(6); e[j] = 1.
            e = e + (u0 @ g0 @ e)*u0
            for q in es:
                e = e - (q @ g0 @ e)*q
            es.append(e/np.sqrt(e @ g0 @ e))
        nv = rng.normal(size=5); nv /= np.linalg.norm(nv)
        null = i % 2 == 0
        k0 = u0 + (1.0 if null else 0.7)*sum(a*e for a, e in zip(nv, es))
        p, P3, l2 = invariants(g0, x0, k0, ff)
        V = V_of(fH, ff, p, P3, l2)
        E = k0[2]**2 + V(rs)
        sol = integrate(fG, x0, k0, 150.)
        rr = sol.y[2]
        cons = np.max(np.abs(sol.y[8]**2 + V(rr) - E))/max(1., abs(E))
        worst_cons = max(worst_cons, cons)
        rmin, rmax = min(rmin, rr.min()), max(rmax, rr.max())
        fell += sol.status == 1
        # the orbit stays where V(r) <= E, and that region is bounded away from 0 and infinity
        inside = bool(np.max(V(rr) - E) <= 1e-8*max(1., abs(E)) and V(0.05) > E and V(200.) > E)
        if not inside:
            print('  outside the allowed band', np.max(V(rr) - E), V(0.05) - E, V(200.) - E)
        ok &= inside
    print(f'40 random null and timelike geodesics from r_s: r in [{rmin:.3f}, {rmax:.3f}], {fell} reached the centre;'
          f' r\'^2 + V(r) conserved to {worst_cons:.1e}; every orbit stays between the turning points of V')
    ok &= bool(fell == 0 and worst_cons < 1e-6)
    print("  [check]", ok) if not ok else None
    print('PASS no leak into r (1+1+4)' if ok else 'FAIL no leak into r (1+1+4)')
    return ok


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
