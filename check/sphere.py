"""The three-sphere of geodesic monism rotating along phi, 1 + 1 + 4, coordinates (t, u, r, psi, theta, phi).

    ds^2 = 2 dt du + (1 + 2H) du^2 + 2 F(r) mu dphi du + dr^2 + r^2 (dpsi^2 + sin^2 psi (dtheta^2 + sin^2 theta dphi^2)),
    mu = sin^2 psi sin^2 theta = (x1^2 + x2^2)/r^2,   Y = 2 mu - 1  (an l = 2 harmonic of the three-sphere).

Part 1 derives F and H from the reduced equations of the null Kaluza form:
    (ii)  rough Laplacian of the current J = 0,       (iii)  Lap(Lap H - F^2/4) = J.J/2 + F^{ik} nabla_i J_k,
and checks them against the closed forms of the post (and of GeodesicMonismSphere.lean).
Part 2 prints R_uu and R_uphi: the metric is not Ricci flat, and R_uu grows like Q^2 r^4.
Part 3 checks the redshift coefficients of both patches by Taylor expansion in mu.
Part 4 integrates null geodesics in the six-dimensional metric and compares the exact redshift of
Proposition 5 with the three-sphere form of Proposition 6,
    z = C_r dr + C_theta q + cos(gamma) [D_r dr],
    q = vartheta dx sin(gamma) cos(beta)/r - dx^2 s/(2 r^2),
    s = sin^2 gamma (equator: the rotation plane),   s = 1 - sin^2 gamma sin^2 beta (pole: the rotation axis),
with gamma the angle of the step from e_phi and beta its angle, around e_phi, from the observer's offset.
"""
import sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from tensors import Geometry

r = sp.symbols('r', positive=True)
ps, th, ph = sp.symbols('psi theta phi')
J, K, B, Q, c0, c1, c2, c3, m_ = sp.symbols('J K B Q c0 c1 c2 c3 mu', real=True)
L = sp.log(r)
F = J/r**2 + K + B*r**2 + Q*r**4
H0 = (J**2/(12*r**6) + J*K/(6*r**4) - K**2*L/(4*r**2) - J*Q*L + B*K*L + B**2*r**2/4
      + K*Q*r**2*L/2 - K*Q*r**2/8 + B*Q*r**4/2 + sp.Rational(19, 96)*Q**2*r**6)
H2 = (J*K*L/(3*r**4) + K**2*L/(4*r**2) + K**2/(16*r**2) + B*J/(2*r**2) + B*K/2
      + K*Q*r**2*L/6 + B*Q*r**4/8 + sp.Rational(27, 160)*Q**2*r**6)
HR = c0/r**2 + c1 + c2*r**2 + c3*L
MU = sp.sin(ps)**2*sp.sin(th)**2


def Lop(y, l):
    return sp.diff(y, r, 2) + 3*sp.diff(y, r)/r - l*(l + 2)*y/r**2


def part1():
    ok = True
    X = [r, ps, th, ph]
    h = sp.diag(1, r**2, r**2*sp.sin(ps)**2, r**2*MU)
    G = Geometry(h, X, simp=sp.simplify)
    good = G.ricci() == sp.zeros(4, 4)
    print(('PASS' if good else 'FAIL') + ' the transverse space dr^2 + r^2 dOmega_3^2 is flat R^4')
    ok &= good
    A = [0, 0, 0, F*MU]
    hi, Gm = G.ginv, G.Gam
    Fm = sp.Matrix(4, 4, lambda i, j: sp.diff(A[j], X[i]) - sp.diff(A[i], X[j]))
    Fup = (hi*Fm*hi).applyfunc(sp.simplify)
    F4 = sp.simplify(sum(Fm[i, j]*Fup[i, j] for i in range(4) for j in range(4))/4)
    good = sp.simplify(F4 - (r**2*sp.diff(F, r)**2*MU + 4*F**2*(1 - MU))/(2*r**4)) == 0
    print(('PASS' if good else 'FAIL') + ' F^2/4 = (r^2 F\'^2 mu + 4 F^2 (1 - mu))/(2 r^4)')
    ok &= good
    DF = [[[sp.diff(Fm[i, j], X[k]) - sum(Gm[e][k][i]*Fm[e, j] + Gm[e][k][j]*Fm[i, e] for e in range(4))
            for j in range(4)] for i in range(4)] for k in range(4)]
    Jv = [sp.simplify(sum(hi[k, l]*DF[l][k][i] for k in range(4) for l in range(4))) for i in range(4)]
    good = all(sp.simplify(Jv[i] - [0, 0, 0, 4*(3*Q*r**4 - K)*MU/r**2][i]) == 0 for i in range(4))
    print(('PASS' if good else 'FAIL') + ' current J_phi = (F\'\' + F\'/r - 4F/r^2) mu = 4(3Qr^4 - K) mu/r^2')
    ok &= good
    Dw = [[sp.diff(Jv[k], X[i]) - sum(Gm[e][i][k]*Jv[e] for e in range(4)) for k in range(4)] for i in range(4)]
    DDw = [[[sp.diff(Dw[i][k], X[j]) - sum(Gm[e][j][i]*Dw[e][k] + Gm[e][j][k]*Dw[i][e] for e in range(4))
             for k in range(4)] for i in range(4)] for j in range(4)]
    good = all(sp.simplify(sum(hi[i, j]*DDw[j][i][k] for i in range(4) for j in range(4))) == 0 for k in range(4))
    print(('PASS' if good else 'FAIL') + ' (ii) the fourth-order Maxwell equation holds for all four modes')
    ok &= good
    Jup = [sum(hi[i, j]*Jv[j] for j in range(4)) for i in range(4)]
    DJ = [[sp.diff(Jv[k], X[i]) - sum(Gm[e][i][k]*Jv[e] for e in range(4)) for k in range(4)] for i in range(4)]
    src = sp.simplify(sum(Jv[i]*Jup[i] for i in range(4))/2
                      + sum(Fup[i, k]*DJ[i][k] for i in range(4) for k in range(4)))
    Hf = HR + H0 + H2*(2*MU - 1)

    def lap(y):
        return sum(hi[i, j]*(sp.diff(y, X[i], X[j]) - sum(Gm[e][i][j]*sp.diff(y, X[e]) for e in range(4)))
                   for i in range(4) for j in range(4))
    good = sp.simplify(lap(sp.simplify(lap(Hf) - F4)) - src) == 0
    print(('PASS' if good else 'FAIL') + ' (iii) back-reaction with H = H_rad + h0 + h2 Y as in the post')
    ok &= good
    # the l = 2 part: Lap(f(r) Y) = L_2 f Y
    good = sp.simplify(lap(H2*(2*MU - 1)) - Lop(H2, 2)*(2*MU - 1)) == 0
    print(('PASS' if good else 'FAIL') + ' Y = 2 sin^2psi sin^2theta - 1 is an l = 2 harmonic: Lap(f Y) = (f\'\' + 3f\'/r - 8f/r^2) Y')
    ok &= good
    return ok


def part2():
    Y = 2*m_ - 1
    F4 = (r**2*sp.diff(F, r)**2*m_ + 4*F**2*(1 - m_))/(2*r**4)
    Ruu = sp.expand(-(Lop(H0 + HR, 0) + Lop(H2, 2)*Y) + F4)
    R0 = sp.expand(Ruu.subs(m_, sp.Rational(1, 2)))
    R2 = sp.expand((Ruu - R0).subs(m_, 1))
    R0c = (-8*c2 - 2*c3/r**2 + 2*J*K/(3*r**6) + K**2/(2*r**4) - 4*K*Q*L - 6*B*Q*r**2 - sp.Rational(9, 2)*Q**2*r**4)
    R2c = 2*B*K/r**2 - 6*J*Q/r**2 + 2*K**2*L/r**4 - 3*K*Q - sp.Rational(15, 4)*Q**2*r**4
    good = sp.simplify(R0 - R0c) == 0 and sp.simplify(R2 - R2c) == 0
    print(('PASS' if good else 'FAIL') + ' R_uu = R0 + R2 Y,')
    print('       R0 = -8c2 - 2c3/r^2 + 2JK/(3r^6) + K^2/(2r^4) - 4KQ ln r - 6BQ r^2 - 9Q^2 r^4/2,')
    print('       R2 = 2BK/r^2 - 6JQ/r^2 + 2K^2 ln r/r^4 - 3KQ - 15Q^2 r^4/4;  R_uphi = -2(3Qr^4 - K) mu/r^2'
          ' (GeodesicMonismSphere.lean, sphere_ricci)')
    return good


def part3():
    """coefficients in both patches: functions of (r, mu); equator mu = 1 - n^2, pole mu = n^2 (n the angle)."""
    ok = True
    n = sp.Symbol('n')
    Hm = HR + H0 + H2*(2*m_ - 1)
    a2 = 1 - 2*Hm
    N2 = 1 - 2*Hm + F**2*m_/r**2
    He, Hp = HR + H0 + H2, HR + H0 - H2
    Ne2 = 1 - 2*He + F**2/r**2
    claims = {
        ('equator', 'static'): (-sp.diff(He, r)/(1 - 2*He), 4*H2/(1 - 2*He)),
        ('equator', 'ZAMO'): (sp.diff(Ne2, r)/(2*Ne2), (4*H2 - F**2/r**2)/Ne2),
        ('pole', 'static'): (-sp.diff(Hp, r)/(1 - 2*Hp), -4*H2/(1 - 2*Hp)),
        ('pole', 'ZAMO'): (-sp.diff(Hp, r)/(1 - 2*Hp), (F**2/r**2 - 4*H2)/(1 - 2*Hp)),
    }
    for (patch, obs), (Cr, Ct) in claims.items():
        fn = sp.log(a2)/2 if obs == 'static' else sp.log(N2)/2
        sub = 1 - n**2 if patch == 'equator' else n**2
        e = fn.subs(m_, sub)
        Cr_ = sp.diff(e, r).subs(n, 0)
        Ct_ = sp.diff(e, n, 2).subs(n, 0)
        lin = sp.diff(e, n).subs(n, 0)
        vals = {J: 0.7, K: -0.13, B: 0.011, Q: 3e-4, c0: -0.4, c1: 0.02, c2: -0.003, c3: 0.05, r: 2.3}
        good = all(abs(float((x - y).subs(vals))) < 1e-12 for x, y in ((Cr_, Cr), (Ct_, Ct))) and lin == 0
        print(('PASS' if good else 'FAIL') + f' {patch:7s} {obs:6s}: C_r and C_theta as in the post')
        ok &= good
    # ZAMO rotation term: omega = F/r^2 does not depend on the angles, D_theta = 0, D_r = (r/N) omega' on the equator
    return ok


# ------------------------------------------------------------- part 4: integrated rays
PARS = dict(J=0.5, K=-0.2, B=0.005, Q=-1e-4, c0=-1.0, c1=0.0, c2=-0.005, c3=0.02)


def build(P):
    Jn, Kn, Bn, Qn = P['J'], P['K'], P['B'], P['Q']
    Fn = Jn/r**2 + Kn + Bn*r**2 + Qn*r**4
    sub = {J: Jn, K: Kn, B: Bn, Q: Qn, c0: P['c0'], c1: P['c1'], c2: P['c2'], c3: P['c3']}
    Hn = (HR + H0 + H2*(2*MU - 1)).subs(sub)
    g = sp.zeros(6, 6)
    g[0, 1] = g[1, 0] = 1
    g[1, 1] = 1 + 2*Hn
    g[1, 5] = g[5, 1] = Fn*MU
    g[2, 2] = 1
    g[3, 3] = r**2
    g[4, 4] = r**2*sp.sin(ps)**2
    g[5, 5] = r**2*MU
    X = [sp.Symbol('x%d' % i) for i in range(6)]
    gX = g.subs({r: X[2], ps: X[3], th: X[4]})
    gi = gX.inv()
    dg = [gX.diff(X[c]) for c in range(6)]
    Gam = [[[sum(gi[a, d]*(dg[b][d, c] + dg[c][d, b] - dg[d][b, c]) for d in range(6))/2
             for c in range(6)] for b in range(6)] for a in range(6)]
    fg = sp.lambdify([X], gX, 'numpy')
    fG = sp.lambdify([X], Gam, 'numpy')
    co = {}
    for k, v in dict(F=Fn, He=(HR + H0 + H2).subs(sub), Hp=(HR + H0 - H2).subs(sub), h2=H2.subs(sub)).items():
        co[k] = sp.lambdify(r, v, 'numpy')
        co[k + "'"] = sp.lambdify(r, sp.diff(v, r), 'numpy')
    return (lambda x: np.array(fg(x), dtype=float)), (lambda x: np.array(fG(x), dtype=float)), co


def coeffs(co, patch, obs, rr):
    F_, Fp_, h2 = co['F'](rr), co["F'"](rr), co['h2'](rr)
    if patch == 'equator':
        He, Hep = co['He'](rr), co["He'"](rr)
        if obs == 'static':
            return -Hep/(1 - 2*He), 4*h2/(1 - 2*He), 0.
        Ne2 = 1 - 2*He + F_**2/rr**2
        dNe2 = -2*Hep + (2*F_*Fp_*rr**2 - 2*rr*F_**2)/rr**4
        Dr = (rr/np.sqrt(Ne2))*(Fp_/rr**2 - 2*F_/rr**3)
        return dNe2/(2*Ne2), (4*h2 - F_**2/rr**2)/Ne2, Dr
    Hp, Hpp = co['Hp'](rr), co["Hp'"](rr)
    Ct = -4*h2 if obs == 'static' else F_**2/rr**2 - 4*h2
    return -Hpp/(1 - 2*Hp), Ct/(1 - 2*Hp), 0.


def cart(x):
    """unit vector in R^4: x1 + i x2 = sin psi sin theta e^{i phi}, x3 = sin psi cos theta, x4 = cos psi"""
    R_, p_, t_, f_ = x[2], x[3], x[4], x[5]
    return np.array([np.sin(p_)*np.sin(t_)*np.cos(f_), np.sin(p_)*np.sin(t_)*np.sin(f_),
                     np.sin(p_)*np.cos(t_), np.cos(p_)])


def run(fg, fG, co, patch, obs, scale, nhat):
    xi = np.array([1., -1, 0, 0, 0, 0])
    eph = np.zeros(6); eph[5] = 1.
    rO = 4.0
    # emitter E: equator = near the rotation plane (psi = theta = pi/2), offset in psi and theta;
    # pole = near the rotation axis (theta = 0), at psi = 1.1, offset in theta
    if patch == 'equator':
        x0 = np.array([0, 0, rO, np.pi/2 + 0.08*scale, np.pi/2, 0.3])
    else:
        x0 = np.array([0, 0, rO, 1.1, 0.03*scale, 0.3])
    g0 = fg(x0)

    def obsU(g):
        if obs == 'static':
            U = xi
        else:
            w = -(xi @ g @ eph)/(eph @ g @ eph)
            U = xi + w*eph
        return U/np.sqrt(-(U @ g @ U))
    U0 = obsU(g0)
    es = []
    for j in (2, 3, 4, 5):
        e = np.zeros(6); e[j] = 1.
        e = e + (U0 @ g0 @ e)*U0
        for f in es:
            e = e - (f @ g0 @ e)*f
        es.append(e/np.sqrt(e @ g0 @ e))
    nv = np.array(nhat, dtype=float); nv[0] *= scale; nv /= np.linalg.norm(nv)
    k0 = U0 + sum(a*e for a, e in zip(nv, es))
    nphi = nv[3]                             # e_phi is the last Gram-Schmidt vector, along d_phi

    def rhs(l, y):
        x, k = y[:6], y[6:]
        return np.concatenate([k, -np.einsum('abc,b,c->a', fG(x), k, k)])
    sol = solve_ivp(rhs, (0, 0.6*scale), np.concatenate([x0, k0]), method='DOP853', rtol=1e-12, atol=1e-14)
    x1 = sol.y[:6, -1]
    g1 = fg(x1)
    E = -(k0 @ g0 @ xi); Lph = k0 @ g0 @ eph
    b = Lph/E

    def NW(g):
        if obs == 'static':
            return np.sqrt(-(xi @ g @ xi)), 0.
        w = -(xi @ g @ eph)/(eph @ g @ eph)
        U = xi + w*eph
        return np.sqrt(-(U @ g @ U)), w
    N0, w0 = NW(g0); N1, w1 = NW(g1)
    z_exact = (N1/N0)*(1 - w0*b)/(1 - w1*b) - 1
    # local coordinates at the observer O = x1 (all lengths at r = r_O, ell = r_O)
    rr = x1[2]
    dr = x1[2] - x0[2]
    vO, vE = cart(x1), cart(x0)
    if patch == 'equator':
        # rotation plane x1 x2; offsets P = r (x3, x4)
        PO, PE = rr*vO[2:], rr*vE[2:]
        ephi = np.array([-np.sin(x1[5]), np.cos(x1[5])])
        dpar = rr*(vO[:2] - vE[:2]) @ ephi
        dP = PO - PE
        dx = np.sqrt(dpar**2 + dP @ dP)
        cg = dpar/dx
        sg = np.sqrt(max(0., 1 - cg**2))
        cb = (dP @ PO)/(np.linalg.norm(dP)*np.linalg.norm(PO))
        vt = np.linalg.norm(PO)/rr
        s_ = sg**2
    else:
        # axis plane x3 x4 (a circle on the sphere); offsets P = r (x1, x2)
        PO, PE = rr*vO[:2], rr*vE[:2]
        dP = PO - PE
        ephi = np.array([-PO[1], PO[0]])/np.linalg.norm(PO)
        ehat_rho = PO/np.linalg.norm(PO)
        dcirc = rr*np.linalg.norm(vO[2:] - vE[2:])
        dx = np.sqrt(dP @ dP + dcirc**2)
        cg = (dP @ ephi)/dx
        sg = np.sqrt(max(0., 1 - cg**2))
        cb = (dP @ ehat_rho)/(dx*sg)
        vt = np.linalg.norm(PO)/rr
        s_ = 1 - sg**2*(1 - cb**2)
    Cr, Ct, Dr = coeffs(co, patch, obs, rr)
    q = vt*dx*sg*cb/rr - dx**2*s_/(2*rr**2)
    z_local = Cr*dr + Ct*q + nphi*Dr*dr
    return z_exact, z_local, dx, dr, np.degrees(np.arccos(np.clip(cg, -1, 1))), np.degrees(np.arccos(np.clip(cb, -1, 1)))


def part4():
    ok = True
    fg, fG, co = build(PARS)
    for patch in ('equator', 'pole'):
        for obs in ('static', 'ZAMO'):
            for nhat in ((0.05, 0.1, 0.6, 0.75), (0.05, -0.2, 0.8, 0.5)):
                errs, d0 = [], None
                for s_ in (1.0, 0.5, 0.25):
                    ze, zl, dx, dr, gam, bet = run(fg, fG, co, patch, obs, s_, nhat)
                    errs.append(abs(ze - zl))
                    if d0 is None:
                        d0 = (ze, zl, dx, dr, gam, bet)
                ratio = errs[1]/errs[2]
                print(f'  {patch:7s} {obs:6s} dx={d0[2]:.4f} gamma={d0[4]:5.1f} beta={d0[5]:5.1f} dr={d0[3]:+.1e}:'
                      f' z exact {d0[0]:+.6e}, Prop. 6 {d0[1]:+.6e}; error ratio on halving {errs[0]/errs[1]:.1f}, {ratio:.1f}')
                ok &= bool(ratio > 12)   # the error falls like eps^4 (ratio 16); a wrong q gives eps^3 or worse
    print(('PASS' if ok else 'FAIL') + ' z(dx, gamma, beta) on the three-sphere: the error falls by 16 per halving (eps^4) in both patches')
    return ok


def main():
    ok = part1()
    ok &= part2()
    ok &= part3()
    ok &= part4()
    print('PASS three-sphere rotating along phi' if ok else 'FAIL three-sphere rotating along phi')
    return ok


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
