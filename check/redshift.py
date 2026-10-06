"""Redshift of a short climbing light ray, both patches, three metrics.

Setting.  A stationary, axisymmetric metric with a timelike Killing vector xi
(the clock) and an axial Killing vector d_phi.  Emitter E, observer O, joined
by a null geodesic with wave vector k.  tau_i is the emitted period, tau_f the
received period, both proper times of the chosen observers.

  tau_f / tau_i = (k.U)_E / (k.U)_O = 1 + z                      (Theorem 3)

Observer families
  static  U = xi / alpha,            alpha^2 = -g(xi,xi)
  ZAMO    U = (xi + w d_phi) / N,    w = -g(xi,d_phi)/g_phiphi,
                                     N^2 = -g(xi,xi) + g(xi,d_phi)^2/g_phiphi
Proposition 5 (Killing charges E = -k.xi, L = k.d_phi):
  static : ln(1+z) = (F(O) - F(E))/2,        F = ln alpha^2
  ZAMO   : ln(1+z) = (G(O) - G(E))/2 + ln(1 - w_E b) - ln(1 - w_O b),
           G = ln N^2,  b = L/E,  b/(1 - w_E b) = sqrt(g_phiphi) n_phi / N at E,
           n_phi = cosine between the photon direction and e_phi in the ZAMO frame.

Regime  r >> Delta x >> delta r, counted with eps = Delta x / r:
  Delta theta, Delta phi, observer offset vartheta ~ eps,   delta r / r ~ eps^2.
Coordinates: O = (r, theta0 + vartheta), E = (r - delta r, theta0 + vartheta - Delta theta).
Equator patch theta0 = pi/2, pole patch theta0 = 0.

Taylor expansion to O(eps^2), with X a function of (r, theta):
  X(O) - X(E) = Dth X_th + delta_r X_r + (vt Dth - Dth^2/2) X_thth + O(eps^3)
and the ZAMO cross term is (sqrt(g_phiphi) n_phi / N)(w(O) - w(E)).  At the
pole sqrt(g_phiphi) = O(eps), so the cross term is O(eps^3) there.
All metric functions are written in c = cos(theta); d/dtheta = -sin(theta) d/dc.
"""
import sympy as sp

r = sp.symbols('r', positive=True)
cth = sp.symbols('c', real=True)                     # c = cos(theta)
dr, dth, vt, nphi = sp.symbols('delta_r Delta_theta vartheta n_phi', real=True)
M, a, mu, J = sp.symbols('M a mu J', real=True)
c0, c1, c2, c3 = sp.symbols('c0 c1 c2 c3', real=True)
Kk, Bb, Qq = sp.symbols('K B Q', real=True)
S2 = 1 - cth**2                                      # sin^2(theta)


def kerr():
    rho2 = r**2 + a**2*cth**2
    return (-(1 - 2*M*r/rho2), -2*M*a*r*S2/rho2,
            (r**2 + a**2 + 2*M*a**2*r*S2/rho2)*S2)


def mp5():
    rho2 = r**2 + a**2*cth**2
    return (-1 + mu/rho2, -mu*a*S2/rho2, (r**2 + a**2)*S2 + mu*a**2*S2**2/rho2)


def gm_H():
    P2 = (3*cth**2 - 1)/2
    return c0/r + c1 + c2*r + c3*r**2 + J**2/(12*r**4) + J**2*P2/(6*r**4)


def gm(lmb=-1):
    """Null Kaluza hydrogen, clock xi = d_t + lambda d_u."""
    H = gm_H()
    A = J*S2/r
    return (2*lmb + lmb**2*(1 + 2*H), lmb*A, r**2*S2)


def rot_parts():
    """The rotating 1+4 solution: f, h0, h2 (H = Hrad + h0 + h2 P2, A_phi = f sin^2)."""
    f = J/r + Kk*r + Bb*r**2 + Qq*r**4
    h0 = (J**2/(12*r**4) + J*Kk/(6*r**2) - 2*J*Qq*r/3 + Kk**2*sp.log(r) + 4*Bb*Kk*r/3 + Bb**2*r**2/3
          + Kk*Qq*r**3 + 2*Bb*Qq*r**4/3 + sp.Rational(11, 42)*Qq**2*r**6)
    h2 = (J**2/(6*r**4) - 2*J*Kk/(3*r**2) - 2*J*Bb/(3*r) + J*Qq*r/3 - Kk**2/2 - 2*Bb*Kk*r/3
          - 2*Kk*Qq*r**3/3 - sp.Rational(4, 21)*Bb*Qq*r**4 - sp.Rational(29, 126)*Qq**2*r**6)
    Hrad = c0/r + c1 + c2*r + c3*r**2
    return f, h0, h2, Hrad


def gmrot(lmb=-1):
    """Rotating 1+4 geodesic-monism solution, clock xi = d_t - d_u."""
    f, h0, h2, Hrad = rot_parts()
    H = Hrad + h0 + h2*(3*cth**2 - 1)/2
    A = f*S2
    return (2*lmb + lmb**2*(1 + 2*H), lmb*A, r**2*S2)


def derivs(X, theta0):
    """X_r, X_theta, X_thetatheta at (r, theta0) for X(r, c)."""
    c_val = 0 if theta0 == 'equator' else 1
    s_val = 1 if theta0 == 'equator' else 0
    Xc, Xcc = sp.diff(X, cth), sp.diff(X, cth, 2)
    X_r = sp.diff(X, r).subs(cth, c_val)
    X_th = (-s_val*Xc).subs(cth, c_val)
    X_thth = (s_val**2*Xcc - c_val*Xc).subs(cth, c_val)   # d2/dth2 = s^2 d2/dc2 - c d/dc
    return [sp.simplify(e) for e in (X_r, X_th, X_thth)]


def difference(X, patch):
    Xr, Xt, Xtt = derivs(X, patch)
    return dth*Xt + dr*Xr + (vt*dth - dth**2/2)*Xtt


def coefficients(metric, patch, **kw):
    """z = C_r dr + C_t (vt*Dth - Dth^2/2) [+ n_phi (D_r dr + D_t (vt*Dth - Dth^2/2))]."""
    gxx, gxp, gpp = [sp.cancel(e) for e in metric(**kw)]
    alpha2 = sp.cancel(-gxx)
    w = sp.cancel(-gxp/gpp)
    N2 = sp.cancel(-gxx + gxp**2/gpp)
    out = {}
    for key, F in (('static', sp.log(alpha2)), ('ZAMO', sp.log(N2))):
        Fr, Ft, Ftt = derivs(F, patch)
        assert Ft == 0
        out[key] = dict(C_r=sp.factor(Fr/2), C_t=sp.factor(Ftt/2))
    if patch == 'equator':
        pref = sp.sqrt(sp.factor(sp.cancel((gpp/N2).subs(cth, 0))))
        wr, wt, wtt = derivs(w, patch)
        assert wt == 0
        out['ZAMO'].update(D_r=sp.factor(sp.simplify(pref*wr)), D_t=sp.factor(sp.simplify(pref*wtt)))
    else:
        out['ZAMO'].update(D_r=0, D_t=0)
    return out


def ell(name, patch):
    """sqrt(g_thetatheta) at the patch centre: the proper length of one radian of theta."""
    if name.startswith('Geodesic'):
        return r      # flat transverse space: g_thetatheta = r^2
    return r if patch == 'equator' else sp.sqrt(r**2 + a**2)


def z_local(co, patch, l, dr_, dx, gam, vt_, nphi=None):
    """Proposition 6 in the observer's local coordinates.
    dx: proper horizontal step from emitter to observer, gam: its angle from the
    rotation direction e_phi (towards e_theta), vt_: observer offset from the patch
    centre in radians of theta, dr_: the climb r_O - r_E.
        equator: Delta theta = dx sin(gam)/l
        pole:    (theta_O^2 - theta_E^2)/2 = vt dx sin(gam)/l - dx^2/(2 l^2)"""
    s2 = sp.sin(gam)**2 if patch == 'equator' else 1
    q = vt_*dx*sp.sin(gam)/l - dx**2*s2/(2*l**2)
    n = sp.cos(gam) if nphi is None else nphi     # photon direction cosine along e_phi (ZAMO frame)
    return (co['C_r']*dr_ + co['C_t']*q
            + n*(co.get('D_r', 0)*dr_ + co.get('D_t', 0)*q))


def leading(e, n=2):
    """first terms of the large-r expansion"""
    x = sp.symbols('x', positive=True)
    return sp.simplify(sp.series(e.subs(r, 1/x), x, 0, n).removeO().subs(x, 1/r))


# Closed forms quoted in the post and certified in Kerr.lean, MyersPerry.lean,
# GeodesicMonism.lean.  Keys: (metric, patch, observer, coefficient).
def closed_forms():
    Dl = r**2 - 2*M*r + a**2
    S = r**3 + a**2*r + 2*M*a**2
    T = r**4 + a**2*r**2 + a**2*mu
    Dm = r**2 + a**2 - mu
    H0 = c0/r + c1 + c2*r + c3*r**2
    H0p = sp.diff(H0, r)
    a0 = 1 - 2*H0
    Ne2 = a0 + J**2/r**4
    hp = H0 + J**2/(4*r**4)
    ap = 1 - 2*hp
    f, h0, h2, Hrad = rot_parts()
    He = Hrad + h0 - h2/2
    Hp = Hrad + h0 + h2
    ae, ap_ = 1 - 2*He, 1 - 2*Hp
    Nr2 = ae + f**2/r**2
    RT = 'Geodesic monism, rotating 1+4'
    rot = {
        (RT, 'equator', 'static'): dict(C_r=-sp.diff(He, r)/ae, C_t=-3*h2/ae),
        (RT, 'equator', 'ZAMO'): dict(C_r=sp.diff(Nr2, r)/(2*Nr2), C_t=(-3*h2 - f**2/r**2)/Nr2,
                                      D_r=r*sp.diff(f/r**2, r)/sp.sqrt(Nr2), D_t=0),
        (RT, 'pole', 'static'): dict(C_r=-sp.diff(Hp, r)/ap_, C_t=3*h2/ap_),
        (RT, 'pole', 'ZAMO'): dict(C_r=-sp.diff(Hp, r)/ap_, C_t=(3*h2 + f**2/r**2)/ap_, D_r=0, D_t=0),
    }
    out = {
        ('Kerr', 'equator', 'static'): dict(C_r=M/(r*(r - 2*M)), C_t=2*M*a**2/(r**2*(r - 2*M))),
        ('Kerr', 'equator', 'ZAMO'): dict(C_r=M*(r**4 + 2*a**2*r**2 + a**4 - 4*M*a**2*r)/(r*S*Dl),
                                         C_t=2*M*a**2*(r**2 + a**2)/(r**2*S),
                                         D_r=-2*M*a*(3*r**2 + a**2)/(r*sp.sqrt(Dl)*S),
                                         D_t=-4*M*a**3*sp.sqrt(Dl)/(r**2*S)),
        ('Kerr', 'pole', 'static'): dict(C_r=M*(r**2 - a**2)/((r**2 + a**2)*Dl), C_t=-2*M*a**2*r/((r**2 + a**2)*Dl)),
        ('Kerr', 'pole', 'ZAMO'): dict(C_r=M*(r**2 - a**2)/((r**2 + a**2)*Dl), C_t=-2*M*a**2*r/(r**2 + a**2)**2,
                                      D_r=0, D_t=0),
        ('Myers-Perry 5D', 'equator', 'static'): dict(C_r=mu/(r*(r**2 - mu)), C_t=a**2*mu/(r**2*(r**2 - mu))),
        ('Myers-Perry 5D', 'equator', 'ZAMO'): dict(C_r=mu*(r**4 + 2*a**2*r**2 + a**4 - a**2*mu)/(r*Dm*T),
                                                   C_t=a**2*mu*(r**2 + a**2)/(r**2*T),
                                                   D_r=-2*a*mu*(2*r**2 + a**2)/(r*sp.sqrt(Dm)*T),
                                                   D_t=-2*a**3*mu*sp.sqrt(Dm)/(r**2*T)),
        ('Myers-Perry 5D', 'pole', 'static'): dict(C_r=mu*r/((r**2 + a**2)*Dm), C_t=-a**2*mu/((r**2 + a**2)*Dm)),
        ('Myers-Perry 5D', 'pole', 'ZAMO'): dict(C_r=mu*r/((r**2 + a**2)*Dm), C_t=-a**2*mu/(r**2 + a**2)**2,
                                                D_r=0, D_t=0),
        ('Geodesic monism (xi = d_t - d_u)', 'equator', 'static'): dict(C_r=-H0p/a0, C_t=-J**2/(2*r**4*a0)),
        ('Geodesic monism (xi = d_t - d_u)', 'equator', 'ZAMO'): dict(C_r=(-H0p - 2*J**2/r**5)/Ne2,
                                                                     C_t=-3*J**2/(2*r**4*Ne2),
                                                                     D_r=-3*J/(r**3*sp.sqrt(Ne2)), D_t=0),
        ('Geodesic monism (xi = d_t - d_u)', 'pole', 'static'): dict(C_r=-sp.diff(hp, r)/ap, C_t=J**2/(2*r**4*ap)),
        ('Geodesic monism (xi = d_t - d_u)', 'pole', 'ZAMO'): dict(C_r=-sp.diff(hp, r)/ap, C_t=3*J**2/(2*r**4*ap),
                                                                  D_r=0, D_t=0),
    }
    out.update(rot)
    return out


def numerically_equal(e1, e2, syms, trials=4):
    """compare at random points in the physical domain (r large, positive params)"""
    import random
    random.seed(7)
    for _ in range(trials):
        vals = {s: sp.Rational(random.randint(1, 9), random.randint(10, 40)) for s in syms}
        vals[r] = sp.Integer(random.randint(20, 60))
        vals[c2] = vals[c2]/1000
        vals[c3] = vals[c3]/10**6
        vals[Kk] = vals[Kk]/100
        vals[Bb] = vals[Bb]/10**4
        vals[Qq] = vals[Qq]/10**8
        d = sp.N((sp.sympify(e1) - sp.sympify(e2)).subs(vals), 30)
        if abs(d) > 1e-20:
            return False
    return True


METRICS = [('Kerr', kerr, 4), ('Myers-Perry 5D', mp5, 5), ('Geodesic monism (xi = d_t - d_u)', gm, 6),
           ('Geodesic monism, rotating 1+4', gmrot, 6)]

if __name__ == '__main__':
    import sys
    ok = True
    CF = closed_forms()
    syms = [M, a, mu, J, c0, c1, c2, c3, Kk, Bb, Qq]
    for name, met, order in METRICS:
        for patch in ('equator', 'pole'):
            co = coefficients(met, patch)
            print(f'== {name}, {patch}')
            for obs in ('static', 'ZAMO'):
                for k, v in co[obs].items():
                    good = numerically_equal(v, CF[(name, patch, obs)][k], syms)
                    ok &= good
                    print(f'  {"PASS" if good else "FAIL"} closed form {name} / {patch} / {obs} / {k}')
                    print(f'  {obs:6s} {k}: {v}')
                    if v != 0 and not name.startswith('Geodesic monism'):
                        print(f'  {"":6s} {k} large r: {leading(v, order)}')
    print('PASS all closed forms' if ok else 'FAIL closed forms')
    sys.exit(0 if ok else 1)
