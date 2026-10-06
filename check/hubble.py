"""The Hubble law of the non-Ricci-flat modes: every source at the centre of its own hydrogen solution.

Each emitter (each hydrogen atom) is the centre of its own solution, and the light climbs out of it to an
observer at distance d.  Along the ray (J = 0):
    ds^2 = 2 dt du + (1 + 2H) du^2 + dr^2,   H = c0/r + c1 + c2 r + c3 r^2,   c2 < 0.
The emitter's own c0/r well gives a fixed, distance-independent factor (like the redshift from a star's
surface); the part that grows with distance comes from c2 d + c3 d^2, so below H_E = c1, H_O = c1 + c2 d + c3 d^2.
A photon has conserved k_t = q, k_u = E_u, and 1 + z = (k.U)_E/(k.U)_O (Theorem 3).
  Killing (static) observers, xi = d_t - d_u:   1 + z = alpha_O/alpha_E,  alpha^2 = 1 - 2H.
  Matter at rest in space drifting in u at u' = sigma (u' is conserved along geodesics), photon with E_u = 0:
      1 + z = (1 + (1 + 2H_E) sigma^2)/(1 + (1 + 2H_O) sigma^2).
Results (D = 1 + (1 + 2c1) sigma^2, a = 1 - 2c1):
  slope   H0/c = -c2/a (static),  -2 sigma^2 c2/D (drifting): positive for c2 < 0
  q0      = -2 - 2 c3 a/c2^2 (static),   1 - c3 D/(sigma^2 c2^2) (drifting)
  drifting, c3 = 0:  1 + z = 1/(1 - (H0/c) d) exactly, a horizon at d = c/H0
  dark sky (shells n L dd/(1+z)^2): drifting, c3 = 0: finite, n L c/(3 H0);
          static, c3 = 0: (1+z)^2 = 1 + 2(H0/c) d, the sum grows like ln d (no horizon)
  sharp images: on the screen of a radial ray the c2, c3 modes are (c2/r + 2 c3) times the identity.
"""
import sys
import sympy as sp

d, c1, c2, c3, sg, q, Eu = sp.symbols('d c1 c2 c3 sigma q E_u', real=True)


def q0_of(z):
    """z(d) = h d + k d^2 + ...;  d_L = d(1+z) = (1/h)[z + (1 - q0) z^2/2] + ..."""
    s = sp.series(z, d, 0, 3).removeO()
    h, k = s.coeff(d, 1), s.coeff(d, 2)
    return sp.simplify(h), sp.simplify(1 - 2*(1 - k/h**2))


def main():
    ok = True
    HE = c1
    HO = c1 + c2*d + c3*d**2
    a = 1 - 2*c1
    D = 1 + (1 + 2*c1)*sg**2
    # static (Killing) observers
    zK = sp.sqrt((1 - 2*HO)/(1 - 2*HE)) - 1
    hK, qK = q0_of(zK)
    # drifting matter: k.U from the normalisation 2 t' sigma + (1 + 2H) sigma^2 = -1, then E_u = 0
    def kU(Hx):
        tdot = -(1 + (1 + 2*Hx)*sg**2)/(2*sg)
        return q*tdot + Eu*sg
    zD = sp.simplify((kU(HE)/kU(HO)).subs(Eu, 0)) - 1
    good = sp.simplify(zD + 1 - (1 + (1 + 2*HE)*sg**2)/(1 + (1 + 2*HO)*sg**2)) == 0
    print(('PASS' if good else 'FAIL') + ' drifting matter with E_u = 0 gives the hydrogen note\'s 1 + z')
    ok &= good
    hD, qD = q0_of(zD)
    good = (sp.simplify(hK + c2/a) == 0 and sp.simplify(hD + 2*sg**2*c2/D) == 0
            and sp.simplify(qK - (-2 - 2*c3*a/c2**2)) == 0 and sp.simplify(qD - (1 - c3*D/(sg**2*c2**2))) == 0)
    print(('PASS' if good else 'FAIL') + f' slopes H0/c: static {sp.factor(hK)}, drifting {sp.factor(hD)} (positive for c2 < 0)')
    print(f'       q0: static {sp.simplify(qK)},  drifting {sp.simplify(qD)}')
    ok &= good
    # drifting, c3 = 0: 1 + z = 1/(1 - h d) exactly, horizon at d = 1/h = c/H0
    h = sp.Symbol('h', positive=True)
    zD0 = sp.simplify(zD.subs(c3, 0))
    good = sp.simplify(zD0 + 1 - 1/(1 - hD.subs(c3, 0)*d)) == 0
    print(('PASS' if good else 'FAIL') + ' drifting, c3 = 0: 1 + z = 1/(1 - (H0/c) d), a horizon at d = c/H0')
    ok &= good
    # static, c3 = 0: (1 + z)^2 = 1 + 2 (H0/c) d
    good = sp.simplify((zK.subs(c3, 0) + 1)**2 - (1 + 2*hK*d)) == 0
    print(('PASS' if good else 'FAIL') + ' static, c3 = 0: (1 + z)^2 = 1 + 2 (H0/c) d, no horizon')
    ok &= good
    # the dark sky (Olbers): shells send n L dd/(1+z)^2 (d_L = d(1+z) in the Euclidean sections)
    Dd = sp.Symbol('Dd', positive=True)
    drift = sp.integrate((1 - h*d)**2, (d, 0, 1/h))
    static = sp.integrate(1/(1 + 2*h*d), (d, 0, Dd))
    good = sp.simplify(drift - 1/(3*h)) == 0 and sp.limit(static, Dd, sp.oo) == sp.oo
    print(('PASS' if good else 'FAIL') + f' dark sky: drifting matter sums to {drift} (finite); static observers to {sp.simplify(static)},'
          ' which grows like ln of the depth')
    ok &= good
    # sharp images: on the screen of a radial ray the c2 and c3 modes are (c2/r + 2 c3) times the identity
    x, y, zz = sp.symbols('x y z', real=True)
    R = sp.sqrt(x**2 + y**2 + zz**2)
    Hc = c2*R + c3*R**2
    Hess = sp.Matrix(3, 3, lambda i, j: sp.diff(Hc, [x, y, zz][i], [x, y, zz][j]))
    pt = {x: sp.Rational(3, 7), y: sp.Rational(-2, 7), zz: sp.Rational(6, 7)}      # |n| = 1, r = 1
    nvec = sp.Matrix([pt[x], pt[y], pt[zz]])
    e1 = sp.Matrix([2, 3, 0])/sp.sqrt(13)
    e2 = nvec.cross(e1)
    S = sp.Matrix(2, 2, lambda a_, b_: ([e1, e2][a_].T*Hess.subs(pt)*[e1, e2][b_])[0])
    good = sp.simplify(S - (c2 + 2*c3)*sp.eye(2)) == sp.zeros(2, 2)
    print(('PASS' if good else 'FAIL') + ' on the screen of a radial ray the c2, c3 modes are (c2/r + 2c3) times the identity: focusing, no shear')
    ok &= good
    print('PASS Hubble law of the non-Ricci-flat modes' if ok else 'FAIL Hubble law')
    return ok


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
