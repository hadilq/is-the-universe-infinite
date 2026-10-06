"""Check that the three metrics of the post solve their field equations.

Kerr and Myers-Perry: Ricci flat (vacuum Einstein).
Geodesic monism hydrogen: E_ab = 0 of S = int R_ab R^ab sqrt(-g).
"""
import sys, time
import sympy as sp
from tensors import Geometry

t, u, r, ph, psi = sp.symbols('t u r phi psi', real=True)
x = sp.symbols('x', real=True)  # x = cos(theta); all components become rational
th = x
M, a, mu, J = sp.symbols('M a mu J', real=True)
c0, c1, c2, c3 = sp.symbols('c0 c1 c2 c3', real=True)
ok = True


def kerr_metric():
    rho2 = r**2 + a**2 * x**2
    Delta = r**2 - 2*M*r + a**2
    g = sp.zeros(4, 4)
    g[0, 0] = -(1 - 2*M*r/rho2)
    g[0, 3] = g[3, 0] = -2*M*a*r*(1-x**2)/rho2
    g[3, 3] = (r**2 + a**2 + 2*M*a**2*r*(1-x**2)/rho2)*(1-x**2)
    g[1, 1] = rho2/Delta
    g[2, 2] = rho2/(1-x**2)
    return g, [t, r, th, ph]


def mp5_metric(post_grr=False):
    rho2 = r**2 + a**2 * x**2
    g = sp.zeros(5, 5)
    g[0, 0] = -1 + mu/rho2
    g[0, 3] = g[3, 0] = -mu*a*(1-x**2)/rho2
    g[3, 3] = (r**2 + a**2)*(1-x**2) + mu*a**2*(1-x**2)**2/rho2
    if post_grr:   # negative control: rho^2/Delta, Delta = r^2 (r^2+a^2-mu)
        g[1, 1] = rho2/(r**2*(r**2 + a**2 - mu))
    else:          # correct: r^2 rho^2/Delta  =  rho^2/(r^2+a^2-mu)
        g[1, 1] = rho2/(r**2 + a**2 - mu)
    g[2, 2] = rho2/(1-x**2)
    g[4, 4] = r**2*x**2
    return g, [t, r, th, ph, psi]


def gm_hydrogen(p2_power=4, coef12=sp.Rational(1, 12), coef6=sp.Rational(1, 6)):
    P2 = (3*x**2 - 1)/2
    H = c0/r + c1 + c2*r + c3*r**2 + coef12*J**2/r**4 + coef6*J**2*P2/r**p2_power
    A = J*(1-x**2)/r
    g = sp.zeros(5, 5)          # (t, u, r, x=cos theta, phi)
    g[0, 1] = g[1, 0] = 1
    g[1, 1] = 1 + 2*H
    g[1, 4] = g[4, 1] = A
    g[2, 2] = 1
    g[3, 3] = r**2/(1-x**2)
    g[4, 4] = r**2*(1-x**2)
    # inverse, written out (checked below): g^{tu}=1, g^{tt}=-(1+2H)+A_i A^i, g^{t phi}=-A^phi
    Aup = A/(r**2*(1-x**2))
    gi = sp.zeros(5, 5)
    gi[0, 1] = gi[1, 0] = 1
    gi[0, 0] = -(1 + 2*H) + A*Aup
    gi[0, 4] = gi[4, 0] = -Aup
    gi[2, 2] = 1
    gi[3, 3] = (1-x**2)/r**2
    gi[4, 4] = 1/(r**2*(1-x**2))
    assert sp.simplify(g*gi - sp.eye(5)) == sp.zeros(5, 5)
    return g, [t, u, r, th, ph], gi


def ricci_at_random_points(g, X, params, npts=3, seed=3):
    """Exact Ricci components at random rational points (Schwartz-Zippel):
    Christoffels and their derivatives are formed symbolically, then every
    symbol is replaced by a random rational and the result evaluated exactly."""
    import random
    random.seed(seed)
    G = Geometry(g, X, simp=lambda e: e)
    n = len(X)
    Gm, x = G.Gam, G.x
    dG = [[[[sp.diff(Gm[a][b][c], x[d]) for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
    worst = 0
    for _ in range(npts):
        v = {X[1]: sp.Rational(random.randint(30, 90), random.randint(1, 3)),     # r
             X[2]: sp.Rational(random.randint(1, 9), random.randint(10, 20))}     # x = cos theta
        for p in params:
            v[p] = sp.Rational(random.randint(1, 9), random.randint(2, 7))
        Gv = [[[Gm[a][b][c].subs(v) for c in range(n)] for b in range(n)] for a in range(n)]
        for b in range(n):
            for d in range(b, n):
                R = sum(dG[a][d][b][a].subs(v) - dG[a][a][b][d].subs(v)
                        + sum(Gv[a][a][e]*Gv[e][d][b] - Gv[a][d][e]*Gv[e][a][b] for e in range(n))
                        for a in range(n))
                worst = max(worst, abs(sp.nsimplify(R)))
    return worst


def report(name, cond):
    global ok
    ok &= bool(cond)
    print(('PASS ' if cond else 'FAIL ') + name)


if __name__ == '__main__':
    t0 = time.time()
    g, X = kerr_metric()
    report('Kerr is Ricci flat (exact zero at random rational points)', ricci_at_random_points(g, X, [M, a]) == 0)
    g, X = mp5_metric()
    report('Myers-Perry (g_rr = rho^2/(r^2+a^2-mu)) is Ricci flat (exact zero at random rational points)',
           ricci_at_random_points(g, X, [mu, a]) == 0)
    g, X = mp5_metric(post_grr=True)
    report('Myers-Perry without the factor r^2 in g_rr (rho^2/(r^2(r^2+a^2-mu))) is NOT Ricci flat',
           ricci_at_random_points(g, X, [mu, a], npts=1) != 0)
    if '--slow' in sys.argv:
        g, X = kerr_metric()
        report('Kerr is Ricci flat (symbolic)', Geometry(g, X, simp=sp.cancel).ricci() == sp.zeros(4, 4))
        g, X = mp5_metric()
        report('Myers-Perry is Ricci flat (symbolic)', Geometry(g, X, simp=sp.cancel).ricci() == sp.zeros(5, 5))

    if '--slow' not in sys.argv:
        print('skipping the geodesic-monism E_ab check (about 20 minutes); GeodesicMonismAction.lean covers it')
        print('time %.1fs' % (time.time() - t0))
        sys.exit(0 if ok else 1)
    simp = sp.cancel
    g, X, gi = gm_hydrogen(4)
    E = Geometry(g, X, simp=simp, ginv=gi).field_equations()
    report('GM hydrogen (J^2 P2/(6 r^4)) solves E_ab = 0', E == sp.zeros(5, 5))
    g, X, gi = gm_hydrogen(2)
    E = Geometry(g, X, simp=simp, ginv=gi).field_equations()
    report('GM hydrogen with J^2 P2/(6 r^2) does NOT solve E_ab = 0', E != sp.zeros(5, 5))
    g, X, gi = gm_hydrogen(4, coef12=sp.Rational(1, 13))
    E = Geometry(g, X, simp=simp, ginv=gi).field_equations()
    report('negative control: 1/13 instead of 1/12 fails', E != sp.zeros(5, 5))
    print('time %.1fs' % (time.time() - t0))
    sys.exit(0 if ok else 1)
