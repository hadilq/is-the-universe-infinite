"""Derivation of the rotating 1+4 geodesic-monism solution (reduced equations).

Null Kaluza ansatz with flat transverse space and A_phi = f(r) sin^2(theta):
  (ii)  Lap^ nabla^k F_ki = 0         ->  f'''' type Euler equation, f = J/r + K r + B r^2 + Q r^4
  (iii) Lap(Lap H - F^2/4) = J_k J^k/2 + F^{ik} nabla_i J_k
Solved here by Legendre components (P0, P2) and the radial Euler operators
L_l y = y'' + 2y'/r - l(l+1) y/r^2, with r^m ln r at resonances.  The result is
compared with redshift.rot_parts(); the full E_ab = 0 is proved in
GeodesicMonismRotating.lean (the reduced equations are only used to find H).
"""
import sympy as sp
from tensors import Geometry
import redshift as R

r = R.r
x, ph = sp.symbols('x phi', real=True)
J, K, B, Q = R.J, R.Kk, R.Bb, R.Qq
P2s = sp.Symbol('P2')


def reduced_sources(f):
    A = [0, 0, f*(1 - x**2)]
    h = sp.diag(1, r**2/(1 - x**2), r**2*(1 - x**2))
    G = Geometry(h, [r, x, ph], simp=sp.cancel)
    X = [r, x, ph]; hi = G.ginv; Gm = G.Gam
    F = sp.Matrix(3, 3, lambda i, j: sp.diff(A[j], X[i]) - sp.diff(A[i], X[j]))
    Fup = (hi*F*hi).applyfunc(sp.cancel)
    F2 = sp.cancel(sum(F[i, j]*Fup[i, j] for i in range(3) for j in range(3)))
    DF = [[[sp.diff(F[i, j], X[k]) - sum(Gm[e][k][i]*F[e, j] + Gm[e][k][j]*F[i, e] for e in range(3))
            for j in range(3)] for i in range(3)] for k in range(3)]
    Jv = [sp.cancel(sum(hi[k, l]*DF[l][k][i] for k in range(3) for l in range(3))) for i in range(3)]
    Jup = [sum(hi[i, j]*Jv[j] for j in range(3)) for i in range(3)]
    JJ = sp.cancel(sum(Jv[i]*Jup[i] for i in range(3)))
    DJ = [[sp.diff(Jv[k], X[i]) - sum(Gm[e][i][k]*Jv[e] for e in range(3)) for k in range(3)] for i in range(3)]
    FDJ = sp.cancel(sum(Fup[i, k]*DJ[i][k] for i in range(3) for k in range(3)))
    return Jv, sp.cancel(F2/4), sp.cancel(JJ/2 + FDJ)


def split(e):
    e = sp.expand(sp.expand(e).subs(x, sp.sqrt((1 + 2*P2s)/3)))
    return {0: sp.expand(e.coeff(P2s, 0)), 2: sp.expand(e.coeff(P2s, 1))}


def Linv(y, ell):
    out = 0
    for term in sp.Add.make_args(sp.expand(y)):
        c, n = term.as_coeff_exponent(r)
        m = n + 2
        k = m*(m + 1) - ell*(ell + 1)
        out += c*r**m/k if k != 0 else c*r**m*sp.log(r)/(2*m + 1)
    return sp.expand(out)


def L(y, ell):
    return sp.expand(sp.diff(y, r, 2) + 2*sp.diff(y, r)/r - ell*(ell + 1)*y/r**2)


if __name__ == '__main__':
    import sys
    f = J/r + K*r + B*r**2 + Q*r**4
    Jv, F4, src = reduced_sources(f)
    print('current J_phi =', sp.factor(Jv[2]))
    # (ii): J_phi = g(r) sin^2 with g'' - 2g/r^2 = 0
    g = sp.cancel(Jv[2]/(1 - x**2))
    ok = sp.simplify(sp.diff(g, r, 2) - 2*g/r**2) == 0
    print(('PASS' if ok else 'FAIL') + ' (ii) fourth-order Maxwell equation for all four modes')
    S, F4s = split(src), split(F4)
    f_, h0, h2, Hrad = R.rot_parts()
    target = {0: h0, 2: h2}
    for ell in (0, 2):
        W = Linv(S[ell], ell)
        H = Linv(sp.expand(W + F4s[ell]), ell)
        good = sp.simplify(L(L(H, ell) - F4s[ell], ell) - S[ell]) == 0 and sp.simplify(H - target[ell]) == 0
        ok &= good
        print(('PASS' if good else 'FAIL') + f' (iii) l = {ell}: h_{ell} as in the post and in Lean')
    print('PASS rotating solution derivation' if ok else 'FAIL rotating solution derivation')
    sys.exit(0 if ok else 1)
