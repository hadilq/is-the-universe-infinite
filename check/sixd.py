"""Derivation of the rotating three-sphere of geodesic monism in 1 + 1 + 4 (reduced equations).

Hopf coordinates on the transverse R^4, theta in [0, pi/2]:
    h = dr^2 + r^2 (dtheta^2 + sin^2 dphi^2 + cos^2 dpsi^2),
rotation in the phi plane and the psi plane at the same rate (the Hopf fibre):
    A = F(r) (sin^2 dphi + cos^2 dpsi),   rotation direction chi = d_phi + d_psi,  |chi| = r.
Then
    current J = (F'' + F'/r - 4F/r^2)(sin^2 dphi + cos^2 dpsi),
    (ii)  rough Laplacian of J = 0  ->  F = J/r^2 + K + B r^2 + Q r^4,
    (iii) Lap(Lap H - F^2/4) = J.J/2 + F^{ik} nabla_i J_k,  Lap y = y'' + 3y'/r on H(r),
    H_rot = J^2/(6r^6) + JK/(3r^4) - K^2 ln r/(2r^2) + KQ r^2 ln r + BQ r^4 + 19 Q^2 r^6/48.
GeodesicMonism6D.lean proves E_ab = 0 for the same metric written in Euler angles,
theta_E = 2 theta, phi_E = psi - phi, psi_E = phi + psi, where
sin^2 dphi + cos^2 dpsi = sigma_3/2, so its f is F/2 (and its J, K, B, Q are half of these).
This script also checks that coordinate map.
"""
import sys
import sympy as sp
from tensors import Geometry

r, th = sp.symbols('r theta', positive=True)
ph, ps = sp.symbols('phi psi')
J, K, B, Q = sp.symbols('J K B Q', real=True)


def main():
    ok = True
    # the coordinate map to the Euler angles used in Lean
    a, b, dt_ = sp.symbols('a b dth')
    thE = 2*th
    euler = (r**2/4)*((2*dt_)**2 + sp.sin(thE)**2*(b - a)**2 + ((a + b) + sp.cos(thE)*(b - a))**2)
    hopf = r**2*(dt_**2 + sp.sin(th)**2*a**2 + sp.cos(th)**2*b**2)
    good = sp.simplify(sp.expand_trig(euler - hopf)) == 0
    sig3 = sp.expand_trig((a + b) + sp.cos(thE)*(b - a))
    good &= sp.simplify(sig3 - 2*(sp.sin(th)**2*a + sp.cos(th)**2*b)) == 0
    print(('PASS' if good else 'FAIL') + ' Euler (theta_E = 2theta, phi_E = psi - phi, psi_E = phi + psi) = Hopf;'
          ' sigma_3 = 2(sin^2 dphi + cos^2 dpsi)')
    ok &= good
    X = [r, th, ph, ps]
    h = sp.diag(1, r**2, r**2*sp.sin(th)**2, r**2*sp.cos(th)**2)
    G = Geometry(h, X, simp=sp.simplify)
    good = G.ricci() == sp.zeros(4, 4)
    print(('PASS' if good else 'FAIL') + ' the transverse space is flat')
    ok &= good
    F_ = J/r**2 + K + B*r**2 + Q*r**4
    A = [0, 0, F_*sp.sin(th)**2, F_*sp.cos(th)**2]
    hi, Gm = G.ginv, G.Gam
    Fm = sp.Matrix(4, 4, lambda i, j: sp.diff(A[j], X[i]) - sp.diff(A[i], X[j]))
    Fup = (hi*Fm*hi).applyfunc(sp.simplify)
    F4 = sp.expand(sp.simplify(sum(Fm[i, j]*Fup[i, j] for i in range(4) for j in range(4))/4))
    DF = [[[sp.diff(Fm[i, j], X[k]) - sum(Gm[e][k][i]*Fm[e, j] + Gm[e][k][j]*Fm[i, e] for e in range(4))
            for j in range(4)] for i in range(4)] for k in range(4)]
    Jv = [sp.simplify(sum(hi[k, l]*DF[l][k][i] for k in range(4) for l in range(4))) for i in range(4)]
    g = sp.diff(F_, r, 2) + sp.diff(F_, r)/r - 4*F_/r**2
    good = all(sp.simplify(Jv[i] - g*[0, 0, sp.sin(th)**2, sp.cos(th)**2][i]) == 0 for i in range(4))
    print(('PASS' if good else 'FAIL') + f' current J = g(r)(sin^2 dphi + cos^2 dpsi), g = {sp.factor(sp.simplify(g))}')
    ok &= good
    Dw = [[sp.diff(Jv[k], X[i]) - sum(Gm[e][i][k]*Jv[e] for e in range(4)) for k in range(4)] for i in range(4)]
    DDw = [[[sp.diff(Dw[i][k], X[j]) - sum(Gm[e][j][i]*Dw[e][k] + Gm[e][j][k]*Dw[i][e] for e in range(4))
             for k in range(4)] for i in range(4)] for j in range(4)]
    lapJ = [sp.simplify(sum(hi[i, j]*DDw[j][i][k] for i in range(4) for j in range(4))) for k in range(4)]
    good = all(v == 0 for v in lapJ)
    print(('PASS' if good else 'FAIL') + ' (ii) the fourth-order Maxwell equation holds for all four modes')
    ok &= good
    Jup = [sum(hi[i, j]*Jv[j] for j in range(4)) for i in range(4)]
    JJ = sum(Jv[i]*Jup[i] for i in range(4))
    DJ = [[sp.diff(Jv[k], X[i]) - sum(Gm[e][i][k]*Jv[e] for e in range(4)) for k in range(4)] for i in range(4)]
    FDJ = sum(Fup[i, k]*DJ[i][k] for i in range(4) for k in range(4))
    src = sp.expand(sp.simplify(JJ/2 + FDJ))
    lg = sp.log(r)
    H = (J**2/(6*r**6) + J*K/(3*r**4) - K**2*lg/(2*r**2) + K*Q*r**2*lg + B*Q*r**4
         + sp.Rational(19, 48)*Q**2*r**6)
    L4 = lambda y: sp.diff(y, r, 2) + 3*sp.diff(y, r)/r
    good = sp.simplify(L4(L4(H) - F4) - src) == 0
    print(('PASS' if good else 'FAIL') + ' (iii) back-reaction with H_rot as in the post')
    ok &= good
    # A.A = F^2/r^2 and the Lean normalisation f = F/2
    AA = sp.simplify(sum(A[i]*hi[i, j]*A[j] for i in range(4) for j in range(4)))
    good = sp.simplify(AA - F_**2/r**2) == 0
    fE = F_/2
    HE = (sp.Rational(2, 3)*(J/2)**2/r**6 + sp.Rational(4, 3)*(J/2)*(K/2)/r**4 - 2*(K/2)**2*lg/r**2
          + 4*(K/2)*(Q/2)*r**2*lg + 4*(B/2)*(Q/2)*r**4 + sp.Rational(19, 12)*(Q/2)**2*r**6)
    good &= sp.simplify(HE - H) == 0
    print(('PASS' if good else 'FAIL') + ' |A|^2 = F^2/r^2, and with f = F/2 the Lean H is the same H')
    ok &= good
    print('PASS 1+1+4 derivation' if ok else 'FAIL 1+1+4 derivation')
    return ok


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
