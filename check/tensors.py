"""Minimal coordinate tensor calculus used by the checks.

Conventions (Misner-Thorne-Wheeler, signature -+++):
  Gamma^a_{bc} = 1/2 g^{ad} (d_b g_{dc} + d_c g_{db} - d_d g_{bc})
  R^a_{bcd}    = d_c Gamma^a_{db} - d_d Gamma^a_{cb}
                 + Gamma^a_{ce} Gamma^e_{db} - Gamma^a_{de} Gamma^e_{cb}
  R_{bd}       = R^a_{bad}
"""
import sympy as sp


class Geometry:
    def __init__(self, g, coords, simp=sp.simplify, ginv=None):
        self.g = sp.Matrix(g)
        self.x = list(coords)
        self.n = len(coords)
        self.simp = simp
        self.ginv = sp.Matrix(ginv) if ginv is not None else self.g.inv().applyfunc(simp)
        n = self.n
        self.Gam = [[[simp(sum(self.ginv[a, d] * (sp.diff(self.g[d, c], self.x[b])
                                                  + sp.diff(self.g[d, b], self.x[c])
                                                  - sp.diff(self.g[b, c], self.x[d]))
                                for d in range(n)) / 2)
                      for c in range(n)] for b in range(n)] for a in range(n)]

    # Riemann R^a_{bcd}
    def riemann(self):
        if hasattr(self, "_R"):
            return self._R
        n, G, x = self.n, self.Gam, self.x
        R = [[[[self.simp(sp.diff(G[a][d][b], x[c]) - sp.diff(G[a][c][b], x[d])
                          + sum(G[a][c][e] * G[e][d][b] - G[a][d][e] * G[e][c][b] for e in range(n)))
                for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
        self._R = R
        return R

    def ricci(self):
        """R_bd = d_a G^a_db - d_d G^a_ab + G^a_ae G^e_db - G^a_de G^e_ab (contracted directly)."""
        if hasattr(self, "_Ric"):
            return self._Ric
        n, G, x = self.n, self.Gam, self.x
        self._Ric = sp.Matrix(n, n, lambda b, d: self.simp(sum(
            sp.diff(G[a][d][b], x[a]) - sp.diff(G[a][a][b], x[d])
            + sum(G[a][a][e]*G[e][d][b] - G[a][d][e]*G[e][a][b] for e in range(n))
            for a in range(n))))
        return self._Ric

    def scalar(self):
        Ric = self.ricci()
        return self.simp(sum(self.ginv[a, b] * Ric[a, b] for a in range(self.n) for b in range(self.n)))

    # covariant derivative of a covariant 2-tensor T_{ab}: returns D[c][a][b] = nabla_c T_ab
    def cov2(self, T):
        n, G, x = self.n, self.Gam, self.x
        return [[[self.simp(sp.diff(T[a, b], x[c])
                            - sum(G[e][c][a] * T[e, b] + G[e][c][b] * T[a, e] for e in range(n)))
                  for b in range(n)] for a in range(n)] for c in range(n)]

    # nabla_d of a covariant 3-tensor S[c][a][b]
    def cov3(self, S):
        n, G, x = self.n, self.Gam, self.x
        return [[[[self.simp(sp.diff(S[c][a][b], x[d])
                             - sum(G[e][d][c] * S[e][a][b] + G[e][d][a] * S[c][e][b]
                                   + G[e][d][b] * S[c][a][e] for e in range(n)))
                   for b in range(n)] for a in range(n)] for c in range(n)] for d in range(n)]

    def field_equations(self):
        """E_ab of S = int R_ab R^ab sqrt(-g):
        E_ab = Box R_ab + 1/2 g_ab Box R - nabla_a nabla_b R
               + 2 R_acbd R^cd - 1/2 g_ab R_cd R^cd
        """
        n, gi, g, x = self.n, self.ginv, self.g, self.x
        Ric = self.ricci()
        Rs = self.scalar()
        DR = self.cov2(Ric)
        DDR = self.cov3(DR)  # DDR[d][c][a][b] = nabla_d nabla_c R_ab
        boxRic = sp.Matrix(n, n, lambda a, b: self.simp(sum(gi[c, d] * DDR[d][c][a][b]
                                                             for c in range(n) for d in range(n))))
        dRs = [sp.diff(Rs, xx) for xx in x]
        hessR = sp.Matrix(n, n, lambda a, b: self.simp(sp.diff(Rs, x[a], x[b])
                                                        - sum(self.Gam[e][a][b] * dRs[e] for e in range(n))))
        boxR = self.simp(sum(gi[a, b] * hessR[a, b] for a in range(n) for b in range(n)))
        Rup = (gi * Ric * gi).applyfunc(self.simp)
        R = self.riemann()
        # R_{acbd} = g_{ae} R^e_{cbd}
        Rlow = lambda a, c, b, d: sum(g[a, e] * R[e][c][b][d] for e in range(n))
        RR = sp.Matrix(n, n, lambda a, b: sum(Rlow(a, c, b, d) * Rup[c, d]
                                               for c in range(n) for d in range(n)))
        Ric2 = sum(Ric[a, b] * Rup[a, b] for a in range(n) for b in range(n))
        E = sp.Matrix(n, n, lambda a, b: self.simp(boxRic[a, b] + g[a, b] * boxR / 2 - hessR[a, b]
                                                    + 2 * RR[a, b] - g[a, b] * Ric2 / 2))
        return E
