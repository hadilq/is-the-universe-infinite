import Patch

/-
# Five-dimensional Myers–Perry, one spin: redshift coefficients, both patches

    ds² = −dt² + (μ/ρ²)(dt − a sin²θ dφ)² + (ρ²/(r² + a² − μ)) dr² + ρ² dθ²
          + (r² + a²) sin²θ dφ² + r² cos²θ dψ²,          ρ² = r² + a² cos²θ.
(The radial component is r²ρ²/Δ with Δ = r²(r² + a² − μ); without the factor r²
the metric is not Ricci flat — see check/solutions.py.)
  equator, ZAMO:   C_r = μ(r⁴ + 2a²r² + a⁴ − a²μ)/(r(r² + a² − μ)T),  C_t = a²μ(r² + a²)/(r²T)
                   P = T/(r²√(r² + a² − μ)),  ∂_r w = −2aμr(2r² + a²)/T²,
                   ∂²_θ w = −2a³μ(r² + a² − μ)/T²
                   so D_r = −2aμ(2r² + a²)/(r √(r² + a² − μ) T),
                      D_t = −2a³μ √(r² + a² − μ)/(r² T)
  pole, static:    C_r = μr/((r² + a²)(r² + a² − μ)),  C_t = −a²μ/((r² + a²)(r² + a² − μ))
  pole, ZAMO:      C_r as static,  C_t = −a²μ/(r² + a²)²,  D = 0
-/

namespace GMA.MPCert
open GMA.Patch

def iR := 0
def iC := 1
def iMu := 2
def iA := 3
def A : Alg := ⟨4, fun m a => [(m, a)]⟩
def v (i : Nat) (e : Int := 1) : Poly := varP A i e
def k (n : Int) : Poly := constP A (Q.ofInt n)
def m := mulP A
def F (n d : Poly) : Frac := ⟨n, d⟩

def r := v iR
def c := v iC
def mu := v iMu
def a := v iA
def r2 := m r r
def a2 := m a a
def S2 := subP (k 1) (m c c)
def rho2 := addP r2 (m a2 (m c c))
def B := addP (m (addP r2 a2) rho2) (m mu (m a2 S2))     -- ρ² g_φφ / sin²θ
def Dl := subP (addP r2 a2) mu                            -- r² + a² − μ
def T := sumP [m r2 r2, m a2 r2, m a2 mu]

def data : Data :=
  { A := A, iR := iR, iC := iC
    X := F (negP (subP rho2 mu)) rho2
    Y := F (negP (m mu (m a S2))) rho2
    Z := F (m S2 B) rho2
    N2r := F (m Dl rho2) B
    wr := F (m mu a) B }

def equatorClaim : Claim :=
  { staticCr := F mu (m r (subP r2 mu))
    staticCt := F (m a2 mu) (m r2 (subP r2 mu))
    zamoCr := F (m mu (sumP [m r2 r2, m (k 2) (m a2 r2), m a2 a2, negP (m a2 mu)])) (m r (m Dl T))
    zamoCt := F (m a2 (m mu (addP r2 a2))) (m r2 T) }

def equatorExtra : EquatorExtra :=
  { P2 := F (m T T) (m (m r2 r2) Dl)
    wr := F (negP (m (k 2) (m a (m mu (m r (addP (m (k 2) r2) a2)))))) (m T T)
    wtt := F (negP (m (k 2) (m a2 (m a (m mu Dl))))) (m T T) }

def poleClaim : Claim :=
  { staticCr := F (m mu r) (m (addP r2 a2) Dl)
    staticCt := F (negP (m a2 mu)) (m (addP r2 a2) Dl)
    zamoCr := F (m mu r) (m (addP r2 a2) Dl)
    zamoCt := F (negP (m a2 mu)) (m (addP r2 a2) (addP r2 a2)) }

end GMA.MPCert

open GMA GMA.Patch GMA.MPCert in
theorem mp_equator : checkClaim data .equator equatorClaim = true := by native_decide
open GMA GMA.Patch GMA.MPCert in
theorem mp_equator_zamo_cross : checkEquatorExtra data equatorExtra = true := by native_decide
open GMA GMA.Patch GMA.MPCert in
theorem mp_pole : checkClaim data .pole poleClaim = true := by native_decide
open GMA GMA.Patch GMA.MPCert in
theorem mp_pole_cross_term_vanishes : poleZvanishes data = true := by native_decide
