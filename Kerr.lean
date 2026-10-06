import Patch

/-
# Kerr: redshift coefficients in the equator and pole patches

Boyer–Lindquist, G = c = 1, clock ξ = ∂_t, ρ² = r² + a² cos²θ, Δ = r² − 2Mr + a²:
    g_tt = −(1 − 2Mr/ρ²),   g_tφ = −2Mar sin²θ/ρ²,
    g_φφ = (r² + a² + 2Ma²r sin²θ/ρ²) sin²θ.
Generators: r (Laurent), c = cos θ, M, a.

Results (S = r³ + a²r + 2Ma²), with z = C_r δr + C_t (ϑΔθ − Δθ²/2) + n_φ(D_r δr + D_t(…)):

  equator, static: C_r = M/(r(r − 2M)),           C_t = 2Ma²/(r²(r − 2M))
  equator, ZAMO:   C_r = M(r⁴ + 2a²r² + a⁴ − 4Ma²r)/(r S Δ),   C_t = 2Ma²(r² + a²)/(r² S)
                   D_r = P ∂_r w = −2Ma(3r² + a²)/(r √Δ S),   D_t = −4Ma³ √Δ/(r² S)
                   with P = S/(r√Δ), ∂_r w = −2Ma(3r² + a²)/S², ∂²_θ w = −4Ma³Δ/(r S²)
  pole, static:    C_r = M(r² − a²)/((r² + a²)Δ),   C_t = −2Ma²r/((r² + a²)Δ)
  pole, ZAMO:      C_r as static,                   C_t = −2Ma²r/(r² + a²)²,   D = 0
-/

namespace GMA.KerrCert
open GMA.Patch

def iR := 0
def iC := 1
def iM := 2
def iA := 3
def A : Alg := ⟨4, fun m a => [(m, a)]⟩
def v (i : Nat) (e : Int := 1) : Poly := varP A i e
def k (n : Int) : Poly := constP A (Q.ofInt n)
def m := mulP A
def F (n d : Poly) : Frac := ⟨n, d⟩
def P (n : Poly) : Frac := Frac.ofP A n

def r := v iR
def c := v iC
def M := v iM
def a := v iA
def r2 := m r r
def a2 := m a a
def S2 := subP (k 1) (m c c)
def rho2 := addP r2 (m a2 (m c c))
def Delta := addP (subP r2 (m (k 2) (m M r))) a2
def S := sumP [m r2 r, m a2 r, m (k 2) (m M a2)]

def data : Data :=
  { A := A, iR := iR, iC := iC
    X := F (negP (subP rho2 (m (k 2) (m M r)))) rho2
    Y := F (negP (m (k 2) (m M (m a (m r S2))))) rho2
    Z := F (m (addP (m (addP r2 a2) rho2) (m (k 2) (m M (m a2 (m r S2))))) S2) rho2
    -- N² = Δρ²/((r²+a²)ρ² + 2Ma²r sin²θ),  w = 2Mar/((r²+a²)ρ² + 2Ma²r sin²θ)
    N2r := F (m Delta rho2) (addP (m (addP r2 a2) rho2) (m (k 2) (m M (m a2 (m r S2)))))
    wr := F (m (k 2) (m M (m a r))) (addP (m (addP r2 a2) rho2) (m (k 2) (m M (m a2 (m r S2))))) }

def equatorClaim : Claim :=
  { staticCr := F M (m r (subP r (m (k 2) M)))
    staticCt := F (m (k 2) (m M a2)) (m r2 (subP r (m (k 2) M)))
    zamoCr := F (m M (sumP [m r2 r2, m (k 2) (m a2 r2), m a2 a2, negP (m (k 4) (m M (m a2 r)))]))
                (m r (m S Delta))
    zamoCt := F (m (k 2) (m M (m a2 (addP r2 a2)))) (m r2 S) }

def equatorExtra : EquatorExtra :=
  { P2 := F (m S S) (m r2 Delta)
    wr := F (negP (m (k 2) (m M (m a (addP (m (k 3) r2) a2))))) (m S S)
    wtt := F (negP (m (k 4) (m M (m a2 (m a Delta))))) (m r (m S S)) }

def poleClaim : Claim :=
  { staticCr := F (m M (subP r2 a2)) (m (addP r2 a2) Delta)
    staticCt := F (negP (m (k 2) (m M (m a2 r)))) (m (addP r2 a2) Delta)
    zamoCr := F (m M (subP r2 a2)) (m (addP r2 a2) Delta)
    zamoCt := F (negP (m (k 2) (m M (m a2 r)))) (m (addP r2 a2) (addP r2 a2)) }

end GMA.KerrCert

open GMA GMA.Patch GMA.KerrCert in
theorem kerr_equator : checkClaim data .equator equatorClaim = true := by native_decide
open GMA GMA.Patch GMA.KerrCert in
theorem kerr_equator_zamo_cross : checkEquatorExtra data equatorExtra = true := by native_decide
open GMA GMA.Patch GMA.KerrCert in
theorem kerr_pole : checkClaim data .pole poleClaim = true := by native_decide
open GMA GMA.Patch GMA.KerrCert in
theorem kerr_pole_cross_term_vanishes : poleZvanishes data = true := by native_decide

open GMA GMA.Patch GMA.KerrCert in
/-- Negative control: the weak-field value C_r = M/r² is not the exact
    coefficient, so the checker rejects it. -/
theorem kerr_control_weak_field :
    checkClaim data .equator { equatorClaim with staticCr := F M r2 } = false := by native_decide

open GMA GMA.Patch GMA.KerrCert in
/-- Negative control: the pole coefficients are not the equator ones. -/
theorem kerr_control_patches_differ : checkClaim data .equator poleClaim = false := by native_decide
