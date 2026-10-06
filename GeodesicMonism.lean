import Patch

/-
# Geodesic monism (null Kaluza hydrogen): redshift coefficients, both patches

    ds² = 2 dt du + (1 + 2H) du² + 2 A_φ dφ du + dr² + r² dθ² + r² sin²θ dφ²,
    H = c₀/r + c₁ + c₂ r + c₃ r² + J²/(12 r⁴) + J² P₂(cos θ)/(6 r⁴),   A_φ = J sin²θ/r,
which solves E_ab = 0 (GeodesicMonismAction.lean).  Clock ξ = ∂_t − ∂_u, so
    g(ξ,ξ) = −(1 − 2H),   g(ξ,∂_φ) = −A_φ,   g_φφ = r² sin²θ,
    N² = 1 − 2H + J² sin²θ/r⁴,   w = J/r³.
The two dipole terms combine: J²/12 + J²P₂/6 = J² cos²θ/4, so
    H = H₀(r) + J² cos²θ/(4 r⁴),    H₀ = c₀/r + c₁ + c₂ r + c₃ r²,
and on the equator H = H₀ carries no J at all.  With a₀ = 1 − 2H₀,
N_e² = a₀ + J²/r⁴, h_p = H₀ + J²/(4r⁴), a_p = 1 − 2h_p:

  equator, static: C_r = −H₀′/a₀,                    C_t = −J²/(2r⁴a₀)
  equator, ZAMO:   C_r = (−H₀′ − 2J²/r⁵)/N_e²,        C_t = −3J²/(2r⁴N_e²)
                   P = r/N_e, ∂_r w = −3J/r⁴, ∂²_θ w = 0, so D_r = −3J/(r³N_e), D_t = 0
  pole, static:    C_r = −h_p′/a_p,                   C_t = J²/(2r⁴a_p)
  pole, ZAMO:      C_r as static,                     C_t = 3J²/(2r⁴a_p),   D = 0
Generators: r (Laurent), c = cos θ, J, c₀, c₁, c₂, c₃.
-/

namespace GMA.GMCert
open GMA.Patch

def iR := 0
def iC := 1
def iJ := 2
def A : Alg := ⟨7, fun m a => [(m, a)]⟩
def v (i : Nat) (e : Int := 1) : Poly := varP A i e
def k (n : Int) : Poly := constP A (Q.ofInt n)
def kq (q : Q) : Poly := constP A q
def m := mulP A
def F (n d : Poly) : Frac := ⟨n, d⟩
def P (n : Poly) : Frac := Frac.ofP A n

def r (e : Int := 1) := v iR e
def c := v iC
def J := v iJ
def J2 := m J J
def S2 := subP (k 1) (m c c)
def P2 := addP (scaleP (Q.frac 3 2) (m c c)) (kq (Q.frac (-1) 2))
def H0 := sumP [m (v 3) (r (-1)), v 4, m (v 5) (r), m (v 6) (r 2)]
def H0' := sumP [negP (m (v 3) (r (-2))), v 5, m (k 2) (m (v 6) (r))]
def H := sumP [H0, scaleP (Q.frac 1 12) (m J2 (r (-4))), scaleP (Q.frac 1 6) (m J2 (m P2 (r (-4))))]

def data : Data :=
  { A := A, iR := iR, iC := iC
    X := P (subP (m (k 2) H) (k 1))
    Y := P (negP (m J (m S2 (r (-1)))))
    Z := P (m (r 2) S2)
    N2r := P (addP (subP (k 1) (m (k 2) H)) (m J2 (m S2 (r (-4)))))
    wr := P (m J (r (-3))) }

def a0 := subP (k 1) (m (k 2) H0)
def Ne := addP a0 (m J2 (r (-4)))
def hp := addP H0 (scaleP (Q.frac 1 4) (m J2 (r (-4))))
def hp' := subP H0' (m J2 (r (-5)))
def ap := subP (k 1) (m (k 2) hp)

def equatorClaim : Claim :=
  { staticCr := F (negP H0') a0
    staticCt := F (negP J2) (m (k 2) (m (r 4) a0))
    zamoCr := F (subP (negP H0') (m (k 2) (m J2 (r (-5))))) Ne
    zamoCt := F (negP (m (k 3) J2)) (m (k 2) (m (r 4) Ne)) }

def equatorExtra : EquatorExtra :=
  { P2 := F (r 2) Ne
    wr := F (negP (m (k 3) (m J (r (-4))))) (k 1)
    wtt := F [] (k 1) }

def poleClaim : Claim :=
  { staticCr := F (negP hp') ap
    staticCt := F J2 (m (k 2) (m (r 4) ap))
    zamoCr := F (negP hp') ap
    zamoCt := F (m (k 3) J2) (m (k 2) (m (r 4) ap)) }

/-- H = H₀ + J² cos²θ/(4r⁴) -/
def Hsimplifies : Bool := isZeroP (subP H (addP H0 (scaleP (Q.frac 1 4) (m J2 (m (m c c) (r (-4)))))))

end GMA.GMCert

open GMA GMA.GMCert in
theorem gm_H_simplifies : Hsimplifies = true := by native_decide
open GMA GMA.Patch GMA.GMCert in
theorem gm_equator : checkClaim data .equator equatorClaim = true := by native_decide
open GMA GMA.Patch GMA.GMCert in
theorem gm_equator_zamo_cross : checkEquatorExtra data equatorExtra = true := by native_decide
open GMA GMA.Patch GMA.GMCert in
theorem gm_pole : checkClaim data .pole poleClaim = true := by native_decide
open GMA GMA.Patch GMA.GMCert in
theorem gm_pole_cross_term_vanishes : poleZvanishes data = true := by native_decide
