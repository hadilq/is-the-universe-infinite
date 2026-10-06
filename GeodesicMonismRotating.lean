import Patch

/-
# Geodesic monism in 1 + 4 dimensions: the general rotating (dipolar) solution

Null Kaluza form in five dimensions, coordinates (t, u, r, θ, φ):
    ds² = 2 dt du + (1 + 2H) du² + 2 A_φ dφ du + dr² + r² dθ² + r² sin²θ dφ²,
with the rotation carried by the dipolar potential
    A_φ = f(r) sin²θ,      f = J/r + K r + B r² + Q r⁴ .
The four modes are the four solutions of the fourth-order Maxwell equation
□̂∇̂ᵏF_kφ = 0 for an ℓ = 1 potential: J (dipole), B (uniform), and the two
partner modes K, Q that carry the geometric current J_φ = 2(5Qr³ − K) sin²θ / r.
The back-reaction gives
    H = c₀/r + c₁ + c₂ r + c₃ r² + h₀(r) + h₂(r) P₂(cos θ),
    h₀ = J²/(12r⁴) + JK/(6r²) − 2JQr/3 + K² ln r + 4BKr/3 + B²r²/3 + KQr³
         + 2BQr⁴/3 + 11Q²r⁶/42,
    h₂ = J²/(6r⁴) − 2JK/(3r²) − 2JB/(3r) + JQr/3 − K²/2 − 2BKr/3 − 2KQr³/3
         − 4BQr⁴/21 − 29Q²r⁶/126.
J alone gives back the hydrogen solution.

Lean computes all 25 components of E_ab = □R_ab + ½g_ab□R − ∇_a∇_bR
+ 2R_acbdR^cd − ½g_abR_cdR^cd and proves they vanish identically, for all
ten constants.  The ring is that of GeodesicMonismAction.lean with one more
generator, L = ln r, with ∂_r L = r⁻¹; ln r is transcendental over the
rational functions, so the normal form stays unique.  Negative controls:
11/42 → 11/43, and dropping the K² ln r term, both break E_ab = 0.
-/

namespace GMA.Rotating

/-- generators: r, s = sin θ, c = cos θ, L = ln r, J, K, B, Q, c₀, c₁, c₂, c₃ -/
def NV := 12
def iR := 0
def iS := 1
def iC := 2
def iL := 3
def iJ := 4
def iK := 5
def iB := 6
def iQ := 7
def iC0 := 8

def reduce (m : Mono) (a : Q) : Poly :=
  let rec go : Nat → Mono → Q → Poly
    | 0, m, a => [(m, a)]
    | f + 1, m, a =>
      let k := m.getD iC 0
      if k < 2 then [(m, a)]
      else
        let m0 := m.set iC (k - 2)
        let m1 := m0.set iS (m0.getD iS 0 + 2)
        addP (go f m0 a) (go f m1 (-a))
  go 8 m a

def A : Alg := ⟨NV, reduce⟩
def one : Poly := constP A (Q.ofInt 1)
def k (q : Q) : Poly := constP A q
def kz (n : Int) : Poly := constP A (Q.ofInt n)
def v (i : Nat) (e : Int := 1) : Poly := varP A i e
def mul := mulP A
def r (e : Int := 1) := v iR e
def s (e : Int := 1) := v iS e
def c := v iC
def L := v iL
def J := v iJ
def K := v iK
def B := v iB
def Qq := v iQ

/-- ∂_r (r ↦ 1, ln r ↦ r⁻¹) and ∂_θ (sin ↦ cos, cos ↦ −sin) on generators -/
def dGen (coord : Nat) (i : Nat) : Poly :=
  if coord = 2 then (if i = iR then one else if i = iL then r (-1) else [])
  else if coord = 3 then (if i = iS then c else if i = iC then negP (s) else [])
  else []

def d (coord : Nat) (p : Poly) : Poly :=
  if coord = 2 ∨ coord = 3 then derivP A (dGen coord) p else []

def term (q : Q) (factors : List Poly) : Poly :=
  scaleP q (factors.foldl mul one)

def f : Poly := sumP [mul J (r (-1)), mul K (r), mul B (r 2), mul Qq (r 4)]

def P2 : Poly := addP (scaleP (Q.frac 3 2) (mul c c)) (k (Q.frac (-1) 2))

/-- h₀ with the coefficient of Q²r⁶ and of K² ln r left as parameters (for the controls) -/
def h0 (qq : Q) (kl : Q) : Poly := sumP
  [ term (Q.frac 1 12) [J, J, r (-4)], term (Q.frac 1 6) [J, K, r (-2)], term (Q.frac (-2) 3) [J, Qq, r]
  , term kl [K, K, L], term (Q.frac 4 3) [B, K, r], term (Q.frac 1 3) [B, B, r 2]
  , term (Q.ofInt 1) [K, Qq, r 3], term (Q.frac 2 3) [B, Qq, r 4], term qq [Qq, Qq, r 6] ]

def h2 : Poly := sumP
  [ term (Q.frac 1 6) [J, J, r (-4)], term (Q.frac (-2) 3) [J, K, r (-2)], term (Q.frac (-2) 3) [J, B, r (-1)]
  , term (Q.frac 1 3) [J, Qq, r], term (Q.frac (-1) 2) [K, K], term (Q.frac (-2) 3) [B, K, r]
  , term (Q.frac (-2) 3) [K, Qq, r 3], term (Q.frac (-4) 21) [B, Qq, r 4]
  , term (Q.frac (-29) 126) [Qq, Qq, r 6] ]

def Hrad : Poly := sumP [mul (v iC0) (r (-1)), v (iC0 + 1), mul (v (iC0 + 2)) (r), mul (v (iC0 + 3)) (r 2)]

def H (qq kl : Q) : Poly := sumP [Hrad, h0 qq kl, mul h2 P2]

def geo (qq kl : Q) : Geo :=
  let Hh := H qq kl
  let guu := addP one (scaleP (Q.ofInt 2) Hh)
  let Aphi := mul f (s 2)
  let g : Nat → Nat → Poly := fun a b =>
    match a, b with
    | 0, 1 => one | 1, 0 => one
    | 1, 1 => guu
    | 1, 4 => Aphi | 4, 1 => Aphi
    | 2, 2 => one
    | 3, 3 => r 2
    | 4, 4 => mul (r 2) (s 2)
    | _, _ => []
  -- g^tu = 1, g^tt = −(1+2H) + f² sin²θ / r², g^tφ = −f/r², spatial block inverse of the flat metric
  let gtt := addP (negP guu) (mul (mul f f) (mul (s 2) (r (-2))))
  let gtp := negP (mul f (r (-2)))
  let gi : Nat → Nat → Poly := fun a b =>
    match a, b with
    | 0, 1 => one | 1, 0 => one
    | 0, 0 => gtt
    | 0, 4 => gtp | 4, 0 => gtp
    | 2, 2 => one
    | 3, 3 => r (-2)
    | 4, 4 => mul (r (-2)) (s (-2))
    | _, _ => []
  ⟨A, 5, g, gi, d⟩

def rotating : Geo := geo (Q.frac 11 42) (Q.ofInt 1)

/-- the geometric current J_φ = ∇̂ᵏF_kφ = 2(5Qr³ − K) sin²θ / r is non-zero for K, Q ≠ 0:
    check of ∂_r² A_φ + (sin²θ/r²) ∂²_{cos θ}A_φ written as r-part f″ − 2f/r² -/
def current : Poly := subP (d 2 (d 2 f)) (scaleP (Q.ofInt 2) (mul f (r (-2))))
def currentClaim : Poly := subP (scaleP (Q.ofInt 10) (mul Qq (r 2))) (scaleP (Q.ofInt 2) (mul K (r (-1))))

end GMA.Rotating

open GMA GMA.Rotating in
theorem rotating_inverse : rotating.inverseOk = true := by native_decide

open GMA GMA.Rotating in
/-- **Theorem (the rotating 1 + 4 solution).**  All 25 components of E_ab
    vanish identically, for all c₀, c₁, c₂, c₃, J, K, B, Q. -/
theorem rotating_solves : rotating.fieldEq.all isZeroP = true := by native_decide

open GMA GMA.Rotating in
/-- the partner modes carry a current (so this is not a Maxwell-vacuum solution) -/
theorem rotating_current : isZeroP (subP current currentClaim) = true := by native_decide

open GMA GMA.Rotating in
theorem rotating_control_Q2 : (geo (Q.frac 11 43) (Q.ofInt 1)).fieldEq.all isZeroP = false := by native_decide

open GMA GMA.Rotating in
theorem rotating_control_log : (geo (Q.frac 11 42) (Q.ofInt 0)).fieldEq.all isZeroP = false := by native_decide

/-! ## Redshift coefficients of the rotating solution, both patches

Clock ξ = ∂_t − ∂_u:  g(ξ,ξ) = −(1 − 2H),  g(ξ,∂_φ) = −f sin²θ,  g_φφ = r² sin²θ,
N² = 1 − 2H + f² sin²θ/r²,  ω = f/r².  On the equator H_e = H_rad + h₀ − h₂/2,
at the pole H_p = H_rad + h₀ + h₂, and N_e² = 1 − 2H_e + f²/r².

  equator, static: C_r = −H_e′/(1 − 2H_e),        C_θ = −3h₂/(1 − 2H_e)
  equator, ZAMO:   C_r = (N_e²)′/(2N_e²),          C_θ = −(3h₂ + f²/r²)/N_e²,
                   D_r = (r/N_e)(f/r²)′,           D_θ = 0   (ω does not depend on θ)
  pole, static:    C_r = −H_p′/(1 − 2H_p),        C_θ = 3h₂/(1 − 2H_p)
  pole, ZAMO:      C_r as static,                  C_θ = (3h₂ + f²/r²)/(1 − 2H_p),  D = 0

Generators: r (Laurent), c = cos θ, L = ln r, J, K, B, Q, c₀, c₁, c₂, c₃. -/

namespace GMA.RotatingCert
open GMA.Patch

def A : Alg := ⟨11, fun m a => [(m, a)]⟩
def v (i : Nat) (e : Int := 1) : Poly := varP A i e
def kq (q : Q) : Poly := constP A q
def kz (n : Int) : Poly := constP A (Q.ofInt n)
def m := mulP A
def F (n d : Poly) : Frac := ⟨n, d⟩
def P (n : Poly) : Frac := Frac.ofP A n
def r (e : Int := 1) := v 0 e
def c := v 1
def L := v 2
def J := v 3
def K := v 4
def B := v 5
def Qq := v 6
def term (q : Q) (fs : List Poly) : Poly := scaleP q (fs.foldl m (kz 1))

def f : Poly := sumP [m J (r (-1)), m K (r), m B (r 2), m Qq (r 4)]
def h0 : Poly := sumP
  [ term (Q.frac 1 12) [J, J, r (-4)], term (Q.frac 1 6) [J, K, r (-2)], term (Q.frac (-2) 3) [J, Qq, r]
  , term (Q.ofInt 1) [K, K, L], term (Q.frac 4 3) [B, K, r], term (Q.frac 1 3) [B, B, r 2]
  , term (Q.ofInt 1) [K, Qq, r 3], term (Q.frac 2 3) [B, Qq, r 4], term (Q.frac 11 42) [Qq, Qq, r 6] ]
def h2 : Poly := sumP
  [ term (Q.frac 1 6) [J, J, r (-4)], term (Q.frac (-2) 3) [J, K, r (-2)], term (Q.frac (-2) 3) [J, B, r (-1)]
  , term (Q.frac 1 3) [J, Qq, r], term (Q.frac (-1) 2) [K, K], term (Q.frac (-2) 3) [B, K, r]
  , term (Q.frac (-2) 3) [K, Qq, r 3], term (Q.frac (-4) 21) [B, Qq, r 4]
  , term (Q.frac (-29) 126) [Qq, Qq, r 6] ]
def Hrad : Poly := sumP [m (v 7) (r (-1)), v 8, m (v 9) (r), m (v 10) (r 2)]
def P2 : Poly := addP (scaleP (Q.frac 3 2) (m c c)) (kq (Q.frac (-1) 2))
def H : Poly := sumP [Hrad, h0, m h2 P2]
def S2 : Poly := subP (kz 1) (m c c)
def f2r : Poly := m (m f f) (r (-2))

def data : Data :=
  { A := A, iR := 0, iC := 1, iL := 2
    X := P (subP (scaleP (Q.ofInt 2) H) (kz 1))
    Y := P (negP (m f S2))
    Z := P (m (r 2) S2)
    N2r := P (addP (subP (kz 1) (scaleP (Q.ofInt 2) H)) (m f2r S2))
    wr := P (m f (r (-2))) }

def dR' := dR data
def He := sumP [Hrad, h0, scaleP (Q.frac (-1) 2) h2]
def Hp := sumP [Hrad, h0, h2]
def ae := subP (kz 1) (scaleP (Q.ofInt 2) He)
def ap := subP (kz 1) (scaleP (Q.ofInt 2) Hp)
def Ne2 := addP ae f2r

def equatorClaim : Claim :=
  { staticCr := F (negP (dR' He)) ae
    staticCt := F (scaleP (Q.ofInt (-3)) h2) ae
    zamoCr := F (dR' Ne2) (scaleP (Q.ofInt 2) Ne2)
    zamoCt := F (negP (addP (scaleP (Q.ofInt 3) h2) f2r)) Ne2 }

def equatorExtra : EquatorExtra :=
  { P2 := F (r 2) Ne2
    wr := P (dR' (m f (r (-2))))
    wtt := F [] (kz 1) }

def poleClaim : Claim :=
  { staticCr := F (negP (dR' Hp)) ap
    staticCt := F (scaleP (Q.ofInt 3) h2) ap
    zamoCr := F (negP (dR' Hp)) ap
    zamoCt := F (addP (scaleP (Q.ofInt 3) h2) f2r) ap }

end GMA.RotatingCert

open GMA GMA.Patch GMA.RotatingCert in
theorem rotating_equator : checkClaim data .equator equatorClaim = true := by native_decide
open GMA GMA.Patch GMA.RotatingCert in
theorem rotating_equator_zamo_cross : checkEquatorExtra data equatorExtra = true := by native_decide
open GMA GMA.Patch GMA.RotatingCert in
theorem rotating_pole : checkClaim data .pole poleClaim = true := by native_decide
open GMA GMA.Patch GMA.RotatingCert in
theorem rotating_pole_cross_term_vanishes : poleZvanishes data = true := by native_decide
