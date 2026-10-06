import Patch

/-
# Geodesic monism in 1 + 1 + 4: a three-sphere rotating along φ

Coordinates (t, u, r, ψ, θ, φ); the transverse space is ℝ⁴ in hyperspherical
coordinates, r = const being a round three-sphere:
    ds² = 2 dt du + (1 + 2H) du² + 2 F(r) sin²ψ sin²θ dφ du
          + dr² + r² (dψ² + sin²ψ (dθ² + sin²θ dφ²)),
    F = J/r² + K + B r² + Q r⁴.
The rotation is along ∂_φ only, in the plane x₁x₂ of ℝ⁴; its axis is the
plane x₃x₄ (sin ψ sin θ = 0).  With μ = sin²ψ sin²θ and the ℓ = 2 harmonic
Y = 2μ − 1 of the three-sphere,
    H = c₀/r² + c₁ + c₂ r² + c₃ ln r + h₀(r) + h₂(r) Y,
    h₀ = J²/(12r⁶) + JK/(6r⁴) − K² ln r/(4r²) − JQ ln r + BK ln r + B²r²/4
         + KQ r² ln r/2 − KQ r²/8 + BQ r⁴/2 + 19Q²r⁶/96,
    h₂ = JK ln r/(3r⁴) + K² ln r/(4r²) + K²/(16r²) + BJ/(2r²) + BK/2
         + KQ r² ln r/6 + BQ r⁴/8 + 27Q²r⁶/160.
check/sphere.py derives h₀ and h₂ from the reduced equations.

Every component lies in
𝓡 = ℚ[consts][r^{±1}, ln r, a^{±1}, b, s^{±1}, c]/(a² + b² − 1, s² + c² − 1),
a = sin ψ, b = cos ψ, s = sin θ, c = cos θ.  Lean computes all 36 components of
E_ab and proves they vanish identically for all eight constants.  Controls:
19/96 → 19/97 and dropping −K² ln r/(4r²) from h₀ both break it.
-/

namespace GMA.Sphere

/-- generators: r, a = sin ψ, b = cos ψ, s = sin θ, c = cos θ, L = ln r, J, K, B, Q, c₀, c₁, c₂, c₃ -/
def NV := 14
def iR := 0
def iA := 1
def iBc := 2
def iS := 3
def iC := 4
def iL := 5
def iJ := 6
def iK := 7
def iB := 8
def iQ := 9
def iC0 := 10

/-- normal form: b² → 1 − a², c² → 1 − s² -/
def reduce (m : Mono) (q : Q) : Poly :=
  let rec go : Nat → Mono → Q → Poly
    | 0, m, q => [(m, q)]
    | f + 1, m, q =>
      let kb := m.getD iBc 0
      let kc := m.getD iC 0
      if kb ≥ 2 then
        let m0 := m.set iBc (kb - 2)
        let m1 := m0.set iA (m0.getD iA 0 + 2)
        addP (go f m0 q) (go f m1 (-q))
      else if kc ≥ 2 then
        let m0 := m.set iC (kc - 2)
        let m1 := m0.set iS (m0.getD iS 0 + 2)
        addP (go f m0 q) (go f m1 (-q))
      else [(m, q)]
  go 16 m q

def A : Alg := ⟨NV, reduce⟩
def one : Poly := constP A (Q.ofInt 1)
def kq (q : Q) : Poly := constP A q
def v (i : Nat) (e : Int := 1) : Poly := varP A i e
def mul := mulP A
def r (e : Int := 1) := v iR e
def sa (e : Int := 1) := v iA e
def cb := v iBc
def s (e : Int := 1) := v iS e
def c := v iC
def L := v iL
def J := v iJ
def K := v iK
def B := v iB
def Qq := v iQ

/-- coordinates (t, u, r, ψ, θ, φ) = (0, 1, 2, 3, 4, 5) -/
def dGen (coord : Nat) (i : Nat) : Poly :=
  if coord = 2 then (if i = iR then one else if i = iL then r (-1) else [])
  else if coord = 3 then (if i = iA then cb else if i = iBc then negP (sa) else [])
  else if coord = 4 then (if i = iS then c else if i = iC then negP (s) else [])
  else []

def d (coord : Nat) (p : Poly) : Poly :=
  if coord = 2 ∨ coord = 3 ∨ coord = 4 then derivP A (dGen coord) p else []

def term (q : Q) (fs : List Poly) : Poly := scaleP q (fs.foldl mul one)

def F : Poly := sumP [mul J (r (-2)), K, mul B (r 2), mul Qq (r 4)]
/-- μ = sin²ψ sin²θ and Y = 2μ − 1 -/
def mu : Poly := mul (sa 2) (s 2)
def Y : Poly := subP (scaleP (Q.ofInt 2) mu) one

def Hrad : Poly := sumP [mul (v iC0) (r (-2)), v (iC0 + 1), mul (v (iC0 + 2)) (r 2), mul (v (iC0 + 3)) L]

/-- h₀, with the Q²r⁶ coefficient and the K² ln r/r² coefficient free (controls) -/
def h0 (qq kl : Q) : Poly := sumP
  [ term (Q.frac 1 12) [J, J, r (-6)], term (Q.frac 1 6) [J, K, r (-4)]
  , term kl [K, K, L, r (-2)], term (Q.ofInt (-1)) [J, Qq, L], term (Q.ofInt 1) [B, K, L]
  , term (Q.frac 1 4) [B, B, r 2], term (Q.frac 1 2) [K, Qq, L, r 2], term (Q.frac (-1) 8) [K, Qq, r 2]
  , term (Q.frac 1 2) [B, Qq, r 4], term qq [Qq, Qq, r 6] ]

def h2 : Poly := sumP
  [ term (Q.frac 1 3) [J, K, L, r (-4)], term (Q.frac 1 4) [K, K, L, r (-2)], term (Q.frac 1 16) [K, K, r (-2)]
  , term (Q.frac 1 2) [B, J, r (-2)], term (Q.frac 1 2) [B, K], term (Q.frac 1 6) [K, Qq, L, r 2]
  , term (Q.frac 1 8) [B, Qq, r 4], term (Q.frac 27 160) [Qq, Qq, r 6] ]

def H (qq kl : Q) : Poly := sumP [Hrad, h0 qq kl, mul h2 Y]

def geo (qq kl : Q) : Geo :=
  let guu := addP one (scaleP (Q.ofInt 2) (H qq kl))
  let Aphi := mul F mu
  let g : Nat → Nat → Poly := fun a b =>
    match a, b with
    | 0, 1 => one | 1, 0 => one
    | 1, 1 => guu
    | 1, 5 => Aphi | 5, 1 => Aphi
    | 2, 2 => one
    | 3, 3 => r 2
    | 4, 4 => mul (r 2) (sa 2)
    | 5, 5 => mul (r 2) mu
    | _, _ => []
  -- g^tu = 1, g^tt = −(1+2H) + F² μ / r², g^tφ = −F/r², transverse block diagonal
  let gtt := addP (negP guu) (mul (mul F F) (mul mu (r (-2))))
  let gtp := negP (mul F (r (-2)))
  let gi : Nat → Nat → Poly := fun a b =>
    match a, b with
    | 0, 1 => one | 1, 0 => one
    | 0, 0 => gtt
    | 0, 5 => gtp | 5, 0 => gtp
    | 2, 2 => one
    | 3, 3 => r (-2)
    | 4, 4 => mul (r (-2)) (sa (-2))
    | 5, 5 => mul (r (-2)) (mul (sa (-2)) (s (-2)))
    | _, _ => []
  ⟨A, 6, g, gi, d⟩

def sphere : Geo := geo (Q.frac 19 96) (Q.frac (-1) 4)


/-! ### The Ricci tensor: the metric is not Ricci flat

R_uu = −Δ̂H + ¼F_{ij}F^{ij},  R_uφ = −½ J_φ,  all other components zero, with Δ̂ the
flat Laplacian of ℝ⁴, ¼F_{ij}F^{ij} = (r²F′² μ + 4F²(1 − μ))/(2r⁴) and the geometric
current J_φ = (F″ + F′/r − 4F/r²) μ = 4(3Qr⁴ − K) μ/r². -/

def lap (p : Poly) : Poly :=
  let ang := addP (d 3 (d 3 p)) (scaleP (Q.ofInt 2) (mul (mul cb (sa (-1))) (d 3 p)))
  let ang2 := mul (sa (-2)) (addP (d 4 (d 4 p)) (mul (mul c (s (-1))) (d 4 p)))
  sumP [d 2 (d 2 p), scaleP (Q.ofInt 3) (mul (r (-1)) (d 2 p)), mul (r (-2)) (addP ang ang2)]

def F4 : Poly :=
  let Fr := d 2 F
  scaleP (Q.frac 1 2) (mul (r (-4)) (addP (mul (r 2) (mul (mul Fr Fr) mu))
    (scaleP (Q.ofInt 4) (mul (mul F F) (subP one mu)))))

def Jphi : Poly := mul (subP (scaleP (Q.ofInt 12) (mul Qq (r 2))) (scaleP (Q.ofInt 4) (mul K (r (-2))))) mu

def ricciClaim (a b : Nat) : Poly :=
  match a, b with
  | 1, 1 => addP (negP (lap (H (Q.frac 19 96) (Q.frac (-1) 4)))) F4
  | 1, 5 => scaleP (Q.frac (-1) 2) Jphi
  | 5, 1 => scaleP (Q.frac (-1) 2) Jphi
  | _, _ => []

def ricciOk : Bool :=
  let Ric := sphere.ricci (sphere.riemann sphere.christoffel)
  (List.range 6).all fun a => (List.range 6).all fun b =>
    isZeroP (subP Ric[sphere.idx2 a b]! (ricciClaim a b))

/-- the current agrees with F″ + F′/r − 4F/r² -/
def currentOk : Bool :=
  isZeroP (subP Jphi (mul (sumP [d 2 (d 2 F), mul (r (-1)) (d 2 F), scaleP (Q.ofInt (-4)) (mul F (r (-2)))]) mu))

/-! ### The ends of the reduced potential

W = −2H + F² μ/r² = N² − 1 (the potential of geodesics with no angular momentum along φ).
Its only r⁻⁶ term is J² (μ − 1/6) r⁻⁶, its only r⁶ term is Q² (13μ/40 − 7/120) r⁶, and every
other term, logarithms included, has a power strictly in between.  So the centre and
infinity are walled off where μ > 1/6 (near the centre) and μ > 7/39 (far out), and open
along the rotation axis. -/

def W : Poly := addP (scaleP (Q.ofInt (-2)) (H (Q.frac 19 96) (Q.frac (-1) 4))) (mul (mul F F) (mul mu (r (-2))))
def rexp (t : Mono × Q) : Int := t.1.getD iR 0
def endsClaim : Poly :=
  addP (mul (mul J (mul J (r (-6)))) (subP mu (kq (Q.frac 1 6))))
       (mul (mul Qq (mul Qq (r 6))) (subP (scaleP (Q.frac 13 40) mu) (kq (Q.frac 7 120))))
def endsOk : Bool :=
  isZeroP (subP (W.filter fun t => rexp t == -6 || rexp t == 6) endsClaim) &&
  W.all fun t => (rexp t == -6 || rexp t == 6) || (-6 < rexp t && rexp t < 6)

end GMA.Sphere

open GMA GMA.Sphere in
theorem sphere_inverse : sphere.inverseOk = true := by native_decide

open GMA GMA.Sphere in
/-- **The three-sphere rotating along φ solves geodesic monism.**  All 36 components
    of E_ab vanish identically, for all eight constants. -/
theorem sphere_solves : sphere.fieldEq.all isZeroP = true := by native_decide

open GMA GMA.Sphere in
theorem sphere_control_Q2 : (geo (Q.frac 19 97) (Q.frac (-1) 4)).fieldEq.all isZeroP = false := by native_decide

open GMA GMA.Sphere in
theorem sphere_control_log : (geo (Q.frac 19 96) (Q.ofInt 0)).fieldEq.all isZeroP = false := by native_decide

open GMA GMA.Sphere in
/-- **Not Ricci flat.**  R_uu = −Δ̂H + ¼F² and R_uφ = −½J_φ, every other component 0. -/
theorem sphere_ricci : ricciOk = true := by native_decide

open GMA GMA.Sphere in
theorem sphere_current : currentOk = true := by native_decide

/-! ## Redshift coefficients, equator and pole

Clock ξ = ∂_t − ∂_u (any ∂_t + λ∂_u only rescales the constants, Proposition 10 of
the post):  g(ξ,ξ) = −(1 − 2H),  g(ξ,∂_φ) = −F μ,  g_φφ = r² μ,
N² = 1 − 2H + F² μ/r²,  ω = F/r²  (independent of the angles).
The metric has the isometries ∂_φ (rotations of the x₁x₂ plane) and the rotations
of the x₃x₄ plane, and depends on the angles only through μ = (x₁² + x₂²)/r².
So it is enough to work on the slice ψ = π/2, where μ = sin²θ and the sphere is
the two-sphere of the 1 + 4 metrics: the equator θ = π/2 is the rotation plane
(μ = 1), the pole θ = 0 lies on the rotation axis (μ = 0).  There Y = 1 − 2cos²θ,
H_e = H_rad + h₀ + h₂, H_p = H_rad + h₀ − h₂, N_e² = 1 − 2H_e + F²/r², and

  equator, static: C_r = −H_e′/(1 − 2H_e),        C_θ = 4h₂/(1 − 2H_e)
  equator, ZAMO:   C_r = (N_e²)′/(2N_e²),          C_θ = (4h₂ − F²/r²)/N_e²,
                   D_r = (r/N_e)(F/r²)′,           D_θ = 0
  pole, static:    C_r = −H_p′/(1 − 2H_p),        C_θ = −4h₂/(1 − 2H_p)
  pole, ZAMO:      C_r as static,                  C_θ = (F²/r² − 4h₂)/(1 − 2H_p),  D = 0 -/

namespace GMA.SphereCert
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

def Fp : Poly := sumP [m J (r (-2)), K, m B (r 2), m Qq (r 4)]
def h0 : Poly := sumP
  [ term (Q.frac 1 12) [J, J, r (-6)], term (Q.frac 1 6) [J, K, r (-4)]
  , term (Q.frac (-1) 4) [K, K, L, r (-2)], term (Q.ofInt (-1)) [J, Qq, L], term (Q.ofInt 1) [B, K, L]
  , term (Q.frac 1 4) [B, B, r 2], term (Q.frac 1 2) [K, Qq, L, r 2], term (Q.frac (-1) 8) [K, Qq, r 2]
  , term (Q.frac 1 2) [B, Qq, r 4], term (Q.frac 19 96) [Qq, Qq, r 6] ]
def h2 : Poly := sumP
  [ term (Q.frac 1 3) [J, K, L, r (-4)], term (Q.frac 1 4) [K, K, L, r (-2)], term (Q.frac 1 16) [K, K, r (-2)]
  , term (Q.frac 1 2) [B, J, r (-2)], term (Q.frac 1 2) [B, K], term (Q.frac 1 6) [K, Qq, L, r 2]
  , term (Q.frac 1 8) [B, Qq, r 4], term (Q.frac 27 160) [Qq, Qq, r 6] ]
def Hrad : Poly := sumP [m (v 7) (r (-2)), v 8, m (v 9) (r 2), m (v 10) L]
def S2 : Poly := subP (kz 1) (m c c)
def Y : Poly := subP (kz 1) (scaleP (Q.ofInt 2) (m c c))
def H : Poly := sumP [Hrad, h0, m h2 Y]
def f2r : Poly := m (m Fp Fp) (r (-2))

def data : Data :=
  { A := A, iR := 0, iC := 1, iL := 2
    X := P (subP (scaleP (Q.ofInt 2) H) (kz 1))
    Y := P (negP (m Fp S2))
    Z := P (m (r 2) S2)
    N2r := P (addP (subP (kz 1) (scaleP (Q.ofInt 2) H)) (m f2r S2))
    wr := P (m Fp (r (-2))) }

def dR' := dR data
def He := sumP [Hrad, h0, h2]
def Hp := sumP [Hrad, h0, negP h2]
def ae := subP (kz 1) (scaleP (Q.ofInt 2) He)
def ap := subP (kz 1) (scaleP (Q.ofInt 2) Hp)
def Ne2 := addP ae f2r

def equatorClaim : Claim :=
  { staticCr := F (negP (dR' He)) ae
    staticCt := F (scaleP (Q.ofInt 4) h2) ae
    zamoCr := F (dR' Ne2) (scaleP (Q.ofInt 2) Ne2)
    zamoCt := F (subP (scaleP (Q.ofInt 4) h2) f2r) Ne2 }

def equatorExtra : EquatorExtra :=
  { P2 := F (r 2) Ne2
    wr := P (dR' (m Fp (r (-2))))
    wtt := F [] (kz 1) }

def poleClaim : Claim :=
  { staticCr := F (negP (dR' Hp)) ap
    staticCt := F (scaleP (Q.ofInt (-4)) h2) ap
    zamoCr := F (negP (dR' Hp)) ap
    zamoCt := F (subP f2r (scaleP (Q.ofInt 4) h2)) ap }

/-- a deliberately wrong claim: 3h₂ in place of 4h₂ (the coefficient of the two-sphere) -/
def equatorWrong : Claim := { equatorClaim with staticCt := F (scaleP (Q.ofInt 3) h2) ae }

end GMA.SphereCert

open GMA GMA.Patch GMA.SphereCert in
theorem sphere_equator : checkClaim data .equator equatorClaim = true := by native_decide
open GMA GMA.Patch GMA.SphereCert in
theorem sphere_equator_zamo_cross : checkEquatorExtra data equatorExtra = true := by native_decide
open GMA GMA.Patch GMA.SphereCert in
theorem sphere_pole : checkClaim data .pole poleClaim = true := by native_decide
open GMA GMA.Patch GMA.SphereCert in
theorem sphere_pole_cross_term_vanishes : poleZvanishes data = true := by native_decide
open GMA GMA.Patch GMA.SphereCert in
theorem sphere_equator_control : checkClaim data .equator equatorWrong = false := by native_decide

open GMA GMA.Sphere in
/-- **The ends of W.**  W ≈ J²(μ − 1/6)/r⁶ at the centre and Q²(13μ/40 − 7/120) r⁶ far out. -/
theorem sphere_potential_ends : endsOk = true := by native_decide
