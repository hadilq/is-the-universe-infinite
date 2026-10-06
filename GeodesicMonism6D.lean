import Patch

/-
# Geodesic monism in 1 + 1 + 4: the three-sphere rotating equally in both planes

(The closing subsection of the 1 + 1 + 4 section of the post, where the angle θ below is
called χ and ψ is called φ′; the main solution there, rotating along φ only, is
GeodesicMonismSphere.lean.)

In the post the metric is written in Hopf coordinates (θ ∈ [0, π/2]):
    ds² = 2 dt du + (1 + 2H(r)) du² + 2 F(r) (sin²θ dφ + cos²θ dψ) du
          + dr² + r² (dθ² + sin²θ dφ² + cos²θ dψ²),
    F = J/r² + K + B r² + Q r⁴,
    H = c₀/r² + c₁ + c₂ r² + c₃ ln r
        + J²/(6r⁶) + JK/(3r⁴) − K² ln r/(2r²) + KQ r² ln r + BQ r⁴ + 19Q²r⁶/48.
The rotation is in the φ plane and the ψ plane at the same rate: its direction is
χ = ∂_φ + ∂_ψ, the Hopf fibre, which is ∂_φ near θ = π/2 and ∂_ψ near θ = 0.

Here the same metric is written in Euler angles θ_E = 2θ, φ_E = ψ − φ, ψ_E = φ + ψ
(check/sixd.py verifies the map), where
    dr² + r²(dθ² + sin²θ dφ² + cos²θ dψ²) = dr² + (r²/4)(σ₁² + σ₂² + σ₃²),
    σ₁² + σ₂² = dθ_E² + sin²θ_E dφ_E²,   σ₃ = dψ_E + cos θ_E dφ_E = 2(sin²θ dφ + cos²θ dψ),
so the rotation term is 2F(sin²θ dφ + cos²θ dψ) du = 2 f σ₃ du with f = F/2.  In
Euler angles every component lies in
𝓡 = ℚ[consts][r^{±1}, ln r, s^{±1}, c]/(c² + s² − 1),  s = sin θ_E, c = cos θ_E.
Lean computes all 36 components of E_ab and proves they vanish identically, for
all eight constants.  Controls: 19/48 → 19/49 and dropping K² ln r/(2r²) both break it.
-/

namespace GMA.SixD

/-- generators: r, s, c, L = ln r, J, K, B, Q, c₀, c₁, c₂, c₃ -/
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
def kq (q : Q) : Poly := constP A q
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

/-- coordinates (t, u, r, θ, φ, ψ) = (0, 1, 2, 3, 4, 5); ∂_r and ∂_θ on generators -/
def dGen (coord : Nat) (i : Nat) : Poly :=
  if coord = 2 then (if i = iR then one else if i = iL then r (-1) else [])
  else if coord = 3 then (if i = iS then c else if i = iC then negP (s) else [])
  else []

def d (coord : Nat) (p : Poly) : Poly :=
  if coord = 2 ∨ coord = 3 then derivP A (dGen coord) p else []

def term (q : Q) (fs : List Poly) : Poly := scaleP q (fs.foldl mul one)

/-- F of the post, and f = F/2 the coefficient of σ₃ -/
def Fp : Poly := sumP [mul J (r (-2)), K, mul B (r 2), mul Qq (r 4)]
def f : Poly := scaleP (Q.frac 1 2) Fp

def Hrad : Poly := sumP [mul (v iC0) (r (-2)), v (iC0 + 1), mul (v (iC0 + 2)) (r 2), mul (v (iC0 + 3)) L]

/-- the rotation part of H, with two coefficients left free for the controls -/
def Hrot (qq kl : Q) : Poly := sumP
  [ term (Q.frac 1 6) [J, J, r (-6)], term (Q.frac 1 3) [J, K, r (-4)]
  , term kl [K, K, L, r (-2)], term (Q.ofInt 1) [K, Qq, L, r 2]
  , term (Q.ofInt 1) [B, Qq, r 4], term qq [Qq, Qq, r 6] ]

def geo (qq kl : Q) : Geo :=
  let H := addP Hrad (Hrot qq kl)
  let guu := addP one (scaleP (Q.ofInt 2) H)
  let q4 := Q.frac 1 4
  let g : Nat → Nat → Poly := fun a b =>
    match a, b with
    | 0, 1 => one | 1, 0 => one
    | 1, 1 => guu
    | 1, 4 => mul f c | 4, 1 => mul f c
    | 1, 5 => f | 5, 1 => f
    | 2, 2 => one
    | 3, 3 => scaleP q4 (r 2)
    | 4, 4 => scaleP q4 (r 2)
    | 5, 5 => scaleP q4 (r 2)
    | 4, 5 => scaleP q4 (mul (r 2) c) | 5, 4 => scaleP q4 (mul (r 2) c)
    | _, _ => []
  -- inverse: g^tu = 1, g^tt = −(1+2H) + 4f²/r², g^tψ = −4f/r², g^tφ = 0,
  -- transverse block: h^rr = 1, h^θθ = 4/r², (φ,ψ): (4/(r² sin²θ)) [[1, −cos θ], [−cos θ, 1]]
  let gtt := addP (negP guu) (scaleP (Q.ofInt 4) (mul (mul f f) (r (-2))))
  let gtpsi := scaleP (Q.ofInt (-4)) (mul f (r (-2)))
  let blk := scaleP (Q.ofInt 4) (mul (r (-2)) (s (-2)))
  let gi : Nat → Nat → Poly := fun a b =>
    match a, b with
    | 0, 1 => one | 1, 0 => one
    | 0, 0 => gtt
    | 0, 5 => gtpsi | 5, 0 => gtpsi
    | 2, 2 => one
    | 3, 3 => scaleP (Q.ofInt 4) (r (-2))
    | 4, 4 => blk
    | 5, 5 => blk
    | 4, 5 => negP (mul blk c) | 5, 4 => negP (mul blk c)
    | _, _ => []
  ⟨A, 6, g, gi, d⟩

def sixD : Geo := geo (Q.frac 19 48) (Q.frac (-1) 2)

end GMA.SixD

open GMA GMA.SixD in
theorem sixD_inverse : sixD.inverseOk = true := by native_decide

open GMA GMA.SixD in
/-- **Theorem (the rotating 1 + 1 + 4 solution).**  All 36 components of E_ab
    vanish identically, for all c₀, c₁, c₂, c₃, J, K, B, Q. -/
theorem sixD_solves : sixD.fieldEq.all isZeroP = true := by native_decide

open GMA GMA.SixD in
theorem sixD_control_Q2 : (geo (Q.frac 19 49) (Q.frac (-1) 2)).fieldEq.all isZeroP = false := by native_decide

open GMA GMA.SixD in
theorem sixD_control_log : (geo (Q.frac 19 48) (Q.ofInt 0)).fieldEq.all isZeroP = false := by native_decide

/-! ## Redshift coefficients and the radial potential

Clock ξ = ∂_t − ∂_u (any other clock ∂_t + λ∂_u only rescales the constants, see
the post), rotation along χ = ∂_φ + ∂_ψ = 2∂_{ψ_E} (the Hopf fibre):
    g(ξ,ξ) = −(1 − 2H),   g(ξ,χ) = −F,   g_χχ = r²,
    N² = 1 − 2H + F²/r²,   ω = F/r²   (the same at every point of the three-sphere).
Here the Killing vector ∂_{ψ_E} = χ/2 is used, with g(ξ,∂_{ψ_E}) = −f = −F/2 and
g_{ψ_Eψ_E} = r²/4; N² and the redshift coefficients do not depend on that scale.
Nothing depends on θ, so every patch of the sphere has
    static: C_r = −H′/(1 − 2H),  C_θ = 0;   ZAMO: C_r = (N²)′/(2N²),  C_θ = 0,
    D_r = (r/N) (F/r²)′,  D_θ = 0.
The radial potential of the zero-angular-momentum geodesics is W = N² − 1 = −2H + F²/r²
(W, not U: U is the observers' four-velocity in the post). -/

namespace GMA.SixDCert
open GMA.Patch

def A : Alg := ⟨11, fun m a => [(m, a)]⟩
def v (i : Nat) (e : Int := 1) : Poly := varP A i e
def kz (n : Int) : Poly := constP A (Q.ofInt n)
def m := mulP A
def F (n d : Poly) : Frac := ⟨n, d⟩
def P (n : Poly) : Frac := Frac.ofP A n
def r (e : Int := 1) := v 0 e
def L := v 2
def J := v 3
def K := v 4
def B := v 5
def Qq := v 6
def term (q : Q) (fs : List Poly) : Poly := scaleP q (fs.foldl m (kz 1))
def Fp : Poly := sumP [m J (r (-2)), K, m B (r 2), m Qq (r 4)]
def f : Poly := scaleP (Q.frac 1 2) Fp
def H : Poly := sumP
  [ m (v 7) (r (-2)), v 8, m (v 9) (r 2), m (v 10) L
  , term (Q.frac 1 6) [J, J, r (-6)], term (Q.frac 1 3) [J, K, r (-4)]
  , term (Q.frac (-1) 2) [K, K, L, r (-2)], term (Q.ofInt 1) [K, Qq, L, r 2]
  , term (Q.ofInt 1) [B, Qq, r 4], term (Q.frac 19 48) [Qq, Qq, r 6] ]
def a0 := subP (kz 1) (scaleP (Q.ofInt 2) H)
def f2 := scaleP (Q.ofInt 4) (m (m f f) (r (-2)))
def N2 := addP a0 f2

def data : Data :=
  { A := A, iR := 0, iC := 1, iL := 2
    X := P (negP a0)
    Y := P (negP f)
    Z := P (scaleP (Q.frac 1 4) (r 2))
    N2r := P N2
    wr := P (scaleP (Q.ofInt 4) (m f (r (-2)))) }   -- ω along ∂_{ψ_E}: 4f/r² = 2F/r² (twice the ω along χ)

def dR' := dR data
def claim : Claim :=
  { staticCr := F (negP (dR' H)) a0
    staticCt := F [] (kz 1)
    zamoCr := F (dR' N2) (scaleP (Q.ofInt 2) N2)
    zamoCt := F [] (kz 1) }
def extra : EquatorExtra :=
  { P2 := F (scaleP (Q.frac 1 4) (r 2)) N2
    wr := P (dR' (scaleP (Q.ofInt 4) (m f (r (-2)))))
    wtt := F [] (kz 1) }

/-- W = −2H + 4f²/r² = −2H + F²/r² -/
def W : Poly := addP (scaleP (Q.ofInt (-2)) H) f2
def rexp (t : Mono × Q) : Int := t.1.getD 0 0
def logFree (t : Mono × Q) : Bool := t.1.getD 2 0 == 0
def coeffOf (p : Poly) (mono : Mono) : Q :=
  match p.find? (fun t => t.1 == mono) with
  | some t => t.2
  | none => Q.ofInt 0
def mono (re : Int) (extra : List (Nat × Int)) : Mono :=
  (List.range 11).map fun i => if i = 0 then re else
    match extra.find? (·.1 == i) with
    | some (_, e) => e
    | none => 0

/-- the ends of W: the only r⁻⁶ term is (2/3) J² r⁻⁶, the only r⁶ term is (5/24) Q² r⁶,
    and every other term, logarithms included, has −6 < power < 6 -/
def endsOk : Bool :=
  (W.filter fun t => rexp t == -6) == [(mono (-6) [(3, 2)], Q.frac 2 3)] &&
  (W.filter fun t => rexp t == 6) == [(mono 6 [(6, 2)], Q.frac 5 24)] &&
  W.all fun t => (rexp t == -6 || rexp t == 6) || (-6 < rexp t && rexp t < 6)

end GMA.SixDCert

open GMA GMA.Patch GMA.SixDCert in
/-- the coefficients hold at both ends of the θ range (θ_E = π/2, the Clifford torus θ = π/4, and θ_E → 0, the circle θ = 0) -/
theorem sixD_patch_equator : checkClaim data .equator claim = true := by native_decide
open GMA GMA.Patch GMA.SixDCert in
theorem sixD_patch_pole : checkClaim data .pole claim = true := by native_decide
open GMA GMA.Patch GMA.SixDCert in
theorem sixD_patch_rotation : checkEquatorExtra data extra = true := by native_decide
open GMA GMA.SixDCert in
/-- **The ends of the radial potential.**  W → +∞ at r → 0 (as (2/3)J²/r⁶) and
    at r → ∞ (as (5/24)Q²r⁶) whenever J ≠ 0 and Q ≠ 0. -/
theorem sixD_potential_ends : endsOk = true := by native_decide
