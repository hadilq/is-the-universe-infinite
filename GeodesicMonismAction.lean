import Poly

/-
# Geodesic monism: the field equations of the action, and the hydrogen solution

No Mathlib.  Everything below is exact arithmetic over ℚ, done by Lean.

The action is
    S[g] = ∫ R_ab R^ab √(-g) dⁿx ,
and the field equations claimed for it (geodesic monism II) are
    E_ab = □R_ab + ½ g_ab □R − ∇_a∇_b R + 2 R_acbd R^cd − ½ g_ab R_cd R^cd = 0 .

Part A  (`fieldEquationsFromAction`).  On the sector of four-dimensional metrics
    ds² = −e^{2f₀}(dt − p dy)² + e^{2f₁} dx² + e^{2f₂} dy² + e^{2f₃} dz²,
with five free functions f₀, f₁, f₂, f₃, p of two coordinates (x, y), Lean
computes the Euler–Lagrange expression of 𝓛 = √|g| R_ab R^ab,
    EL_φ = Σ_{|α| ≤ 2} (−D)^α ∂𝓛/∂φ_α ,
for each field φ, and proves EL_φ = −√|g| E^ab ∂g_ab/∂φ, with E_ab computed
from the formula above.  That is δS/δg^ab = √|g| E_ab on that sector, for
diagonal and off-diagonal (g_ty) components.  Changing any one of the five
coefficients 1, ½, −1, 2, −½ in E_ab breaks the identity, so the sector fixes
the whole formula.  The general derivation, valid for every metric, is
written out step by step in the post; this is its machine check.

Part B  (`hydrogen_solves`).  For the five-dimensional null Kaluza hydrogen
metric of geodesic monism II,
    ds² = 2 dt du + (1+2H) du² + 2 A_φ dφ du + dr² + r² dθ² + r² sin²θ dφ²,
    H   = c₀/r + c₁ + c₂ r + c₃ r² + J²/(12 r⁴) + J² P₂(cos θ)/(6 r⁴),
    A_φ = J sin²θ / r,
Lean computes all 25 components of E_ab and proves they vanish identically,
for all values of the constants c₀ … c₃, J.  Two negative controls show the
check is not vacuous: replacing 1/12 by 1/13, or r⁻⁴ by r⁻² in the P₂ term,
makes E_ab non-zero.

## Why a computation in a polynomial ring is a proof about functions

Part B works in the ring
    𝓡 = ℚ[c₀,c₁,c₂,c₃,J][r^{±1}, s^{±1}, c] / (c² + s² − 1),
with s = sin θ, c = cos θ, and with the two derivations
    ∂_r r = 1, ∂_r s = ∂_r c = 0;      ∂_θ s = c, ∂_θ c = −s, ∂_θ r = 0.
Both derivations preserve the ideal (∂_θ(c² + s² − 1) = −2cs + 2sc = 0), so
they are derivations of 𝓡.  Every element of 𝓡 has a unique normal form with
c-degree ≤ 1 (𝓡 is free with basis {1, c} over the Laurent ring in r, s);
`Poly` below stores exactly that normal form, so "is zero" is "is the empty
list".  Evaluation at (r, θ) with r > 0, 0 < θ < π is a ring homomorphism
from 𝓡 to smooth functions that intertwines ∂_r, ∂_θ with the partial
derivatives.  The metric, its inverse (checked: g·g⁻¹ = 1 in 𝓡), the
Christoffel symbols, Riemann, Ricci, and every covariant derivative in E_ab
are built from metric components by ring operations and these derivations
only, so E_ab = 0 in 𝓡 implies E_ab = 0 as functions.

Part A works the same way in the differential polynomial ring of the jets
∂ₓⁱ∂ᵧʲ of the five fields f₀ … f₃, p and the exponentials e^{±f_a}, with the total
derivatives D_x, D_y and the vertical partials ∂/∂f_{a,α} as derivations.

## Trust base

The polynomial identities are decided with `native_decide`, which trusts the
Lean compiler in addition to the kernel.  The library `Poly.lean` (≈ 300 lines)
is the only other thing to read: exact rationals, sorted monomial lists, and a
derivation defined by its values on generators (Leibniz rule).  Part B is
cross-checked independently with SymPy in `check/solutions.py`.
-/


namespace GMA

/-! ## Part B: the hydrogen metric of geodesic monism II

Generators of 𝓡, in this order:  r, s = sin θ, c = cos θ, J, c₀, c₁, c₂, c₃.
Coordinates: (t, u, r, θ, φ) = (0, 1, 2, 3, 4). -/

namespace Hydrogen

def NV := 8
def iR := 0
def iS := 1
def iC := 2
def iJ := 3
def iC0 := 4

/-- c² ↦ 1 − s²: a monomial with c-exponent k ≥ 2 is rewritten repeatedly. -/
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
def v (i : Nat) (e : Int := 1) : Poly := varP A i e
def mul := mulP A
def r (e : Int := 1) := v iR e
def s (e : Int := 1) := v iS e
def c := v iC

/-- ∂_r and ∂_θ on generators. -/
def dGen (coord : Nat) (i : Nat) : Poly :=
  if coord = 2 then (if i = iR then one else [])
  else if coord = 3 then
    (if i = iS then c else if i = iC then negP (s) else [])
  else []

def d (coord : Nat) (p : Poly) : Poly :=
  if coord = 2 ∨ coord = 3 then derivP A (dGen coord) p else []

/-- H with the two dipole coefficients and the power of r in the P₂ term
    left as parameters, so that the negative controls use the same code. -/
def H (k12 k6 : Q) (p2pow : Int) : Poly :=
  let P2 := addP (scaleP (Q.frac 3 2) (mul c c)) (k (Q.frac (-1) 2))
  let J2 := mul (v iJ) (v iJ)
  sumP [ mul (v iC0) (r (-1)), v (iC0 + 1), mul (v (iC0 + 2)) (r), mul (v (iC0 + 3)) (r 2)
       , scaleP k12 (mul J2 (r (-4)))
       , scaleP k6 (mul J2 (mul P2 (r (-p2pow)))) ]

def Aphi : Poly := mul (v iJ) (mul (s 2) (r (-1)))

def geo (k12 k6 : Q) (p2pow : Int) : Geo :=
  let Hh := H k12 k6 p2pow
  let guu := addP one (scaleP (Q.ofInt 2) Hh)
  let g : Nat → Nat → Poly := fun a b =>
    match a, b with
    | 0, 1 => one | 1, 0 => one
    | 1, 1 => guu
    | 1, 4 => Aphi | 4, 1 => Aphi
    | 2, 2 => one
    | 3, 3 => r 2
    | 4, 4 => mul (r 2) (s 2)
    | _, _ => []
  -- inverse: g^tu = 1, g^tt = −(1+2H) + A_φ A^φ, g^tφ = −A^φ, A^φ = A_φ/(r² s²)
  let Aup := mul (v iJ) (r (-3))                      -- J s² r⁻¹ / (r² s²)
  let gtt := addP (negP guu) (mul Aphi Aup)
  let gi : Nat → Nat → Poly := fun a b =>
    match a, b with
    | 0, 1 => one | 1, 0 => one
    | 0, 0 => gtt
    | 0, 4 => negP Aup | 4, 0 => negP Aup
    | 2, 2 => one
    | 3, 3 => r (-2)
    | 4, 4 => mul (r (-2)) (s (-2))
    | _, _ => []
  ⟨A, 5, g, gi, d⟩

def hydrogen : Geo := geo (Q.frac 1 12) (Q.frac 1 6) 4

/-! ### Not Ricci flat: the c₂ and c₃ modes fill empty space with R_uu

The only non-zero component of the Ricci tensor of the hydrogen metric is
    R_uu = −2c₂/r − 6c₃,
for every J (the J² terms of ¼F² and of −Δ̂H cancel exactly).  So the c₂ r and c₃ r²
modes, the ones that grow at infinity, are a Ricci source spread through all of space,
and the metric is not Ricci flat at infinity.  In Gauss form, with H₀ = c₀/r + c₁ + c₂ r + c₃ r²,
    (r² H₀′)′ = −r² R_uu,
so the redshift gradient H₀′(r) = (−c₀ + ∫₀ʳ (−R_uu) s² ds)/r² is the source enclosed by the
sphere of radius r over r².  A Ricci-flat exterior would leave H₀′ ∝ 1/r² and a bounded
redshift; a non-zero R_uu at every radius makes it grow with distance. -/

def c2 := v (iC0 + 2)
def c3 := v (iC0 + 3)
def Ruu : Poly := addP (scaleP (Q.ofInt (-2)) (mul c2 (r (-1)))) (scaleP (Q.ofInt (-6)) c3)
def H0 : Poly := sumP [mul (v iC0) (r (-1)), v (iC0 + 1), mul c2 (r), mul c3 (r 2)]

def ricciOk : Bool :=
  let G := hydrogen
  let Ric := G.ricci (G.riemann G.christoffel)
  (List.range 5).all fun a => (List.range 5).all fun b =>
    isZeroP (subP Ric[G.idx2 a b]! (if a = 1 ∧ b = 1 then Ruu else []))

def fluxOk : Bool :=
  isZeroP (addP (d 2 (mul (r 2) (d 2 H0))) (mul (r 2) Ruu))

end Hydrogen

/-- The inverse metric used below really is the inverse: g·g⁻¹ = 1 in 𝓡. -/
theorem hydrogen_inverse : Hydrogen.hydrogen.inverseOk = true := by native_decide

/-- **Theorem (hydrogen solves the field equations).**  All 25 components of
    E_ab vanish identically in 𝓡, for every c₀, c₁, c₂, c₃, J. -/
theorem hydrogen_solves : Hydrogen.hydrogen.fieldEq.all isZeroP = true := by native_decide

/-- Negative control: J²/(13 r⁴) in place of J²/(12 r⁴) is not a solution. -/
theorem control_coefficient : (Hydrogen.geo (Q.frac 1 13) (Q.frac 1 6) 4).fieldEq.all isZeroP = false := by
  native_decide

/-- Negative control: J² P₂/(6 r²) in place of J² P₂/(6 r⁴) is not a solution. -/
theorem control_power : (Hydrogen.geo (Q.frac 1 12) (Q.frac 1 6) 2).fieldEq.all isZeroP = false := by
  native_decide

/-- **The hydrogen is not Ricci flat.**  Its only non-zero Ricci component is
    R_uu = −2c₂/r − 6c₃, for every c₀, c₁, J. -/
theorem hydrogen_ricci : Hydrogen.ricciOk = true := by native_decide

/-- **Gauss form.**  (r² H₀′)′ = −r² R_uu: the redshift gradient is the enclosed R_uu over r². -/
theorem hydrogen_flux : Hydrogen.fluxOk = true := by native_decide

/-! ## Part A: the field equations are the Euler–Lagrange equations of the action

Sector: four-dimensional metrics in coordinates (t, x, y, z),
    ds² = −e^{2f₀}(dt − p dy)² + e^{2f₁} dx² + e^{2f₂} dy² + e^{2f₃} dz²,
with five free functions f₀, f₁, f₂, f₃, p of (x, y).  The shift p makes g_ty
non-zero, so off-diagonal components of E_ab are tested too.  The inverse is
polynomial: g^tt = −e^{−2f₀} + p² e^{−2f₂}, g^ty = p e^{−2f₂}, g^yy = e^{−2f₂},
and √|g| = e^{f₀+f₁+f₂+f₃}.

Generators: E₀ … E₃ with E_a = e^{f_a} (Laurent), then the jets
φ_{(i,j)} = ∂ₓⁱ ∂ᵧʲ φ, i + j ≤ K, of the five fields φ = f₀ … f₃, p, then one
sentinel generator that every jet of order > K is sent to (an order overflow
can only make a test fail, never pass). -/

namespace Action

def D := 4
def NF := 5          -- fields f₀, f₁, f₂, f₃, p
def K := 5
def jets : List (Nat × Nat) :=
  cmap (List.range (K + 1)) fun o => (List.range (o + 1)).map fun i => (i, o - i)
def NJ := jets.length
def NV := D + NF * NJ + 1
def sentinel := NV - 1

def findIndex (l : List (Nat × Nat)) (x : Nat × Nat) : Option Nat :=
  let rec go : List (Nat × Nat) → Nat → Option Nat
    | [], _ => none
    | y :: ys, k => if y == x then some k else go ys (k + 1)
  go l 0

def jetIdx (φ i j : Nat) : Nat :=
  match findIndex jets (i, j) with
  | some k => D + φ * NJ + k
  | none => sentinel

def A : Alg := ⟨NV, fun m a => [(m, a)]⟩
def one : Poly := constP A (Q.ofInt 1)
def v (i : Nat) (e : Int := 1) : Poly := varP A i e
def mul := mulP A

/-- decode a generator index into (field, i, j) if it is a jet -/
def decode (g : Nat) : Option (Nat × Nat × Nat) :=
  if g < D ∨ g ≥ sentinel then none
  else
    let φ := (g - D) / NJ
    let k := (g - D) % NJ
    match jets[k]? with
    | some (i, j) => some (φ, i, j)
    | none => none

/-- total derivative D_x (dir = 0) or D_y (dir = 1) on generators:
    D E_a = (D f_a) E_a,   D φ_{(i,j)} = φ_{(i+1,j)} or φ_{(i,j+1)} -/
def totGen (dir : Nat) (g : Nat) : Poly :=
  if g < D then mul (v (jetIdx g (if dir = 0 then 1 else 0) (if dir = 0 then 0 else 1))) (v g)
  else if g = sentinel then v sentinel
  else match decode g with
    | some (φ, i, j) => if dir = 0 then v (jetIdx φ (i + 1) j) else v (jetIdx φ i (j + 1))
    | none => []

def tot (dir : Nat) (p : Poly) : Poly := derivP A (totGen dir) p

def d (coord : Nat) (p : Poly) : Poly :=
  if coord = 1 then tot 0 p else if coord = 2 then tot 1 p else []

def P : Poly := v (jetIdx 4 0 0)
def E (a : Nat) (e : Int) : Poly := v a e

def geo : Geo :=
  let g : Nat → Nat → Poly := fun a b =>
    match a, b with
    | 0, 0 => negP (E 0 2)
    | 0, 2 => mul (E 0 2) P | 2, 0 => mul (E 0 2) P
    | 2, 2 => subP (E 2 2) (mul (E 0 2) (mul P P))
    | 1, 1 => E 1 2
    | 3, 3 => E 3 2
    | _, _ => []
  let gi : Nat → Nat → Poly := fun a b =>
    match a, b with
    | 0, 0 => addP (negP (E 0 (-2))) (mul (mul P P) (E 2 (-2)))
    | 0, 2 => mul P (E 2 (-2)) | 2, 0 => mul P (E 2 (-2))
    | 2, 2 => E 2 (-2)
    | 1, 1 => E 1 (-2)
    | 3, 3 => E 3 (-2)
    | _, _ => []
  ⟨A, D, g, gi, d⟩

/-- √|g| = E₀E₁E₂E₃ -/
def sqrtg : Poly := (List.range D).foldl (fun acc a => mul acc (v a)) one

/-- the Lagrangian density 𝓛 = √|g| R_ab R^ab -/
def lagrangian : Poly :=
  let G := geo
  let Ric := G.ricci (G.riemann G.christoffel)
  let R2 (a b : Nat) := Ric[G.idx2 a b]!
  let Rup (c d : Nat) : Poly := G.sum fun a => G.sum fun b => G.mul (G.gi c a) (G.mul (R2 a b) (G.gi b d))
  let Ric2 := G.sum fun a => G.sum fun b => G.mul (R2 a b) (Rup a b)
  mul sqrtg Ric2

/-- vertical partial ∂/∂φ_{(i,j)}; for a metric exponent f_a at order (0,0)
    it acts on E_a = e^{f_a} by ∂E_a/∂f_a = E_a -/
def vert (φ i j : Nat) (p : Poly) : Poly :=
  derivP A (fun g => if φ < D ∧ i = 0 ∧ j = 0 then (if g = φ then v φ else [])
                     else if g = jetIdx φ i j then one else []) p

def iterate (f : Poly → Poly) : Nat → Poly → Poly
  | 0, p => p
  | n + 1, p => iterate f n (f p)

/-- EL_φ = Σ_{i+j ≤ 2} (−1)^{i+j} D_xⁱ D_yʲ ∂𝓛/∂φ_{(i,j)} -/
def eulerLagrange (L : Poly) (φ : Nat) : Poly :=
  sumP <| [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2)].map fun (i, j) =>
    let t := iterate (tot 1) j (iterate (tot 0) i (vert φ i j L))
    if (i + j) % 2 = 0 then t else negP t

/-- the prediction of δS = ∫ √|g| E_ab δg^ab = −∫ √|g| E^ab δg_ab :
    EL_φ = −√|g| E^ab ∂g_ab/∂φ -/
def predicted (Ecov : Array Poly) (φ : Nat) : Poly :=
  let G := geo
  let Eup (a b : Nat) : Poly := G.sum fun c => G.sum fun d =>
    G.mul (G.gi a c) (G.mul (Ecov[G.idx2 c d]!) (G.gi d b))
  negP <| mul sqrtg <| G.sum fun a => G.sum fun b =>
    let dg := vert φ 0 0 (G.g a b)
    if dg.isEmpty then [] else G.mul (Eup a b) dg

def checkWith (κ : List Q) : Bool :=
  let L := lagrangian
  let Ecov := geo.fieldEqWith κ
  (List.range NF).all fun φ => isZeroP (subP (eulerLagrange L φ) (predicted Ecov φ))

def check : Bool := checkWith Geo.claimed

/-- change the i-th coefficient by +1 -/
def bumped (i : Nat) : List Q :=
  (List.range 5).map fun j => if i = j then Geo.claimed.getD j (Q.ofInt 0) + Q.ofInt 1
                               else Geo.claimed.getD j (Q.ofInt 0)

/-- every one of the five coefficients is pinned by the action -/
def pinned : Bool := (List.range 5).all fun i => !checkWith (bumped i)

/-- non-vacuity: E_ab has non-zero diagonal and off-diagonal (t,y) components -/
def nontrivial : Bool :=
  let E := geo.fieldEq
  !(E[geo.idx2 0 0]!).isEmpty && !(E[geo.idx2 0 2]!).isEmpty && !(E[geo.idx2 1 2]!).isEmpty

end Action

theorem action_sector_inverse : Action.geo.inverseOk = true := by native_decide

/-- **Theorem (field equations of the action).**  On the sector above, for
    each of the five free functions φ, the Euler–Lagrange expression of
    √|g| R_ab R^ab equals −√|g| E^ab ∂g_ab/∂φ, with E_ab given by
    □R_ab + ½ g_ab □R − ∇_a∇_b R + 2 R_acbd R^cd − ½ g_ab R_cd R^cd.
    That is δS/δg^ab = √|g| E_ab, signs and factors included. -/
theorem fieldEquationsFromAction : Action.check = true := by native_decide

theorem fieldEquations_nontrivial : Action.nontrivial = true := by native_decide

/-- Changing any single coefficient of E_ab (1, ½, −1, 2, −½) by one breaks
    the identity: the sector fixes all five numbers. -/
theorem fieldEquations_coefficients_pinned : Action.pinned = true := by native_decide

end GMA
