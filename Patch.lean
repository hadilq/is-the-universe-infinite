import Poly

/-
# Redshift coefficients of a short climbing ray, equator and pole patches

A stationary axisymmetric metric enters only through three functions of
(r, θ), written in c = cos θ (sin²θ = 1 − c²):
    X = g(ξ, ξ),   Y = g(ξ, ∂_φ),   Z = g_φφ,
with ξ the clock (timelike Killing vector).  From them
    α² = −X                     (static observers  U = ξ/α)
    N² = −X + Y²/Z,  w = −Y/Z   (ZAMOs             U = (ξ + w ∂_φ)/N).

Proposition 5 of the post gives ln(1+z) exactly; its Taylor expansion (Proposition 6) in the regime
r ≫ Δx ≫ δr (Δθ, ϑ ~ ε, δr/r ~ ε²) is
    z = C_r δr + C_t (ϑ Δθ − Δθ²/2) + n_φ [D_r δr + D_t (ϑ Δθ − Δθ²/2)] + O(ε³)
with
    static: C_r = ½ ∂_r ln α²,  C_t = ½ ∂²_θ ln α²,  D = 0
    ZAMO:   C_r = ½ ∂_r ln N²,  C_t = ½ ∂²_θ ln N²,
            D_r = P ∂_r w,  D_t = P ∂²_θ w,  P = √(Z/N²)   (equator; D = O(ε) at the pole)
all evaluated at the patch centre, provided the first θ-derivatives vanish
there (checked below).  On c: ∂_θ = −s ∂_c and ∂²_θ = s² ∂²_c − c ∂_c, so
    equator (c = 0, s = 1):  ∂_θ F = −F_c,   ∂²_θ F = F_cc
    pole    (c = 1, s = 0):  ∂_θ F = 0,      ∂²_θ F = −F_c.
In the observer's local coordinates (Δx, γ) — √g_θθ Δθ = Δx sin γ,
√g_φφ Δφ = Δx cos γ, ℓ = √g_θθ — the bracket ϑΔθ − Δθ²/2 becomes
ϑΔx sinγ/ℓ − Δx² sin²γ/(2ℓ²) on the equator and ϑΔx sinγ/ℓ − Δx²/(2ℓ²) at
the pole, and n_φ = cos γ (check/expansion_check.py verifies this on
integrated rays).  The coefficients C, D are the same in both forms.
This file checks claimed closed forms of every coefficient as identities of
rational functions (cross multiplication, non-zero denominators).
-/

namespace GMA.Patch

structure Data where
  A : Alg
  iR : Nat
  iC : Nat
  X : Frac
  Y : Frac
  Z : Frac
  /-- N² and w with the common factor sin²θ cancelled by hand, so that they
      can be evaluated on the axis; `reducedOk` checks them against the
      definitions as rational functions. -/
  N2r : Frac
  wr : Frac
  /-- index of a generator standing for ln r (∂_r ln r = r⁻¹); none if absent -/
  iL : Nat := 1000

variable (D : Data)

def one : Poly := constP D.A (Q.ofInt 1)
def dR : Poly → Poly :=
  derivP D.A (fun i => if i = D.iR then one D else if i = D.iL then varP D.A D.iR (-1) else [])
def dC : Poly → Poly := derivP D.A (fun i => if i = D.iC then one D else [])
def fdR (x : Frac) : Frac := Frac.deriv D.A (dR D) x
def fdC (x : Frac) : Frac := Frac.deriv D.A (dC D) x

def alpha2 : Frac := Frac.neg D.X
def N2def : Frac := Frac.add D.A (Frac.neg D.X) (Frac.div D.A (Frac.mul D.A D.Y D.Y) D.Z)
def wdef : Frac := Frac.neg (Frac.div D.A D.Y D.Z)
def reducedOk : Bool := Frac.eq D.A (N2def D) D.N2r && Frac.eq D.A (wdef D) D.wr
def N2 : Frac := D.N2r
def w : Frac := D.wr

/-- ∂F/∂r, ∂F/∂c, ∂²F/∂c² of F = ln x -/
def logR (x : Frac) : Frac := Frac.div D.A (fdR D x) x
def logC (x : Frac) : Frac := Frac.div D.A (fdC D x) x
def logCC (x : Frac) : Frac := fdC D (logC D x)

def at0 (x : Frac) : Frac := Frac.evalAt D.iC (Q.ofInt 0) x
def at1 (x : Frac) : Frac := Frac.evalAt D.iC (Q.ofInt 1) x

inductive Patch | equator | pole

/-- first θ-derivative vanishes at the patch centre (pole: s = 0, but the
    c-derivative must be finite, i.e. its denominator non-zero at c = 1) -/
def thetaFlat (p : Patch) (fc : Frac) : Bool :=
  match p with
  | .equator => (at0 D fc).isZero
  | .pole => !(at1 D fc).den.isEmpty

def logTheta2 (p : Patch) (x : Frac) : Frac :=
  match p with
  | .equator => at0 D (logCC D x)
  | .pole => Frac.neg (at1 D (logC D x))

def at_ (p : Patch) (x : Frac) : Frac := match p with | .equator => at0 D x | .pole => at1 D x

def half (x : Frac) : Frac := Frac.scale (Q.frac 1 2) x

/-- C_r and C_t of ln x at the patch -/
def Cr (p : Patch) (x : Frac) : Frac := half (at_ D p (logR D x))
def Ct (p : Patch) (x : Frac) : Frac := half (logTheta2 D p x)

structure Claim where
  staticCr : Frac
  staticCt : Frac
  zamoCr : Frac
  zamoCt : Frac

structure EquatorExtra where
  P2 : Frac      -- P² = Z/N² at c = 0
  wr : Frac      -- ∂_r w at c = 0
  wtt : Frac     -- ∂²_θ w at c = 0

def checkClaim (p : Patch) (cl : Claim) : Bool :=
  let eq := Frac.eq D.A
  reducedOk D && thetaFlat D p (logC D (alpha2 D)) && thetaFlat D p (logC D (N2 D)) &&
  eq (Cr D p (alpha2 D)) cl.staticCr && eq (Ct D p (alpha2 D)) cl.staticCt &&
  eq (Cr D p (N2 D)) cl.zamoCr && eq (Ct D p (N2 D)) cl.zamoCt

def checkEquatorExtra (ex : EquatorExtra) : Bool :=
  let eq := Frac.eq D.A
  (at0 D (fdC D (w D))).isZero &&                               -- ∂_θ w = 0 on the equator
  eq (at0 D (Frac.div D.A D.Z (N2 D))) ex.P2 &&
  eq (at0 D (fdR D (w D))) ex.wr &&
  eq (at0 D (fdC D (fdC D (w D)))) ex.wtt

/-- At the pole √Z → 0, so the ZAMO cross term is O(ε³): Z vanishes at c = 1. -/
def poleZvanishes : Bool := (at1 D D.Z).isZero

end GMA.Patch
