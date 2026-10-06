/-
# A small exact differential-algebra engine (no Mathlib)

* `Q`     : exact rationals, always normalised (gcd 1, positive denominator).
* `Poly`  : a finite sum of monomials Π vᵢ^{eᵢ} (eᵢ ∈ ℤ, so generators may be
            Laurent), stored as a list sorted by exponent vector with no zero
            coefficients.  This is a normal form: two `Poly` are equal as
            elements of the ring iff they are equal as lists, so "= 0" is
            "is the empty list".
* `Alg`   : the generators and an optional rewriting of monomials (used for
            c² = 1 − s² in the ring of trigonometric functions).
* `derivP`: the derivation with prescribed values on the generators, extended
            by the Leibniz rule.
* `Geo`   : coordinate tensor calculus (Christoffel, Riemann, Ricci, covariant
            derivatives) and the field equations E_ab of ∫ R_ab R^ab √(−g).
* `Frac`  : pairs (numerator, denominator) for rational-function identities,
            decided by cross multiplication.
-/

namespace GMA

/-- concatMap, written out so the file does not depend on the name the
    standard library gives it in a given Lean release -/
def cmap {α β : Type} (l : List α) (f : α → List β) : List β :=
  match l with
  | [] => []
  | x :: xs => f x ++ cmap xs f

/-! ## Exact rationals -/

structure Q where
  num : Int
  den : Nat
deriving BEq, Inhabited, Repr

namespace Q
def norm (n : Int) (d : Nat) : Q :=
  let g := Nat.gcd n.natAbs d
  if g = 0 then ⟨0, 1⟩ else ⟨n / (g : Int), d / g⟩
def ofInt (n : Int) : Q := ⟨n, 1⟩
def frac (n : Int) (d : Nat) : Q := norm n d
def add (a b : Q) : Q := norm (a.num * b.den + b.num * a.den) (a.den * b.den)
def mul (a b : Q) : Q := norm (a.num * b.num) (a.den * b.den)
def neg (a : Q) : Q := ⟨-a.num, a.den⟩
def isZero (a : Q) : Bool := a.num == 0
instance : Add Q := ⟨add⟩
instance : Mul Q := ⟨mul⟩
instance : Neg Q := ⟨neg⟩
end Q

/-! ## Polynomials: sorted lists of (exponent vector, non-zero coefficient) -/

abbrev Mono := List Int
abbrev Poly := List (Mono × Q)

def mcmp : Mono → Mono → Ordering
  | [], [] => .eq
  | [], _ => .lt
  | _, [] => .gt
  | a :: as, b :: bs => if a < b then .lt else if b < a then .gt else mcmp as bs

def addAux : Nat → Poly → Poly → Poly
  | 0, p, q => p ++ q
  | _ + 1, [], q => q
  | _ + 1, p, [] => p
  | f + 1, (m, a) :: p, (n, b) :: q =>
    match mcmp m n with
    | .lt => (m, a) :: addAux f p ((n, b) :: q)
    | .gt => (n, b) :: addAux f ((m, a) :: p) q
    | .eq => let c := a + b
             if c.isZero then addAux f p q else (m, c) :: addAux f p q

def addP (p q : Poly) : Poly := addAux (p.length + q.length + 1) p q
def scaleP (k : Q) (p : Poly) : Poly :=
  if k.isZero then [] else p.map fun (m, a) => (m, k * a)
def negP (p : Poly) : Poly := scaleP (Q.ofInt (-1)) p
def subP (p q : Poly) : Poly := addP p (negP q)
def sumP (ps : List Poly) : Poly :=
  -- balanced summation
  let rec go : Nat → List Poly → List Poly
    | 0, l => l
    | f + 1, l =>
      match l with
      | [] => []
      | [p] => [p]
      | _ =>
        let rec pairs : List Poly → List Poly
          | p :: q :: rest => addP p q :: pairs rest
          | l => l
        go f (pairs l)
  match go (ps.length + 1) ps with
  | [] => []
  | p :: _ => p

/-- An algebra: number of generators and the reduction of a monomial
    product to normal form (the relation c² = 1 − s² for 𝓡, nothing for jets). -/
structure Alg where
  nv : Nat
  reduce : Mono → Q → Poly

def monoAdd (m n : Mono) : Mono := List.zipWith (· + ·) m n
def unitMono (nv i : Nat) (e : Int) : Mono := (List.range nv).map fun j => if j = i then e else 0

def sortP (terms : List (Mono × Q)) : Poly :=
  sumP (terms.map fun t => if t.2.isZero then [] else [t])

def mulTerm (A : Alg) (m : Mono) (a : Q) (q : Poly) : Poly :=
  let raw := q.map fun (n, b) => (monoAdd m n, a * b)
  sumP (raw.map fun (mm, cc) => A.reduce mm cc)

def mulP (A : Alg) (p q : Poly) : Poly :=
  sumP (p.map fun (m, a) => mulTerm A m a q)

def constP (A : Alg) (k : Q) : Poly := if k.isZero then [] else [(unitMono A.nv 0 0, k)]
def varP (A : Alg) (i : Nat) (e : Int := 1) : Poly := [(unitMono A.nv i e, Q.ofInt 1)]

/-- The derivation with prescribed values `dv i` on the generators, extended
    by the Leibniz rule: ∂(Π vᵢ^{eᵢ}) = Σᵢ eᵢ vᵢ^{eᵢ−1} (Π_{j≠i} vⱼ^{eⱼ}) ∂vᵢ.
    The formula is also right for negative eᵢ (Laurent generators). -/
def derivP (A : Alg) (dv : Nat → Poly) (p : Poly) : Poly :=
  sumP <| cmap p fun (m, a) =>
    (List.range A.nv).filterMap fun i =>
      let e := m.getD i 0
      if e = 0 then none
      else
        let dvi := dv i
        if dvi.isEmpty then none
        else some (mulTerm A (monoAdd m (unitMono A.nv i (-1))) (a * Q.ofInt e) dvi)

def isZeroP (p : Poly) : Bool := p.isEmpty

/-! ## Coordinate tensor calculus (MTW conventions, signature − + + …)

  Γ^a_bc   = ½ g^ad (∂_b g_dc + ∂_c g_db − ∂_d g_bc)
  R^a_bcd  = ∂_c Γ^a_db − ∂_d Γ^a_cb + Γ^a_ce Γ^e_db − Γ^a_de Γ^e_cb
  R_bd     = R^a_bad
-/

structure Geo where
  A : Alg
  n : Nat
  g : Nat → Nat → Poly
  gi : Nat → Nat → Poly
  d : Nat → Poly → Poly          -- ∂_c acting on 𝓡

namespace Geo
variable (G : Geo)

def idx2 (a b : Nat) : Nat := a * G.n + b
def idx3 (a b c : Nat) : Nat := (a * G.n + b) * G.n + c
def idx4 (a b c d : Nat) : Nat := ((a * G.n + b) * G.n + c) * G.n + d
def rng : List Nat := List.range G.n
def mul (p q : Poly) : Poly := mulP G.A p q
def sum (f : Nat → Poly) : Poly := sumP (G.rng.map f)

def dg : Array Poly :=   -- dg[idx3 c a b] = ∂_c g_ab
  Id.run do
    let mut out := #[]
    for c in G.rng do
      for a in G.rng do
        for b in G.rng do
          out := out.push (G.d c (G.g a b))
    return out

def christoffel : Array Poly :=  -- Γ[idx3 a b c] = Γ^a_bc
  let dg := G.dg
  let D (c a b : Nat) := dg[G.idx3 c a b]!
  Id.run do
    let mut out := #[]
    for a in G.rng do
      for b in G.rng do
        for c in G.rng do
          out := out.push <| scaleP (Q.frac 1 2) <| G.sum fun d =>
            G.mul (G.gi a d) (subP (addP (D b d c) (D c d b)) (D d b c))
    return out

def riemann (Γ : Array Poly) : Array Poly :=  -- R[idx4 a b c d] = R^a_bcd
  let Gm (a b c : Nat) := Γ[G.idx3 a b c]!
  Id.run do
    let mut out := #[]
    for a in G.rng do
      for b in G.rng do
        for c in G.rng do
          for d in G.rng do
            let t1 := subP (G.d c (Gm a d b)) (G.d d (Gm a c b))
            let t2 := G.sum fun e => subP (G.mul (Gm a c e) (Gm e d b)) (G.mul (Gm a d e) (Gm e c b))
            out := out.push (addP t1 t2)
    return out

def ricci (Rm : Array Poly) : Array Poly :=
  Id.run do
    let mut out := #[]
    for b in G.rng do
      for d in G.rng do
        out := out.push (G.sum fun a => Rm[G.idx4 a b a d]!)
    return out

/-- E_ab of S = ∫ R_ab R^ab √(−g):
    E_ab = □R_ab + ½ g_ab □R − ∇_a∇_b R + 2 R_acbd R^cd − ½ g_ab R_cd R^cd. -/
def fieldEqWith (κ : List Q) : Array Poly :=
  let Γ := G.christoffel
  let Gm (a b c : Nat) := Γ[G.idx3 a b c]!
  let Rm := G.riemann Γ
  let Ric := G.ricci Rm
  let R2 (a b : Nat) := Ric[G.idx2 a b]!
  let Rs := G.sum fun a => G.sum fun b => G.mul (G.gi a b) (R2 a b)
  -- ∇_c R_ab
  let DR : Array Poly := Id.run do
    let mut out := #[]
    for c in G.rng do
      for a in G.rng do
        for b in G.rng do
          out := out.push <| subP (G.d c (R2 a b))
            (G.sum fun e => addP (G.mul (Gm e c a) (R2 e b)) (G.mul (Gm e c b) (R2 a e)))
    return out
  let D3 (c a b : Nat) := DR[G.idx3 c a b]!
  -- □R_ab = g^cd ∇_d ∇_c R_ab
  let boxRic (a b : Nat) : Poly := G.sum fun c => G.sum fun d =>
    let gcd := G.gi c d
    if gcd.isEmpty then [] else
    G.mul gcd <| subP (G.d d (D3 c a b)) <| G.sum fun e =>
      addP (G.mul (Gm e d c) (D3 e a b)) (addP (G.mul (Gm e d a) (D3 c e b)) (G.mul (Gm e d b) (D3 c a e)))
  let dRs : Array Poly := (G.rng.map fun e => G.d e Rs).toArray
  let hess (a b : Nat) : Poly :=
    subP (G.d b (dRs[a]!)) (G.sum fun e => G.mul (Gm e a b) (dRs[e]!))
  let boxR := G.sum fun a => G.sum fun b => G.mul (G.gi a b) (hess a b)
  let Rup (c d : Nat) : Poly := G.sum fun a => G.sum fun b => G.mul (G.gi c a) (G.mul (R2 a b) (G.gi b d))
  let RupA : Array Poly := Id.run do
    let mut out := #[]
    for c in G.rng do
      for d in G.rng do
        out := out.push (Rup c d)
    return out
  let Ric2 := G.sum fun a => G.sum fun b => G.mul (R2 a b) (RupA[G.idx2 a b]!)
  -- R_acbd = g_ae R^e_cbd
  let Rlow (a c b d : Nat) : Poly := G.sum fun e => G.mul (G.g a e) (Rm[G.idx4 e c b d]!)
  let RR (a b : Nat) : Poly := G.sum fun c => G.sum fun d =>
    let ru := RupA[G.idx2 c d]!
    if ru.isEmpty then [] else G.mul (Rlow a c b d) ru
  Id.run do
    let mut out := #[]
    for a in G.rng do
      for b in G.rng do
        let gab := G.g a b
        out := out.push <| sumP
          [ scaleP (κ.getD 0 (Q.ofInt 0)) (boxRic a b)
          , scaleP (κ.getD 1 (Q.ofInt 0)) (G.mul gab boxR)
          , scaleP (κ.getD 2 (Q.ofInt 0)) (hess a b)
          , scaleP (κ.getD 3 (Q.ofInt 0)) (RR a b)
          , scaleP (κ.getD 4 (Q.ofInt 0)) (G.mul gab Ric2) ]
    return out

/-- the coefficients of the five terms as claimed: 1, ½, −1, 2, −½ -/
def claimed : List Q := [Q.ofInt 1, Q.frac 1 2, Q.ofInt (-1), Q.ofInt 2, Q.frac (-1) 2]

def fieldEq : Array Poly := G.fieldEqWith claimed

def inverseOk : Bool :=
  G.rng.all fun a => G.rng.all fun b =>
    let p := G.sum fun c => G.mul (G.g a c) (G.gi c b)
    if a = b then p == constP G.A (Q.ofInt 1) else p.isEmpty

end Geo


/-! ## Substitution of a constant for a generator, and rational functions -/

def powQ (k : Q) : Nat → Q
  | 0 => Q.ofInt 1
  | n + 1 => k * powQ k n

/-- substitute the constant `k` for generator `i` (exponents of `i` must be ≥ 0) -/
def evalAt (i : Nat) (k : Q) (p : Poly) : Poly :=
  sortP <| p.map fun (m, a) =>
    let e := m.getD i 0
    (m.set i 0, a * powQ k e.toNat)

structure Frac where
  num : Poly
  den : Poly

namespace Frac
variable (A : Alg)
def ofP (p : Poly) : Frac := ⟨p, constP A (Q.ofInt 1)⟩
def add (x y : Frac) : Frac := ⟨addP (mulP A x.num y.den) (mulP A y.num x.den), mulP A x.den y.den⟩
def neg (x : Frac) : Frac := ⟨negP x.num, x.den⟩
def sub (x y : Frac) : Frac := add A x (neg y)
def mul (x y : Frac) : Frac := ⟨mulP A x.num y.num, mulP A x.den y.den⟩
def div (x y : Frac) : Frac := ⟨mulP A x.num y.den, mulP A x.den y.num⟩
def scale (k : Q) (x : Frac) : Frac := ⟨scaleP k x.num, x.den⟩
/-- quotient rule -/
def deriv (d : Poly → Poly) (x : Frac) : Frac :=
  ⟨subP (mulP A (d x.num) x.den) (mulP A x.num (d x.den)), mulP A x.den x.den⟩
def evalAt (i : Nat) (k : Q) (x : Frac) : Frac := ⟨GMA.evalAt i k x.num, GMA.evalAt i k x.den⟩
/-- x = y as rational functions: denominators non-zero and cross products equal -/
def eq (x y : Frac) : Bool :=
  !x.den.isEmpty && !y.den.isEmpty && isZeroP (subP (mulP A x.num y.den) (mulP A y.num x.den))
def isZero (x : Frac) : Bool := !x.den.isEmpty && x.num.isEmpty
end Frac

end GMA
