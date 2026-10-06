# Proofs and checks for "Is the universe infinite?"

A light ray climbs a short way (`r >> Δx >> δr`) across a patch of a rotating
metric, near the equator and near the pole.  Light is a bundle of null
geodesics; the period is the length of the timelike separation between two
neighbouring rays, Lie transported along the bundle (Theorems 1–3 and
Proposition 5 of the post).  The redshift is written as z(Δx, γ) in the
observer's local coordinates (Proposition 6), with its leading order in Corollary 7.

## Lean (no Mathlib), `lake build`

| file | what it proves |
|---|---|
| `Poly.lean` | the engine: exact rationals, Laurent polynomials in normal form, derivations, coordinate tensor calculus, rational-function identities |
| `Patch.lean` | the coefficients of Proposition 6 as a checker |
| `Kerr.lean` | every redshift coefficient of Kerr, equator and pole, static and ZAMO observers, plus negative controls |
| `MyersPerry.lean` | the same for five-dimensional singly spinning Myers–Perry |
| `GeodesicMonism.lean` | the same for the geodesic-monism hydrogen metric |
| `GeodesicMonismRotating.lean` | the general rotating dipole in 1+4 (J, K, B, Q modes, ln r back-reaction) solves E_ab = 0; its redshift coefficients |
| `GeodesicMonismSphere.lean` | the three-sphere rotating along φ in 1+1+4, (t, u, r, ψ, θ, φ): solves E_ab = 0 (36 components); its Ricci tensor (not Ricci flat); its redshift coefficients on the equator (rotation plane) and the pole (rotation axis); the ends of the reduced potential W |
| `GeodesicMonism6D.lean` | the three-sphere rotating equally in both planes (rotation ∂_φ + ∂_φ′): solves E_ab = 0; no angular dependence, walls in every direction. Written in Euler angles θ_E = 2χ, φ_E = φ′ − φ, ψ_E = φ + φ′, where F = 2f |
| `GeodesicMonismAction.lean` | the field equations of `∫ R_ab R^ab √(−g)` are the Euler–Lagrange equations (five-function sector, all five coefficients pinned); the hydrogen metric solves them (all 25 components); negative controls; the hydrogen is not Ricci flat: R_uu = −2c₂/r − 6c₃, all other components 0, and (r²H₀′)′ = −r²R_uu |

Identities are decided with `native_decide` (trusts the Lean compiler).
Tested with Lean 4.10 (the version in the pinned nixpkgs) and 4.12.

## Python, `sh check/run_all.sh [--slow] [PANTHEON_DIR]`

| file | what it checks |
|---|---|
| `check/solutions.py` | Kerr and Myers–Perry are Ricci flat; Myers–Perry without r² in `g_rr` is not; with `--slow`, `E_ab = 0` for geodesic monism in SymPy (~20 min) |
| `check/redshift.py` | derives all coefficients and compares with the closed forms of the post |
| `check/lie_transport.py` | shoots two neighbouring null geodesics, reads the period on the observer's worldline, compares with Proposition 5 |
| `check/expansion_check.py` | z(Δx, γ) of Proposition 6 against the exact Proposition 5 on integrated rays, with ε³ scaling |
| `check/rotating.py` | derivation of the rotating 1+4 dipole family from the reduced equations |
| `check/sphere.py` | the three-sphere rotating along φ: derivation, R_uu, the patch coefficients, and z(Δx, γ, β) against integrated rays |
| `check/sphere_orbits.py` | its stable orbit, the barrier, the band r_min–r_max of every geodesic from the constants (checked on integrated geodesics), and the leak through the rotation axis |
| `check/sphere_fit.py` | the three-sphere against all the data: Pantheon+ and DES-Dovekie fits, the dipole fixed toward Shapley with the rotation axis scanned, CMB dipole and quadrupole limits, curvature bound on r, bulk-flow comparison (needs pandas and astropy) |
| `check/hubble.py` | the Hubble law of the non-Ricci-flat modes for static and drifting observers (slope, q0), and the shear-free focusing that keeps images sharp |
| `check/sixd.py`, `check/trapping.py` | the three-sphere rotating equally in both planes: derivation and walls |
| `check/pantheon_fit.py` | dipole, quadrupole, Corollary 7 (leading order) and full z(Δx, γ) fits to Pantheon+ with the full covariance |

`PANTHEON_DIR` must contain `Pantheon+SH0ES.dat` and `Pantheon+SH0ES_STAT+SYS.cov`
from https://github.com/PantheonPlusSH0ES/DataRelease, and `des/repo/4_DISTANCES_COVMAT/`
from https://github.com/des-science/DES-SN5YR (`git clone --depth 1 --filter=blob:none --sparse` then
`git sparse-checkout set 4_DISTANCES_COVMAT`).

## Nix

`nix build` runs `lake build` and the fast checks; `nix build .#full` adds `--slow`.
