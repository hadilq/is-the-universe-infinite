import Lake
open Lake DSL

package redshift where

@[default_target]
lean_lib Proofs where
  roots := #[`Poly, `Patch, `GeodesicMonismAction, `GeodesicMonismRotating, `GeodesicMonism6D, `GeodesicMonismSphere, `Kerr, `MyersPerry, `GeodesicMonism]
