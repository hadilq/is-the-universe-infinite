{
  description = "Redshift of a climbing light ray in Kerr, Myers-Perry and geodesic monism: Lean certificates and checks";
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.11";
  outputs = { self, nixpkgs }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" "x86_64-darwin" "aarch64-darwin" ];
      build = system: slow:
        let pkgs = import nixpkgs { inherit system; };
        in pkgs.stdenv.mkDerivation {
          pname = "redshift-proofs";
          version = "0.3.0";
          src = ./.;
          nativeBuildInputs = [
            pkgs.lean4
            (pkgs.python3.withPackages (ps: [ ps.sympy ps.numpy ps.scipy ps.pandas ps.astropy ]))
          ];
          buildPhase = ''
            export HOME=$TMPDIR
            lake build
            sh check/run_all.sh ${if slow then "--slow" else ""}
          '';
          installPhase = "mkdir -p $out; echo ok > $out/VERIFIED";
        };
    in {
      packages = nixpkgs.lib.genAttrs systems (system: {
        default = build system false;
        full = build system true;
      });
      devShells = nixpkgs.lib.genAttrs systems (system:
        let pkgs = import nixpkgs { inherit system; };
        in {
          default = pkgs.mkShell {
            packages = [ pkgs.lean4 (pkgs.python3.withPackages (ps: [ ps.sympy ps.numpy ps.scipy ps.pandas ps.astropy ])) ];
          };
        });
    };
}
