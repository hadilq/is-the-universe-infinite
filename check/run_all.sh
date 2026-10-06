#!/usr/bin/env sh
# Runs every check.  Usage: sh check/run_all.sh [--slow] [PANTHEON_DATA_DIR]
#   --slow  also recomputes the geodesic-monism field equations in SymPy (~20 min);
#           Lean already proves them (GeodesicMonismAction.lean).
#   PANTHEON_DATA_DIR  directory with Pantheon+SH0ES.dat and Pantheon+SH0ES_STAT+SYS.cov;
#           the data fit is skipped when it is not given.
set -e
cd "$(dirname "$0")"
SLOW=
DATA=
for a in "$@"; do
  case "$a" in
    --slow) SLOW=--slow ;;
    *) DATA="$a" ;;
  esac
done
echo "== metrics solve their field equations";        python3 solutions.py $SLOW
echo "== redshift coefficients, both patches";        python3 redshift.py
echo "== two neighbouring rays (Lie transport)";      python3 lie_transport.py
echo "== z(dx, gamma) against the exact redshift";     python3 expansion_check.py
echo "== rotating 1+4 solution: derivation";            python3 rotating.py
echo "== the Hubble law of the non-Ricci-flat modes";   python3 hubble.py
echo "== three-sphere rotating along phi: solution, Ricci, patches";  python3 sphere.py
echo "== three-sphere rotating along phi: orbits and leaks";           python3 sphere_orbits.py
echo "== three-sphere rotating equally in both planes: derivation";   python3 sixd.py
echo "== three-sphere rotating equally in both planes: walls";        python3 trapping.py
if [ -n "$DATA" ]; then
  echo "== Pantheon+ fits";                           python3 pantheon_fit.py "$DATA"
  echo "== the three-sphere against all data";        python3 sphere_fit.py "$DATA"
fi
echo "ALL CHECKS PASSED"
