#!/bin/bash
# Validate `anarci --pair` (structural VH/VL pairing + antibody-format inference) against
# a curated set of real PDB structures spanning formats easy→hard. Ground truth is from
# RCSB entry metadata + the deposition's own chain labels; the exotic-format IDs were
# format-verified against the RCSB data API.
#
# Build first:  cargo build --release -p anarci-cli
# Run:          scripts/validate_pairing.sh
set -u
BIN="${BIN:-$(dirname "$0")/../target/release/anarci}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

# id | expected format substring | description (ground truth)
cases=(
  "7FAB|Fab/Fv|Fab — one VH–VL pair on separate chains"
  "1DZB|scFv|scFv (VH+VL on one chain), 2 copies + antigen"
  "1MOE|diabody|diabody — cross-paired"
  "9LSH|diabody|bispecific diabody, deposition-labeled cross-pairing"
  "6VUO|diabody|diabody"
  "1MEL|VHH|camelid VHH + lysozyme (the original nanobody)"
  "5G5R|VHH|llama nanobody"
  "21EN|triabody|triabody — cyclic 3-Fv cross-pairing"
  "4HJJ|DVD-Ig|DVD-Ig — stacked dual-variable, 2 Fvs"
  "6O8D|CODV|CODV — crossover dual-variable (naive-pairer trap)"
  "7KN6|Fab/Fv + VHH|Fab + unpaired nanobody in one file"
  "8CII|Fab/Fv + VHH|two Fabs + one nanobody (5 V domains)"
  "8IDM|VHH|two independent nanobodies (must NOT pair to each other)"
)

pass=0; total=0
for c in "${cases[@]}"; do
  IFS='|' read -r id exp desc <<< "$c"
  total=$((total + 1))
  f="$TMP/$id.pdb"
  curl -sf "https://files.rcsb.org/download/$id.pdb" -o "$f" 2>/dev/null ||
    { curl -sf "https://files.rcsb.org/download/$id.cif" -o "$TMP/$id.cif" 2>/dev/null && f="$TMP/$id.cif"; }
  got=$("$BIN" --pair "$f" 2>/dev/null | head -1 | sed 's/format: //')
  if [[ "$got" == *"$exp"* ]]; then ok="PASS"; pass=$((pass + 1)); else ok="FAIL"; fi
  printf "%-4s %-5s  %-38s %s\n" "$ok" "$id" "${got:-<none>}" "$desc"
done

echo "-- antigen detection (--pair reports the chain(s) each Fv/VHH binds) --"
# id | grep pattern (or __none__) | description
ag_cases=(
  "1MEL|binds antigen: 1MEL_L|VHH → lysozyme (chain L)"
  "1DZB|binds antigen: 1DZB_X|scFv → turkey lysozyme (chain X)"
  "6MFT|binds antigen|Fab → single antigen chain"
  "5FUU|multi-subunit epitope|PGT151 → quaternary epitope (≥2 Env chains)"
  "7FAB|__none__|Fab alone → no antigen in file"
  "1DVF|anti-idiotype|anti-idiotype Fv–Fv (E5.2 binds D1.3)"
  "8CII|__not:anti-idiotype|3 antibodies, same antigen, NO spurious anti-idiotype"
)
for c in "${ag_cases[@]}"; do
  IFS='|' read -r id pat desc <<< "$c"
  total=$((total + 1))
  f="$TMP/$id.pdb"; [ -f "$f" ] || f="$TMP/$id.cif"
  if [ ! -f "$f" ]; then
    curl -sf "https://files.rcsb.org/download/$id.pdb" -o "$TMP/$id.pdb" 2>/dev/null && f="$TMP/$id.pdb" ||
      { curl -sf "https://files.rcsb.org/download/$id.cif" -o "$TMP/$id.cif" 2>/dev/null && f="$TMP/$id.cif"; }
  fi
  out=$("$BIN" --pair "$f" 2>/dev/null)
  if [ "$pat" = "__none__" ]; then
    if echo "$out" | grep -q "binds"; then ok="FAIL"; else ok="PASS"; pass=$((pass + 1)); fi
  elif [[ "$pat" == __not:* ]]; then
    if echo "$out" | grep -q "${pat#__not:}"; then ok="FAIL"; else ok="PASS"; pass=$((pass + 1)); fi
  else
    if echo "$out" | grep -q "$pat"; then ok="PASS"; pass=$((pass + 1)); else ok="FAIL"; fi
  fi
  printf "%-4s %-5s  %s\n" "$ok" "$id" "$desc"
done

echo "----"
echo "$pass/$total correct"
