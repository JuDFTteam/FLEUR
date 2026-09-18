#!/usr/bin/env bash
#
# Run the FLEUR metaGGA (SCAN) validation suite and compare the resulting
# magnetic moments against the SCAN reference values in reference_values.csv
# (Tran et al., Phys. Rev. B 102, 024407 (2020), Tables III and VI).
#
# Usage:  ./run_mgga_tests.sh [options] [case ...]
# See     ./run_mgga_tests.sh --help
#
# Requires a FLEUR binary built WITH libxc: SCAN is only available through the
# libxc interface, there is no in-built metaGGA in FLEUR.  Build one with e.g.
#     ./configure.sh -libxc TRUE -make <label>
#
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"

FLEUR=""
WORKDIR="$HERE/work"
COMPARE_ONLY=0
KEEP=0
declare -a WANTED=()

usage() {
   cat <<EOF
Run the FLEUR metaGGA (SCAN) validation suite.

Usage: $(basename "$0") [options] [case ...]

  case            Substring of a case path to select, e.g. 'Fe', 'fm/', 'af/CrSb2'.
                  May be repeated.  Default: every case that has an inp.xml.

Options:
  -b, --binary PATH   FLEUR executable to use.  Default: autodetected from
                      \$REPO/build*/fleur, preferring one built with libxc.
  -w, --workdir DIR   Directory to run in.  Default: $WORKDIR
  -c, --compare-only  Do not run FLEUR; only re-evaluate existing results in
                      the work directory.
  -k, --keep          Keep an existing work directory (default: wipe per case
                      before running, so a stale cdn does not leak in).
  -h, --help          This message.

Exit status is 0 if every selected case ran to completion, 1 otherwise.
Agreement with the paper is reported but never fails the run: these are
physics reference values, not regression baselines (see README.md).
EOF
}

while [[ $# -gt 0 ]]; do
   case "$1" in
      -b|--binary)       FLEUR="$2"; shift 2 ;;
      -w|--workdir)      WORKDIR="$2"; shift 2 ;;
      -c|--compare-only) COMPARE_ONLY=1; shift ;;
      -k|--keep)         KEEP=1; shift ;;
      -h|--help)         usage; exit 0 ;;
      -*)                echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
      *)                 WANTED+=("$1"); shift ;;
   esac
done

# --- locate a python3 (used for XML/CSV parsing only) ----------------------
PYTHON="${juDFT_PYTHON:-}"
if [[ -z "$PYTHON" ]]; then
   for c in python3 python; do
      if command -v "$c" >/dev/null 2>&1; then PYTHON="$(command -v "$c")"; break; fi
   done
fi
if [[ -z "$PYTHON" ]]; then
   echo "error: no python3 found; set juDFT_PYTHON to one" >&2; exit 2
fi

# --- locate the FLEUR binary ----------------------------------------------
# There is no reliable way to tell from the build tree whether libxc is usable
# (FLEUR_USE_LIBXC in CMakeCache.txt does not reflect reality, and CPP_LIBXC is
# defined in every build), so just take the newest non-debug binary and print
# its version - a stale binary is the usual cause of a schema-validation abort.
if [[ $COMPARE_ONLY -eq 0 ]]; then
   if [[ -z "$FLEUR" ]]; then
      newest=""
      for cand in "$REPO"/build*/fleur; do
         [[ -x "$cand" ]] || continue
         [[ "$cand" == *".debug/fleur" ]] && continue
         if [[ -z "$newest" || "$cand" -nt "$newest" ]]; then newest="$cand"; fi
      done
      FLEUR="$newest"
   fi
   if [[ -z "$FLEUR" || ! -x "$FLEUR" ]]; then
      echo "error: no FLEUR binary found. Pass one with -b, or build with:" >&2
      echo "       ./configure.sh -libxc TRUE -make <label>" >&2
      exit 2
   fi
   version=$(strings "$FLEUR" 2>/dev/null | grep -oE 'MaX-Release [0-9.]+' | head -1)
   echo "FLEUR binary : $FLEUR  (${version:-version unknown})"
   if ! strings "$FLEUR" 2>/dev/null | grep -qx 'mgga_x_scan'; then
      echo "WARNING: no libxc metaGGA symbols in this binary; SCAN will not work." >&2
      echo "         Rebuild with: ./configure.sh -libxc TRUE -make <label>" >&2
   fi
fi
echo "work dir     : $WORKDIR"
echo

# --- collect cases ---------------------------------------------------------
declare -a CASES=()
for d in "$HERE"/fm/*/ "$HERE"/af/*/; do
   [[ -d "$d" ]] || continue
   rel="${d#"$HERE"/}"; rel="${rel%/}"
   CASES+=("$rel")
done

if [[ ${#WANTED[@]} -gt 0 ]]; then
   declare -a SEL=()
   for w in "${WANTED[@]}"; do
      # exact match on the full path or on the material name wins, so that
      # 'Fe' selects fm/Fe and not also fm/FeCo; fall back to substring.
      exact=0
      for rel in "${CASES[@]}"; do
         if [[ "$rel" == "$w" || "$(basename "$rel")" == "$w" ]]; then
            SEL+=("$rel"); exact=1
         fi
      done
      if [[ $exact -eq 0 ]]; then
         for rel in "${CASES[@]}"; do
            [[ "$rel" == *"$w"* ]] && SEL+=("$rel")
         done
      fi
   done
   # de-duplicate, preserving order
   declare -a UNIQ=()
   for rel in "${SEL[@]}"; do
      dup=0
      for u in ${UNIQ[@]+"${UNIQ[@]}"}; do [[ "$u" == "$rel" ]] && dup=1; done
      [[ $dup -eq 0 ]] && UNIQ+=("$rel")
   done
   CASES=(${UNIQ[@]+"${UNIQ[@]}"})
fi

if [[ ${#CASES[@]} -eq 0 ]]; then
   echo "no cases selected" >&2; exit 2
fi

# --- run -------------------------------------------------------------------
rc=0
if [[ $COMPARE_ONLY -eq 0 ]]; then
   for rel in "${CASES[@]}"; do
      src="$HERE/$rel"
      if [[ ! -f "$src/inp.xml" ]]; then
         printf '%-12s SKIP  no inp.xml (run inpgen -f structure.inpgen first)\n' "$rel"
         continue
      fi
      dst="$WORKDIR/$rel"
      [[ $KEEP -eq 0 ]] && rm -rf "$dst"
      mkdir -p "$dst"
      for f in inp.xml kpts.xml sym.xml relax.xml; do
         [[ -f "$src/$f" ]] && cp "$src/$f" "$dst/"
      done
      printf '%-12s running ... ' "$rel"
      start=$SECONDS
      ( cd "$dst" && "$FLEUR" ) > "$dst/run.log" 2>&1
      status=$?
      took=$((SECONDS - start))
      if [[ $status -ne 0 ]]; then
         msg=$(grep -A1 'FLEUR-Error' "$dst/run.log" | sed -n '2s/Error message: //p' | head -1)
         printf 'FAILED (%ds) %s\n' "$took" "${msg:-see $dst/run.log}"
         rc=1
      else
         last=$(grep -oE 'Iteration: *[0-9]+ *Distance: *[0-9.]+' "$dst/run.log" | tail -1)
         printf 'done (%ds)  %s\n' "$took" "$last"
      fi
   done
   echo
fi

# --- compare ---------------------------------------------------------------
"$PYTHON" "$HERE/compare_moments.py" --workdir "$WORKDIR" --suite "$HERE" "${CASES[@]}"
cmp_rc=$?

[[ $cmp_rc -ne 0 ]] && rc=1
exit $rc
