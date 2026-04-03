#!/usr/bin/env bash
# Rename carrier_screening artifact folders (<work>/<old_sample> → <work>/<new_sample>)
# and files whose basename starts with the old sample name (e.g. "dark genes 2.cleaned.vcf").
#
# Usage:
#   ./scripts/carrier-rename-sample-workdir.sh <work_dir> <old_sample> <new_sample> [work_root]
#
# Example:
#   ./scripts/carrier-rename-sample-workdir.sh 2601 "dark genes 2" "dark_genes_2"
#
# Optional:
#   DRY_RUN=1  — print actions only
#
set -euo pipefail

DRY_RUN="${DRY_RUN:-0}"

die() { echo "ERROR: $*" >&2; exit 1; }

[[ $# -ge 3 ]] || die "usage: $0 <work_dir> <old_sample> <new_sample> [work_root]"

WORK_DIR="$1"
OLD_SAMPLE="$2"
NEW_SAMPLE="$3"
WORK_ROOT="${4:-/data/gx-exome-work}"

[[ -n "$WORK_DIR" ]] || die "work_dir empty"
[[ -n "$OLD_SAMPLE" ]] || die "old_sample empty"
[[ -n "$NEW_SAMPLE" ]] || die "new_sample empty"
[[ "$OLD_SAMPLE" == "$NEW_SAMPLE" ]] && die "old and new sample names are the same"

if [[ "$NEW_SAMPLE" =~ [[:space:]] ]]; then
  echo "WARN: new_sample contains whitespace — pipeline paths may break; prefer underscores." >&2
fi

rename_prefixed_files() {
  local dir="$1"
  [[ -d "$dir" ]] || return 0
  local old="$OLD_SAMPLE"
  local new="$NEW_SAMPLE"
  # shellcheck disable=SC2016
  find "$dir" -type f -name "${old}*" \( ! -path '*/.git/*' \) 2>/dev/null | sort -r | while IFS= read -r f; do
    [[ -f "$f" ]] || continue
    base=$(basename -- "$f")
    parent=$(dirname -- "$f")
    suffix="${base#"${old}"}"
    [[ "$suffix" != "$base" ]] || continue
    dest="${parent}/${new}${suffix}"
    if [[ -e "$dest" ]]; then
      echo "SKIP (target exists): $f → $dest" >&2
      continue
    fi
    if [[ "$DRY_RUN" == "1" ]]; then
      echo "FILE  $f → $dest"
    else
      mv -- "$f" "$dest"
      echo "FILE  $f → $dest"
    fi
  done
}

rename_topdir() {
  local kind="$1"
  local old_path="${WORK_ROOT}/${kind}/${WORK_DIR}/${OLD_SAMPLE}"
  local new_path="${WORK_ROOT}/${kind}/${WORK_DIR}/${NEW_SAMPLE}"

  [[ -d "$old_path" ]] || return 0

  if [[ -e "$new_path" ]] || [[ -L "$new_path" ]]; then
    die "${kind}: target already exists: $new_path"
  fi

  # Rename per-sample files first while the directory still uses the old name (stable paths).
  rename_prefixed_files "$old_path"

  if [[ "$DRY_RUN" == "1" ]]; then
    echo "DIR   $old_path → $new_path"
  else
    mv -- "$old_path" "$new_path"
    echo "DIR   $old_path → $new_path"
  fi
}

echo "work_root=$WORK_ROOT work_dir=$WORK_DIR"
echo "old_sample=$OLD_SAMPLE"
echo "new_sample=$NEW_SAMPLE"
[[ "$DRY_RUN" == "1" ]] && echo "DRY_RUN=1 (no changes)"

for kind in fastq analysis output log; do
  rename_topdir "$kind"
done

echo "Done. Update the daemon order: sample_name (and order_id if you use the same label) must equal '$NEW_SAMPLE' for paths to match."
