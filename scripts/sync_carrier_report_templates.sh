#!/usr/bin/env bash
# Copy Jinja carrier PDF templates into this repo's data/carrier_report/.
# Your Carrier_result project may use a folder named "templates" (not "carrier_report"):
#   e.g. /path/to/Carrier_result/templates/
#
# Requires read permission on the source files (chmod a+r … if needed).
#
# Usage:
#   ./scripts/sync_carrier_report_templates.sh /path/to/templates_or_carrier_report
#   CARRIER_REPORT_TEMPLATE_SOURCE=/path/to/templates ./scripts/sync_carrier_report_templates.sh

set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="${1:-${CARRIER_REPORT_TEMPLATE_SOURCE:-}}"
if [[ -z "$SRC" ]]; then
  echo "Usage: $0 /path/to/carrier_report" >&2
  echo "Or: CARRIER_REPORT_TEMPLATE_SOURCE=/path/to/carrier_report $0" >&2
  exit 1
fi
if [[ ! -d "$SRC" ]]; then
  echo "Not a directory: $SRC" >&2
  exit 1
fi
DST="$ROOT/data/carrier_report"
mkdir -p "$DST"
count=0
shopt -s nullglob
for f in "$SRC"/*.html; do
  cp -v "$f" "$DST/"
  count=$((count + 1))
done
if [[ "$count" -eq 0 ]]; then
  echo "No *.html files found in $SRC" >&2
  exit 1
fi
echo "Copied $count file(s) into $DST"
