#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/../00-common/00-library.sh"

"$BENCHMARK_ENV/bin/python" "$SCRIPT_DIR/00-build-comparability-masks.py" \
  --root "$BENCHMARK_ROOT"
