#!/bin/bash
# Launch one bamstorm benchmark instance per instance type, in parallel,
# to compare throughput/cost across sizes/families.
#
# Usage:
#   ./run_sweep.sh <ecr-image-uri> [instance-type ...]
#
# Default sweep: i4i.2xlarge i4i.4xlarge i4i.8xlarge i3en.3xlarge
set -euo pipefail

ECR_IMAGE_URI="${1:?usage: run_sweep.sh <ecr-image-uri> [instance-type ...]}"
shift || true
TYPES=("$@")
if [[ ${#TYPES[@]} -eq 0 ]]; then
    TYPES=(i4i.2xlarge i4i.4xlarge i4i.8xlarge i3en.3xlarge)
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
for T in "${TYPES[@]}"; do
    "$SCRIPT_DIR/launch.sh" -i "$ECR_IMAGE_URI" -t "$T"
done
