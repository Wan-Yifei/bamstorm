#!/bin/bash
# Sync benchmark results from S3 to a local directory for plotting.
#
# Usage:
#   ./fetch_results.sh [dest-dir]
set -euo pipefail

PROFILE="${AWS_PROFILE:-admin}"
REGION="${AWS_REGION:-us-east-1}"
S3_BUCKET="dfci-bioinformatics-dev"
DEST="${1:-$(cd "$(dirname "$0")/.." && pwd)/bench/result/aws}"

mkdir -p "$DEST"
aws s3 sync "s3://${S3_BUCKET}/bamstorm-bench-results/" "$DEST" --profile "$PROFILE" --region "$REGION"

echo "Synced to $DEST"
echo "Plot one run with:"
echo "  python3 bench/plot_report.py \"$DEST/<run-dir>/benchmark.csv\""
