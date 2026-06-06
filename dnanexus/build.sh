#!/bin/bash
# Build the Docker image, push to a registry, then upload the applet to DNAnexus.
#
# Usage:
#   ./build.sh -r docker.io/yourname -t latest [-p project-xxxx] [-d /path/in/dx]
#
# Prerequisites:
#   docker login <registry>
#   dx login && dx select <project>
set -euo pipefail

REGISTRY="${REGISTRY:-docker.io/$(whoami)}"
TAG="${TAG:-latest}"
DX_PROJECT="${DX_PROJECT:-}"
DX_DEST="${DX_DEST:-/}"

usage() {
    cat <<EOF
Usage: $0 [OPTIONS]

  -r REGISTRY   Docker registry prefix  (default: docker.io/\$(whoami))
  -t TAG        Image tag               (default: latest)
  -p PROJECT    DNAnexus project ID     (default: currently selected project)
  -d DEST       Destination path in DX  (default: /)
  -h            Show this help

Steps:
  1. docker build  (project root Dockerfile)
  2. docker push
  3. dx build      (applet upload)
EOF
}

while getopts "r:t:p:d:h" opt; do
    case "$opt" in
        r) REGISTRY="$OPTARG" ;;
        t) TAG="$OPTARG" ;;
        p) DX_PROJECT="$OPTARG" ;;
        d) DX_DEST="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

IMAGE="$REGISTRY/bamstorm-bench:$TAG"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

echo "=== 1/3  Build Docker image ==="
docker build -t "$IMAGE" "$PROJECT_ROOT"
echo "Built: $IMAGE"

echo ""
echo "=== 2/3  Push to registry ==="
docker push "$IMAGE"

echo ""
echo "=== 3/3  Build DNAnexus applet ==="
DX_ARGS=()
if [[ -n "$DX_PROJECT" ]]; then
    DX_ARGS+=(--destination "${DX_PROJECT}:${DX_DEST}")
fi
dx build -f "$SCRIPT_DIR" "${DX_ARGS[@]+"${DX_ARGS[@]}"}"

echo ""
echo "Build complete."
echo ""
echo "To run a benchmark job:"
echo ""
echo "  dx run bamstorm_bench \\"
echo "    -i bam_file=<file-id>        \\"
echo "    -i bai_file=<file-id>        \\"
echo "    -i docker_image='$IMAGE' \\"
echo "    -i threads='2,4,8,16,32,64' \\"
echo "    -i repeats=3                 \\"
echo "    --instance-type mem1_ssd1_v2_x16"
echo ""
echo "See run_sweep.sh to submit multiple instance types in one shot."
