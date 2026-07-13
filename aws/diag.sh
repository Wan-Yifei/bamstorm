#!/bin/bash
# Standalone diagnostic instance: verify (1) real local-NVMe parallel
# bandwidth with a properly-configured fio, and (2) whether
# /proc/sys/vm/drop_caches actually evicts the page cache (both directly on
# the host and via a --privileged container, matching bench.py's code path).
#
# Cheaper/faster than a full launch.sh run -- no Docker image pull, no
# thread sweep. Reuses the IAM instance profile created by setup.sh.
#
# Usage:
#   ./diag.sh [-t i4i.4xlarge] [-p admin] [-r us-east-1]
set -euo pipefail

# Git Bash (MSYS) rewrites leading-slash args (SSM parameter names) into
# Windows paths unless this is set.
export MSYS_NO_PATHCONV=1

PROFILE="${AWS_PROFILE:-admin}"
REGION="${AWS_REGION:-us-east-1}"
INSTANCE_TYPE="i4i.4xlarge"
S3_BUCKET="dfci-bioinformatics-dev"
S3_BAM_KEY="test_15gb.bam"
S3_RESULTS_PREFIX="bamstorm-bench-results"
INSTANCE_PROFILE="bamstorm-bench-ec2-profile"

usage() {
    cat <<EOF
Usage: $0 [OPTIONS]

  -t TYPE       Instance type, must have local NVMe instance store (default: i4i.4xlarge)
  -p PROFILE    AWS CLI profile (default: admin)
  -r REGION     AWS region (default: us-east-1)
  -h            Show this help
EOF
}

while getopts "t:p:r:h" opt; do
    case "$opt" in
        t) INSTANCE_TYPE="$OPTARG" ;;
        p) PROFILE="$OPTARG" ;;
        r) REGION="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# aws.exe (native Windows binary) can't resolve MSYS-style file:// paths,
# so route every file:// argument through cygpath -w.
winpath() { cygpath -w "$1" 2>/dev/null || echo "$1"; }

AMI_ID=$(aws ssm get-parameters --profile "$PROFILE" --region "$REGION" \
    --names /aws/service/canonical/ubuntu/server/22.04/stable/current/amd64/hvm/ebs-gp2/ami-id \
    --query 'Parameters[0].Value' --output text)

TS=$(date -u +%Y%m%dT%H%M%SZ)
RESULT_PREFIX="${S3_RESULTS_PREFIX}/${TS}-diag-${INSTANCE_TYPE}"

USER_DATA="$SCRIPT_DIR/.tmp-diag-user-data-${TS}.sh"
trap 'rm -f "$USER_DATA"' EXIT
sed \
    -e "s|__REGION__|${REGION}|g" \
    -e "s|__S3_BUCKET__|${S3_BUCKET}|g" \
    -e "s|__S3_BAM_KEY__|${S3_BAM_KEY}|g" \
    -e "s|__RESULT_PREFIX__|${RESULT_PREFIX}|g" \
    "$SCRIPT_DIR/diag-user-data.sh.tmpl" > "$USER_DATA"

INSTANCE_ID=$(aws ec2 run-instances \
    --profile "$PROFILE" --region "$REGION" \
    --image-id "$AMI_ID" \
    --instance-type "$INSTANCE_TYPE" \
    --iam-instance-profile "Name=${INSTANCE_PROFILE}" \
    --instance-initiated-shutdown-behavior terminate \
    --user-data "file://$(winpath "$USER_DATA")" \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=bamstorm-diag-${INSTANCE_TYPE}},{Key=project,Value=bamstorm-bench}]" \
    --query 'Instances[0].InstanceId' --output text)

echo "Launched ${INSTANCE_ID} (${INSTANCE_TYPE}) -- diagnostic run"
echo "Results will land at: s3://${S3_BUCKET}/${RESULT_PREFIX}/"
echo ""
echo "Watch boot/diag log:"
echo "  aws ec2 get-console-output --profile $PROFILE --region $REGION --instance-id $INSTANCE_ID --output text"
echo ""
echo "Instance self-terminates when the diagnostic finishes (~5-10 min)."
