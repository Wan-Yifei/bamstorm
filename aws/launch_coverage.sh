#!/bin/bash
# Launch one EC2 instance that runs the bamstorm vs pysam coverage benchmark
# on a single contig and self-terminates.
#
# Usage:
#   ./launch_coverage.sh -i <ecr-image-uri> [-c chr22] [-T 600]
set -euo pipefail

export MSYS_NO_PATHCONV=1

PROFILE="${AWS_PROFILE:-admin}"
REGION="${AWS_REGION:-us-east-1}"
INSTANCE_TYPE="i4i.4xlarge"
ECR_IMAGE_URI=""
S3_BUCKET="dfci-bioinformatics-dev"
S3_BAM_KEY="test_15gb.bam"
S3_BAI_KEY="test_15gb.bam.bai"
S3_RESULTS_PREFIX="bamstorm-coverage-results"
INSTANCE_PROFILE="bamstorm-bench-ec2-profile"
CONTIG="chr4"
START=""
STOP=""
TIMEOUT=3600

usage() {
    cat <<EOF
Usage: $0 -i <ecr-image-uri> [OPTIONS]

  -i URI        ECR image URI (required; printed by cfn-deploy.sh)
  -t TYPE       Instance type, must have local NVMe (default: i4i.4xlarge)
  -c CONTIG     Contig/chromosome name (default: chr4)
  -s START      0-based region start, optional (default: full contig)
  -e END        0-based region end,   optional (default: full contig)
  -T SECONDS    pysam pileup hard timeout per run (default: 3600)
  -p PROFILE    AWS CLI profile (default: admin)
  -r REGION     AWS region (default: us-east-1)
  -h            Show this help
EOF
}

while getopts "i:t:c:s:e:T:p:r:h" opt; do
    case "$opt" in
        i) ECR_IMAGE_URI="$OPTARG" ;;
        t) INSTANCE_TYPE="$OPTARG" ;;
        c) CONTIG="$OPTARG" ;;
        s) START="$OPTARG" ;;
        e) STOP="$OPTARG" ;;
        T) TIMEOUT="$OPTARG" ;;
        p) PROFILE="$OPTARG" ;;
        r) REGION="$OPTARG" ;;
        h) usage; exit 0 ;;
        *) usage; exit 1 ;;
    esac
done

if [[ -z "$ECR_IMAGE_URI" ]]; then
    echo "ERROR: -i <ecr-image-uri> is required (run cfn-deploy.sh first)" >&2
    usage
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
winpath() { cygpath -w "$1" 2>/dev/null || echo "$1"; }

AMI_ID=$(aws ssm get-parameters --profile "$PROFILE" --region "$REGION" \
    --names /aws/service/canonical/ubuntu/server/22.04/stable/current/amd64/hvm/ebs-gp2/ami-id \
    --query 'Parameters[0].Value' --output text)

TS=$(date -u +%Y%m%dT%H%M%SZ)
# Build a region label for the result prefix (e.g. chr4 or chr4_0_10000000)
if [[ -n "$START" && -n "$STOP" ]]; then
    REGION_LABEL="${CONTIG}_${START}_${STOP}"
else
    REGION_LABEL="${CONTIG}"
fi
RESULT_PREFIX="${S3_RESULTS_PREFIX}/${TS}-${INSTANCE_TYPE}-${REGION_LABEL}"

# Build optional --start / --stop flags for bench_coverage.py
START_STOP_ARGS=""
[[ -n "$START" ]] && START_STOP_ARGS="$START_STOP_ARGS --start $START"
[[ -n "$STOP"  ]] && START_STOP_ARGS="$START_STOP_ARGS --stop $STOP"

USER_DATA="$SCRIPT_DIR/.tmp-coverage-user-data-${TS}.sh"
trap 'rm -f "$USER_DATA"' EXIT
sed \
    -e "s|__REGION__|${REGION}|g" \
    -e "s|__S3_BUCKET__|${S3_BUCKET}|g" \
    -e "s|__S3_BAM_KEY__|${S3_BAM_KEY}|g" \
    -e "s|__S3_BAI_KEY__|${S3_BAI_KEY}|g" \
    -e "s|__RESULT_PREFIX__|${RESULT_PREFIX}|g" \
    -e "s|__ECR_IMAGE_URI__|${ECR_IMAGE_URI}|g" \
    -e "s|__CONTIG__|${CONTIG}|g" \
    -e "s|__TIMEOUT__|${TIMEOUT}|g" \
    -e "s|__START_STOP_ARGS__|${START_STOP_ARGS}|g" \
    "$SCRIPT_DIR/coverage-user-data.sh.tmpl" > "$USER_DATA"

INSTANCE_ID=$(aws ec2 run-instances \
    --profile "$PROFILE" --region "$REGION" \
    --image-id "$AMI_ID" \
    --instance-type "$INSTANCE_TYPE" \
    --iam-instance-profile "Name=${INSTANCE_PROFILE}" \
    --instance-initiated-shutdown-behavior terminate \
    --user-data "file://$(winpath "$USER_DATA")" \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=bamstorm-coverage-${REGION_LABEL}},{Key=project,Value=bamstorm-bench}]" \
    --query 'Instances[0].InstanceId' --output text)

echo "Launched ${INSTANCE_ID} (${INSTANCE_TYPE})"
echo "Region   : ${REGION_LABEL}  (pysam timeout: ${TIMEOUT}s per run)"
echo "Results  : s3://${S3_BUCKET}/${RESULT_PREFIX}/"
echo ""
echo "Watch log (available ~2 min after launch):"
echo "  aws ec2 get-console-output --profile $PROFILE --region $REGION --instance-id $INSTANCE_ID --output text"
echo ""
echo "Fetch results after completion (~35-45 min):"
echo "  aws s3 sync s3://${S3_BUCKET}/${RESULT_PREFIX}/ bench/result/coverage/ --profile $PROFILE"
echo ""
echo "Instance self-terminates when benchmark finishes."
