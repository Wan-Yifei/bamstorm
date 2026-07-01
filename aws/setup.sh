#!/bin/bash
# One-time AWS setup for the bamstorm benchmark: ECR repo + pushed image,
# IAM role/instance-profile scoped to the test S3 bucket + this ECR repo.
#
# Usage:
#   ./setup.sh [-p admin] [-r us-east-1] [-b dfci-bioinformatics-dev]
set -euo pipefail

PROFILE="${AWS_PROFILE:-admin}"
REGION="${AWS_REGION:-us-east-1}"
S3_BUCKET="dfci-bioinformatics-dev"
ECR_REPO="bamstorm-bench"
ROLE_NAME="bamstorm-bench-ec2-role"
INSTANCE_PROFILE="bamstorm-bench-ec2-profile"

while getopts "p:r:b:h" opt; do
    case "$opt" in
        p) PROFILE="$OPTARG" ;;
        r) REGION="$OPTARG" ;;
        b) S3_BUCKET="$OPTARG" ;;
        h) echo "Usage: $0 [-p profile] [-r region] [-b s3-bucket]"; exit 0 ;;
        *) exit 1 ;;
    esac
done

ACCOUNT_ID=$(aws sts get-caller-identity --profile "$PROFILE" --query Account --output text)
ECR_URI="${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com/${ECR_REPO}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
TMP_DIR="$SCRIPT_DIR/.tmp"
mkdir -p "$TMP_DIR"
trap 'rm -rf "$TMP_DIR"' EXIT

# aws.exe (native Windows binary) can't resolve MSYS-style file:// paths,
# so route every --*-document file:// argument through cygpath -w.
winpath() { cygpath -w "$1" 2>/dev/null || echo "$1"; }

echo "=== 1/4  ECR repository ==="
if ! aws ecr describe-repositories --profile "$PROFILE" --region "$REGION" \
        --repository-names "$ECR_REPO" >/dev/null 2>&1; then
    aws ecr create-repository --profile "$PROFILE" --region "$REGION" --repository-name "$ECR_REPO"
fi
echo "Repo: $ECR_URI"

echo ""
echo "=== 2/4  Build + push image (linux/amd64) ==="
aws ecr get-login-password --profile "$PROFILE" --region "$REGION" \
    | docker login --username AWS --password-stdin "${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com"
docker buildx build --platform linux/amd64 -t "${ECR_URI}:latest" --push "$PROJECT_ROOT"

echo ""
echo "=== 3/4  IAM role + instance profile ==="
cat > "$TMP_DIR/trust.json" <<JSON
{
  "Version": "2012-10-17",
  "Statement": [{"Effect": "Allow", "Principal": {"Service": "ec2.amazonaws.com"}, "Action": "sts:AssumeRole"}]
}
JSON
if ! aws iam get-role --profile "$PROFILE" --role-name "$ROLE_NAME" >/dev/null 2>&1; then
    aws iam create-role --profile "$PROFILE" --role-name "$ROLE_NAME" \
        --assume-role-policy-document "file://$(winpath "$TMP_DIR/trust.json")"
fi

cat > "$TMP_DIR/policy.json" <<JSON
{
  "Version": "2012-10-17",
  "Statement": [
    {"Effect": "Allow", "Action": ["s3:GetObject"], "Resource": "arn:aws:s3:::${S3_BUCKET}/test_15gb.bam*"},
    {"Effect": "Allow", "Action": ["s3:PutObject"], "Resource": "arn:aws:s3:::${S3_BUCKET}/bamstorm-bench-results/*"},
    {"Effect": "Allow", "Action": ["s3:ListBucket"], "Resource": "arn:aws:s3:::${S3_BUCKET}"},
    {"Effect": "Allow", "Action": ["ecr:GetAuthorizationToken"], "Resource": "*"},
    {"Effect": "Allow", "Action": ["ecr:BatchGetImage", "ecr:GetDownloadUrlForLayer", "ecr:BatchCheckLayerAvailability"],
     "Resource": "arn:aws:ecr:${REGION}:${ACCOUNT_ID}:repository/${ECR_REPO}"}
  ]
}
JSON
aws iam put-role-policy --profile "$PROFILE" --role-name "$ROLE_NAME" \
    --policy-name bamstorm-bench-policy --policy-document "file://$(winpath "$TMP_DIR/policy.json")"

if ! aws iam get-instance-profile --profile "$PROFILE" --instance-profile-name "$INSTANCE_PROFILE" >/dev/null 2>&1; then
    aws iam create-instance-profile --profile "$PROFILE" --instance-profile-name "$INSTANCE_PROFILE"
    aws iam add-role-to-instance-profile --profile "$PROFILE" \
        --instance-profile-name "$INSTANCE_PROFILE" --role-name "$ROLE_NAME"
    echo "Waiting for instance profile propagation..."
    sleep 15
fi

echo ""
echo "=== 4/4  Done ==="
echo "ECR image        : ${ECR_URI}:latest"
echo "Instance profile : ${INSTANCE_PROFILE}"
echo ""
echo "Next:"
echo "  ./launch.sh -t i4i.4xlarge -i ${ECR_URI}:latest"
echo "  ./run_sweep.sh ${ECR_URI}:latest"
