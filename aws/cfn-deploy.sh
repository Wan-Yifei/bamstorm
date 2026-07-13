#!/bin/bash
# Deploy the bamstorm benchmark stack via CloudFormation, then build and push
# the Docker image to the newly created ECR repo.
# CloudFormation handles IAM + ECR creation; Docker build/push is done here
# because CFN cannot run local shell commands.
#
# Usage:
#   ./cfn-deploy.sh [-p admin] [-r us-east-1] [-s bamstorm-bench-stack]
set -euo pipefail

export MSYS_NO_PATHCONV=1

# On Windows, Docker Desktop may not be injected into the Git Bash PATH.
# Add the known Docker Desktop bin directory as a fallback.
if ! command -v docker &>/dev/null; then
    export PATH="$PATH:/c/Program Files/Docker/Docker/resources/bin"
fi

PROFILE="${AWS_PROFILE:-admin}"
REGION="${AWS_REGION:-us-east-1}"
STACK_NAME="bamstorm-bench-stack"
TEMPLATE_FILE="$(cd "$(dirname "$0")" && pwd)/cloudformation.yaml"
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

while getopts "p:r:s:h" opt; do
    case "$opt" in
        p) PROFILE="$OPTARG" ;;
        r) REGION="$OPTARG" ;;
        s) STACK_NAME="$OPTARG" ;;
        h) echo "Usage: $0 [-p profile] [-r region] [-s stack-name]"; exit 0 ;;
        *) exit 1 ;;
    esac
done

echo "=== 1/3  CloudFormation stack: ${STACK_NAME} ==="
STACK_STATUS=$(aws cloudformation describe-stacks \
    --profile "$PROFILE" --region "$REGION" \
    --stack-name "$STACK_NAME" \
    --query 'Stacks[0].StackStatus' --output text 2>/dev/null || echo "DOES_NOT_EXIST")

if [[ "$STACK_STATUS" == "DOES_NOT_EXIST" ]]; then
    echo "Creating stack..."
    aws cloudformation create-stack \
        --profile "$PROFILE" --region "$REGION" \
        --stack-name "$STACK_NAME" \
        --template-body "file://$(cygpath -w "$TEMPLATE_FILE" 2>/dev/null || echo "$TEMPLATE_FILE")" \
        --capabilities CAPABILITY_NAMED_IAM
    echo "Waiting for CREATE_COMPLETE..."
    aws cloudformation wait stack-create-complete \
        --profile "$PROFILE" --region "$REGION" --stack-name "$STACK_NAME"
elif [[ "$STACK_STATUS" == *"ROLLBACK"* ]] || [[ "$STACK_STATUS" == "DELETE_COMPLETE" ]]; then
    echo "Stack is in ${STACK_STATUS} — delete it first and re-run."
    exit 1
else
    echo "Stack exists (${STACK_STATUS}), updating..."
    UPDATE_OUTPUT=$(aws cloudformation update-stack \
        --profile "$PROFILE" --region "$REGION" \
        --stack-name "$STACK_NAME" \
        --template-body "file://$(cygpath -w "$TEMPLATE_FILE" 2>/dev/null || echo "$TEMPLATE_FILE")" \
        --capabilities CAPABILITY_NAMED_IAM 2>&1 || true)
    if echo "$UPDATE_OUTPUT" | grep -q "No updates are to be performed"; then
        echo "No infrastructure changes."
    else
        echo "Waiting for UPDATE_COMPLETE..."
        aws cloudformation wait stack-update-complete \
            --profile "$PROFILE" --region "$REGION" --stack-name "$STACK_NAME"
    fi
fi

ECR_URI=$(aws cloudformation describe-stacks \
    --profile "$PROFILE" --region "$REGION" \
    --stack-name "$STACK_NAME" \
    --query "Stacks[0].Outputs[?OutputKey=='ECRRepositoryUri'].OutputValue" \
    --output text)
echo "ECR image URI: ${ECR_URI}"

echo ""
echo "=== 2/3  Docker login + build + push ==="
ACCOUNT_ID=$(aws sts get-caller-identity --profile "$PROFILE" --query Account --output text)
aws ecr get-login-password --profile "$PROFILE" --region "$REGION" \
    | docker login --username AWS --password-stdin "${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com"
docker buildx build --platform linux/amd64 -t "$ECR_URI" --push \
    "$(cygpath -w "$PROJECT_ROOT" 2>/dev/null || echo "$PROJECT_ROOT")"

echo ""
echo "=== 3/3  Done ==="
echo ""
echo "Stack   : ${STACK_NAME}"
echo "ECR     : ${ECR_URI}"
echo ""
echo "Next:"
echo "  ./launch.sh    -i ${ECR_URI} -t i4i.4xlarge"
echo "  ./run_sweep.sh ${ECR_URI}"
echo ""
echo "Tear down:"
echo "  aws cloudformation delete-stack --profile ${PROFILE} --region ${REGION} --stack-name ${STACK_NAME}"
