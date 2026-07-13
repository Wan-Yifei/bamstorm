# Bamstorm AWS Benchmark

Runs the same `bench/bench.py` comparison (bamstorm vs samtools vs rabbitbam
vs pysam) on EC2 instead of DNAnexus. DNAnexus was abandoned for this
purpose because its workers have no root access — `/proc/sys/vm/drop_caches`
is read-only and `posix_fadvise(DONTNEED)` doesn't evict pages on its
369 GB-RAM machines, so cold-cache reads couldn't be measured. EC2 gives
real root, so `bench.py`'s `drop_caches()` now performs a genuine sysctl
drop (falling back to the old copy/fadvise trick only when root isn't
available).

## Architecture

- One instance type with local NVMe **instance store** (i4i/i3en family) —
  local disk avoids EBS network-storage variance and gives root access to
  truly evict the page cache.
- The instance downloads 4 isolated copies of the test BAM (one per tool,
  same isolation scheme as `dnanexus/src/code.sh`), runs the full
  `bench.toml` thread sweep via the existing Docker image, uploads
  `benchmark.csv` + logs to S3, then **self-terminates** (no SSH, no
  lingering cost).
- Account: `733024369092`, profile `admin`, region `us-east-1`.
- Test data: `s3://dfci-bioinformatics-dev/test_15gb.bam` (+ `.bai`).
- Results land under `s3://dfci-bioinformatics-dev/bamstorm-bench-results/<timestamp>-<instance-type>/`.

## File overview

| File | Purpose |
|---|---|
| `cloudformation.yaml` | IaC template — ECR repo + IAM role/instance-profile |
| `cfn-deploy.sh` | Deploy CFN stack then build + push Docker image |
| `setup.sh` | Same one-time setup via raw AWS CLI (no CFN state file) |
| `launch.sh` | Launch one benchmark instance |
| `run_sweep.sh` | Launch one instance per instance type in parallel |
| `diag.sh` | Lightweight diagnostic: verify disk bandwidth + drop_caches |
| `fetch_results.sh` | Sync results from S3 to `bench/result/aws/` |
| `user-data.sh.tmpl` | EC2 user-data script (templated by launch.sh) |
| `diag-user-data.sh.tmpl` | User-data for diag.sh |

---

## One-time setup

Two equivalent options — pick one.

### Option A: CloudFormation (recommended)

CloudFormation tracks all created resources in a stack and supports
one-command teardown. `cfn-deploy.sh` deploys the stack then builds and
pushes the Docker image (the part CFN cannot do).

```bash
./cfn-deploy.sh
```

This creates/updates the `bamstorm-bench-stack` CloudFormation stack
(ECR repo + IAM role + instance profile), then runs
`docker buildx build --push`. Idempotent — safe to re-run after code
changes to push a fresh image.

Tear down everything when done:

```bash
aws cloudformation delete-stack --profile admin --region us-east-1 \
  --stack-name bamstorm-bench-stack
```

### Option B: raw AWS CLI

```bash
./setup.sh
```

Idempotent CLI calls — creates the same ECR repo + IAM resources without a
CFN state file. Resources must be cleaned up individually (see Cleanup).

---

## Run benchmarks

```bash
# Single instance
./launch.sh -i <ecr-image-uri> -t i4i.4xlarge

# Multi-instance-type sweep (parallel)
./run_sweep.sh <ecr-image-uri> i4i.2xlarge i4i.4xlarge i4i.8xlarge i3en.3xlarge

# Lightweight diagnostic only (fio bandwidth + drop_caches check, ~5-10 min)
./diag.sh -t i4i.4xlarge
```

`<ecr-image-uri>` is printed by `cfn-deploy.sh` / `setup.sh`, or retrieved with:

```bash
export MSYS_NO_PATHCONV=1
aws cloudformation describe-stacks --profile admin --region us-east-1 \
  --stack-name bamstorm-bench-stack \
  --query "Stacks[0].Outputs[?OutputKey=='ECRRepositoryUri'].OutputValue" \
  --output text
```

Each benchmark instance takes ~20-40 min and self-terminates after uploading
results. Monitor progress:

```bash
aws ec2 get-console-output --profile admin --region us-east-1 \
  --instance-id <id> --latest --output text
```

## Pull + plot results

```bash
./fetch_results.sh
python3 ../bench/plot_report.py bench/result/aws/<run-dir>/benchmark.csv
```

## Cost

i4i.4xlarge on-demand is roughly $1.3-1.5/hr in us-east-1; a single ~30 min
run is well under $1. A 4-type sweep is a few dollars total. All instances
self-terminate after uploading results. The ECR repo and IAM resources have
no ongoing cost.

## Cleanup (Option B only — Option A uses `delete-stack` above)

```bash
aws iam remove-role-from-instance-profile --profile admin \
  --instance-profile-name bamstorm-bench-ec2-profile \
  --role-name bamstorm-bench-ec2-role
aws iam delete-instance-profile --profile admin \
  --instance-profile-name bamstorm-bench-ec2-profile
aws iam delete-role-policy --profile admin \
  --role-name bamstorm-bench-ec2-role --policy-name bamstorm-bench-policy
aws iam delete-role --profile admin --role-name bamstorm-bench-ec2-role
aws ecr delete-repository --profile admin --region us-east-1 \
  --repository-name bamstorm-bench --force
```
