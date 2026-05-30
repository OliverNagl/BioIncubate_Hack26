#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH -J t3_sweep_driver
#SBATCH -o t3_capsids/logs/sweep_driver_%j.out
#SBATCH -e t3_capsids/logs/sweep_driver_%j.err
#
# --- hack26 overrides ---
# (later #SBATCH lines override earlier ones for the same directive)
#SBATCH --job-name=hack26_fold
#SBATCH --output=logs/fold_%j.out
#SBATCH --error=logs/fold_%j.err
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
# JUSTIFIED DEVIATIONS (the copied verbatim block + the spec's overrides don't
# run as-is on this cluster; the spec permits deviations with justification):
#  1. `-p cpu` can't satisfy a GPU request, but ESMFold2 runs
#     `apptainer exec --nv`. Move to a GPU partition. (AF3 is dispatched by the
#     pipelines engine to gpu-train as its own child SLURM job.)
#  2. This cluster's `gpu` partition rejects an untyped `--gres=gpu:1`
#     ("did not specify a GPU type"); GRES is `gpu:small`/`gpu:large`. Use
#     gpu:small:1 (enough VRAM for a ~400-residue ESMFold2 co-fold).
#  3. The copied `--mem=8G` OOM-kills ESMFold2 (SLURM OUT_OF_MEMORY at ~98s,
#     job 2314798) — its weights + the container need far more host RAM. 64G.
#SBATCH -p gpu
#SBATCH --gres=gpu:small:1
#SBATCH --mem=64G

set -euo pipefail

# run_sweep.sbatch has no module-load / conda / venv activation block to copy
# verbatim -- it relies on sweep_pipeline.py's shebang to self-activate the
# pipelines .venv. Our pipeline body calls `python ...` directly, so we
# activate that same venv here. Resolve the repo root from the submit dir
# (this script lives in <repo>/t3_capsids/data) or honour $PIPELINES_REPO.
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$SUBMIT_DIR"
REPO_ROOT="${PIPELINES_REPO:-$(cd "$SUBMIT_DIR/../.." && pwd)}"
# shellcheck disable=SC1091
source "$REPO_ROOT/.venv/bin/activate"

mkdir -p logs

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-local} dir=$SUBMIT_DIR"
echo "[$(date)] Preparing inputs"
python prepare_inputs.py

echo "[$(date)] Running fold pipeline"
python fold_pipeline.py \
    --variants all \
    --engines  all \
    ${EXTRA_ARGS:-}

echo "[$(date)] Done — results in results/"

# Usage:
#   sbatch run_fold.sbatch
#   EXTRA_ARGS="--variants Cetuximab --engines esm2" sbatch run_fold.sbatch
#   EXTRA_ARGS="--dry-run" sbatch run_fold.sbatch
#   EXTRA_ARGS="--seeds 42,123" sbatch run_fold.sbatch
