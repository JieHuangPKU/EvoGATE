#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
REMOTE_REPO="jiehuang@162.105.248.143:/home/jiehuang/software/fungi/ProGATE_v2"
DRY_RUN=0

usage() {
  cat <<'EOF'
Usage:
  scripts/sync_phase1_bingo_refactor_to_server.sh [--dry-run]

Description:
  Sync the minimal ESM2 cache preparation workflow files and required processed
  inputs from the local ProGATE_v2 checkout to the remote ProGATE_v2 checkout.

Options:
  --dry-run   Preview rsync actions without modifying the remote side.
  -h, --help  Show this help message.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

FILES=(
  "configs/prepare_esm2_cache.yaml"
  "configs/phase1_esm2_mlp.yaml"
  "workflow/prepare_esm2_cache.smk"
  "workflow/phase1_esm2_mlp.smk"
  "scripts/run_prepare_esm2_cache.sh"
  "scripts/run_phase1_esm2_mlp.sh"
  "scripts/sync_phase1_bingo_refactor_to_server.sh"
  "src/features/extract_esm2_pooled.py"
  "data/processed/essential_gene/human/protein.fa"
  "data/processed/essential_gene/human/labels.standard.tsv"
  "data/processed/essential_gene/celegans/protein.fa"
  "data/processed/essential_gene/celegans/labels.standard.tsv"
  "data/processed/essential_gene/scerevisiae/protein.fa"
  "data/processed/essential_gene/scerevisiae/labels.standard.tsv"
  "data/processed/essential_gene/melanogaster/protein.fa"
  "data/processed/essential_gene/melanogaster/labels.standard.tsv"
  "data/processed/essential_gene/fgraminearum/protein.fa"
  "data/processed/essential_gene/fgraminearum/oldlabel/labels.tsv"
  "data/processed/essential_gene/fgraminearum/newlabel/labels.tsv"
)

RSYNC_ARGS=(-avzL --progress --relative)
if [[ "${DRY_RUN}" -eq 1 ]]; then
  RSYNC_ARGS+=(-n)
fi

cd "${REPO_ROOT}"

for relpath in "${FILES[@]}"; do
  if [[ ! -e "${relpath}" ]]; then
    echo "Required path does not exist: ${relpath}" >&2
    exit 1
  fi
done

echo "[sync_phase1_bingo_refactor_to_server] repo_root=${REPO_ROOT}"
echo "[sync_phase1_bingo_refactor_to_server] remote_repo=${REMOTE_REPO}"
echo "[sync_phase1_bingo_refactor_to_server] dry_run=${DRY_RUN}"
echo "[sync_phase1_bingo_refactor_to_server] file_count=${#FILES[@]}"

for relpath in "${FILES[@]}"; do
  echo "  - ${relpath}"
done

for relpath in "${FILES[@]}"; do
  rsync "${RSYNC_ARGS[@]}" "./${relpath}" "${REMOTE_REPO}/"
done
