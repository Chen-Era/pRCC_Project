#!/bin/bash
set -euo pipefail

# Usage: ./SCENIC.sh <expression_csv> <tf_list.txt> <output_dir> [workers=8] [seed=2022]
if [ $# -lt 3 ]; then
  echo "Usage: $0 <expression_csv> <tf_list.txt> <output_dir> [workers=8] [seed=2022]"
  exit 1
fi

matrix="$1"
tf="$2"
outdir="$3"
workers="${4:-8}"
seed="${5:-2022}"
method="grnboost2"

mkdir -p "$outdir"

echo "Running GRN inference..."
pyscenic grn -o "${outdir}/grn_output.tsv" -t -m "${method}" --num_workers "${workers}" "${matrix}" "${tf}" --seed "${seed}"
echo "GRN step finished"

echo "Running motif enrichment (ctx)..."
pyscenic ctx "${outdir}/grn_output.tsv" "10kb_up_and_down_tss.feather" "500bp_up_and_100bp_down_tss.feather" -t -o "${outdir}/ctx_output.tsv" --num_workers "${workers}" --nes_threshold 3 --auc_threshold 0.05 --rank_threshold 5000 --annotations_fname "Anno.txt" --expression_mtx_fname "${matrix}"
echo "CTX step finished"

echo "Running AUCell..."
pyscenic aucell "${matrix}" "${outdir}/ctx_output.tsv" -t -o "${outdir}/aucell_output.tsv" --seed "${seed}" --num_workers "${workers}"
echo "AUCell step finished. Results written to ${outdir}"