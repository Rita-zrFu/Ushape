#!/bin/bash
# dispatch_reprocess_oracle_v2.sh — oracle-only reprocessing (2026-06-25)
#
# Reprocesses existing per-rep .rds files (produced by dispatch_full.sh) to add
# a new `oracle_v2` field computed with the corrected Z-specific 7-parameter
# (B_R) oracle Cox spec. The original `oracle` field is preserved for audit.
#
# Does NOT touch MCE, bootstrap, or any non-oracle competitor.
#
# Run after dispatch_full.sh completes:
#   bash cluster/dispatch_reprocess_oracle_v2.sh
#
# Environment overrides:
#   RESULTS_DIR  per-rep RDS root, updated in place (default $SIM_ROOT/results)
#   WALL         SLURM wall time (default 2:00:00)
#   MEM          SLURM memory (default 4G)

SIM_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "${SIM_ROOT}"

SCRIPT=${SIM_ROOT}/reprocess_oracle_v2.R
RESULTS_DIR=${RESULTS_DIR:-${SIM_ROOT}/results}
WALL=${WALL:-2:00:00}
MEM=${MEM:-4G}
LOGDIR=${SIM_ROOT}/cluster/logs

mkdir -p "${LOGDIR}"

submit() {
  local SID=$1
  if [ ! -d "${RESULTS_DIR}/${SID}" ]; then
    echo "  ${SID}: SKIP — no results directory ${RESULTS_DIR}/${SID}"
    return
  fi
  local NRDS=$(ls ${RESULTS_DIR}/${SID}/*.rds 2>/dev/null | wc -l)
  if [ "${NRDS}" -eq 0 ]; then
    echo "  ${SID}: SKIP — no .rds files in ${RESULTS_DIR}/${SID}"
    return
  fi
  JOB=$(sbatch --parsable \
    --job-name=ora_v2_${SID} \
    --cpus-per-task=1 --mem=${MEM} --time=${WALL} \
    --output="${LOGDIR}/oracle_v2_${SID}_%j.out" \
    --error="${LOGDIR}/oracle_v2_${SID}_%j.err" \
    --wrap="source /etc/profile && module load conda_R/4.4 && \
      cd ${SIM_ROOT} && \
      Rscript ${SCRIPT} ${SID} ${RESULTS_DIR}")
  echo "  ${SID} (${NRDS} reps): Job ${JOB}"
}

echo "=============================================="
echo "Oracle v2 reprocessing dispatch"
echo "SIM_ROOT:    ${SIM_ROOT}"
echo "RESULTS_DIR: ${RESULTS_DIR}"
echo "Wall / mem:  ${WALL} / ${MEM}"
echo "Date:        $(date)"
echo "=============================================="
echo ""
echo "Submitting one job per scenario (sequential reprocessing within scenario)..."

for sid in S01 S02 S03 S04 S05 S06 S07 \
           S08 S09 S10 \
           S11 S12 S13 S14 S15 S16 S17 \
           S18 S19 S20 \
           S21 S22 S23; do
  submit ${sid}
done

echo ""
echo "=============================================="
echo "Monitor:   sacct -u \$USER --name=ora_v2_* --format=JobID,JobName,State,Elapsed,ExitCode"
echo "Logs:      ${LOGDIR}/oracle_v2_S*.{out,err}"
echo "Aggregate: Rscript ${SIM_ROOT}/summarize_all_to_csv.R ${RESULTS_DIR} ${SIM_ROOT}/csv_output"
echo "=============================================="
