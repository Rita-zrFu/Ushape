#!/bin/bash -l
# dispatch.sh
# Submit the full USR simulation as a set of SLURM array jobs.
#
# Reads scenarios from ../SCENARIOS.csv (single source of truth for scenario
# definitions) and submits one array job per scenario row.
#
# Usage:
#   bash dispatch.sh                  # submit all 23 scenarios at default reps
#   N_REPS=200 bash dispatch.sh       # use a smaller rep count for a pilot
#
# Environment overrides (all have sensible defaults):
#   N_REPS              number of replicates per scenario (default 1100; target 1000 with ~10% buffer)
#   WALL                SLURM wall time (default 4:00:00)
#   MEM                 SLURM memory request (default 4G)
#   RESULTS_DIR         directory for per-rep RDS files (default $REPRO_ROOT/results)
#   SCENARIOS_CSV       scenario config (default $REPRO_ROOT/SCENARIOS.csv)
#   MYCINDEX_CPP_PATH   path to C++ kernel (default $REPRO_ROOT/code/myCindex.cpp)
#   RUN_COMPETITORS     "TRUE" or "FALSE" (default TRUE)
#
# Cluster expectations:
#   - SLURM with QOS allowing array jobs
#   - module load conda_R/4.4 (or equivalent providing R >= 4.4 + Rcpp toolchain)
#   - C++11 compiler (gcc 11+) for Rcpp::sourceCpp on first run

# Locate repro root (parent of cluster/) so the script works from any cwd
REPRO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "${REPRO_ROOT}"

N_REPS=${N_REPS:-1100}
WALL=${WALL:-4:00:00}
MEM=${MEM:-4G}
RESULTS_DIR=${RESULTS_DIR:-${REPRO_ROOT}/results}
SCENARIOS_CSV=${SCENARIOS_CSV:-${REPRO_ROOT}/SCENARIOS.csv}
MYCINDEX_CPP_PATH=${MYCINDEX_CPP_PATH:-${REPRO_ROOT}/code/myCindex.cpp}
RUN_COMPETITORS=${RUN_COMPETITORS:-TRUE}
SCRIPT=${REPRO_ROOT}/code/simulation_one_rep.R
LOGDIR=${REPRO_ROOT}/cluster/logs

if [ ! -f "${SCENARIOS_CSV}" ]; then
  echo "ERROR: SCENARIOS_CSV not found at ${SCENARIOS_CSV}"
  exit 1
fi
if [ ! -f "${SCRIPT}" ]; then
  echo "ERROR: simulation script not found at ${SCRIPT}"
  exit 1
fi
if [ ! -f "${MYCINDEX_CPP_PATH}" ]; then
  echo "ERROR: C++ kernel not found at ${MYCINDEX_CPP_PATH}"
  exit 1
fi

mkdir -p "${LOGDIR}" "${RESULTS_DIR}"

submit_one() {
  local SID=$1 N=$2 EI=$3 GT=$4 CLO=$5 CHI=$6 B1=$7
  mkdir -p "${RESULTS_DIR}/${SID}"
  JOB=$(sbatch --parsable \
    --array=1-${N_REPS} \
    --job-name=usr_${SID} \
    --cpus-per-task=1 --mem=${MEM} --time=${WALL} \
    --output="${LOGDIR}/${SID}_%A_%a.out" \
    --error="${LOGDIR}/${SID}_%A_%a.err" \
    --wrap="source /etc/profile && module load conda_R/4.4 && \
      cd ${REPRO_ROOT} && \
      export MYCINDEX_CPP_PATH=${MYCINDEX_CPP_PATH} && \
      export RESULTS_DIR=${RESULTS_DIR} && \
      Rscript ${SCRIPT} ${SID} \$SLURM_ARRAY_TASK_ID ${N} ${EI} ${GT} ${CLO} ${CHI} ${RUN_COMPETITORS} ${B1}")
  echo "${SID} (n=${N}, ei=${EI}, G=${GT}, b1=${B1}): job ${JOB}"
}

echo "=============================================="
echo "USR Simulation Dispatch"
echo "Date:         $(date)"
echo "REPRO_ROOT:   ${REPRO_ROOT}"
echo "Scenarios:    ${SCENARIOS_CSV}"
echo "Reps each:    ${N_REPS}"
echo "Wall time:    ${WALL}"
echo "Results dir:  ${RESULTS_DIR}"
echo "Competitors:  ${RUN_COMPETITORS}"
echo "=============================================="

# Parse SCENARIOS.csv. Expected header columns include:
#   scenario_id, n, ei, G_type, censor_lo, censor_hi, beta1
# (other columns are ignored at dispatch time)
header_line=$(head -n 1 "${SCENARIOS_CSV}")
IFS=',' read -ra cols <<< "${header_line}"

# Find column indices for the fields we need
col_idx() {
  local target=$1
  for i in "${!cols[@]}"; do
    if [ "${cols[$i]}" = "${target}" ]; then echo "$i"; return; fi
  done
  echo -1
}
I_SID=$(col_idx scenario_id)
I_N=$(col_idx n)
I_EI=$(col_idx ei)
I_G=$(col_idx G_type)
I_CLO=$(col_idx censor_lo)
I_CHI=$(col_idx censor_hi)
I_B1=$(col_idx beta1)

if [ "${I_SID}" -lt 0 ] || [ "${I_N}" -lt 0 ] || [ "${I_EI}" -lt 0 ] || \
   [ "${I_G}" -lt 0 ] || [ "${I_CLO}" -lt 0 ] || [ "${I_CHI}" -lt 0 ] || [ "${I_B1}" -lt 0 ]; then
  echo "ERROR: SCENARIOS.csv header missing required columns."
  echo "Required: scenario_id, n, ei, G_type, censor_lo, censor_hi, beta1"
  echo "Found:    ${header_line}"
  exit 1
fi

N_SUBMITTED=0
N_FAILED=0
tail -n +2 "${SCENARIOS_CSV}" | while IFS=',' read -ra row; do
  SID=${row[$I_SID]}
  N=${row[$I_N]}
  EI=${row[$I_EI]}
  GT=${row[$I_G]}
  CLO=${row[$I_CLO]}
  CHI=${row[$I_CHI]}
  B1=${row[$I_B1]}
  [ -z "${SID}" ] && continue
  if submit_one "${SID}" "${N}" "${EI}" "${GT}" "${CLO}" "${CHI}" "${B1}"; then
    N_SUBMITTED=$((N_SUBMITTED + 1))
  else
    N_FAILED=$((N_FAILED + 1))
  fi
done

echo "=============================================="
echo "Total submission requests issued."
echo "Monitor with: squeue -u \$USER | grep usr_"
echo "Aggregate after completion:"
echo "  Rscript ${REPRO_ROOT}/code/summarize_all_to_csv.R ${RESULTS_DIR} ${REPRO_ROOT}/csv_output ${SCENARIOS_CSV}"
echo "=============================================="
