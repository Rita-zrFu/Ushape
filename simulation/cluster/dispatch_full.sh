#!/bin/bash
# dispatch_full.sh — full USR simulation dispatch (revised 2026-03-31)
# 23 scenarios × 1050 replicates, DEoptim bootstrap, 10h wall
# Design: beta1=2 main, min EV replaces Gumbel, G_exp c=1,
#         t0_CR: logistic→2.0, exp→0.3, t0_S: exp→1.0
#
# Run from the simulation/ directory of the repo:
#   bash cluster/dispatch_full.sh
# or from any cwd:
#   bash /path/to/repo/simulation/cluster/dispatch_full.sh
#
# Args passed to simulation_one_rep.R:
#   scenario_id seed n ei G_type censor_lo censor_hi run_competitors [beta1]
#
# Environment overrides:
#   N_REPS       replicates per scenario (default 1050; keeps ~1000 after loss)
#   WALL         SLURM wall time (default 10:00:00)
#   MEM          SLURM memory (default 4G)
#   RESULTS_DIR  per-rep RDS output root (default $SIM_ROOT/results)
#   MYCINDEX_CPP_PATH  path to C++ kernel (default $SIM_ROOT/myCindex.cpp)

# Locate simulation/ root (parent of cluster/) so this script runs from any cwd.
SIM_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "${SIM_ROOT}"

N_REPS=${N_REPS:-1050}
SCRIPT=${SIM_ROOT}/simulation_one_rep.R
CPP_PATH=${MYCINDEX_CPP_PATH:-${SIM_ROOT}/myCindex.cpp}
WALL=${WALL:-10:00:00}
MEM=${MEM:-4G}
RESULTS_DIR=${RESULTS_DIR:-${SIM_ROOT}/results}
LOGDIR=${SIM_ROOT}/cluster/logs

mkdir -p "${LOGDIR}" "${RESULTS_DIR}"

submit() {
  local SID=$1 N=$2 EI=$3 GT=$4 CLO=$5 CHI=$6 B1=${7:-2}
  mkdir -p "${RESULTS_DIR}/${SID}"
  JOB=$(sbatch --parsable \
    --array=1-${N_REPS} \
    --job-name=usr_${SID} \
    --cpus-per-task=1 --mem=${MEM} --time=${WALL} \
    --output="${LOGDIR}/${SID}_%A_%a.out" \
    --error="${LOGDIR}/${SID}_%A_%a.err" \
    --wrap="source /etc/profile && module load conda_R/4.4 && \
      cd ${SIM_ROOT} && \
      export MYCINDEX_CPP_PATH=${CPP_PATH} && \
      export RESULTS_DIR=${RESULTS_DIR} && \
      Rscript ${SCRIPT} ${SID} \$SLURM_ARRAY_TASK_ID ${N} ${EI} ${GT} ${CLO} ${CHI} TRUE ${B1}")
  echo "${SID} (n=${N}, ${EI}, ${GT}, b1=${B1}): Job ${JOB}"
}

echo "=============================================="
echo "USR Simulation Full Dispatch"
echo "SIM_ROOT:    ${SIM_ROOT}"
echo "RESULTS_DIR: ${RESULTS_DIR}"
echo "Reps:        ${N_REPS} | Wall: ${WALL} | Mem: ${MEM}"
echo "Date:        $(date)"
echo "=============================================="

echo ""
echo "=== Main scenarios (beta1=2, Xc=1.50) ==="

echo "--- Logistic / Normal ---"
submit S01 200  norm logistic 1.9800 9.9800 2
submit S02 200  norm logistic 0.3000 8.3000 2
submit S03 200  norm logistic 0.1100 4.9100 2
submit S04 500  norm logistic 1.9800 9.9800 2
submit S05 500  norm logistic 0.3000 8.3000 2
submit S06 500  norm logistic 0.1100 4.9100 2
submit S07 1000 norm logistic 0.3000 8.3000 2

echo "--- Logistic / Min EV ---"
submit S08 200  ev logistic 0.3300 8.3300 2
submit S09 500  ev logistic 0.3300 8.3300 2
submit S10 1000 ev logistic 0.3300 8.3300 2

echo "--- Exp(c=1) / Normal ---"
submit S11 200  norm exp 0.0400 4.0400 2
submit S12 200  norm exp 0.0100 1.9100 2
submit S13 200  norm exp 0.0300 0.8300 2
submit S14 500  norm exp 0.0400 4.0400 2
submit S15 500  norm exp 0.0100 1.9100 2
submit S16 500  norm exp 0.0300 0.8300 2
submit S17 1000 norm exp 0.0100 1.9100 2

echo "--- Exp(c=1) / Min EV ---"
submit S18 200  ev exp 0.0400 1.8400 2
submit S19 500  ev exp 0.0400 1.8400 2
submit S20 1000 ev exp 0.0400 1.8400 2

echo ""
echo "=== Sensitivity scenarios (beta1=1, Xc=2.25) ==="

echo "--- Logistic / Min EV ---"
submit S21 200  ev logistic 1.9000 7.9000 1
submit S22 500  ev logistic 1.9000 7.9000 1
submit S23 1000 ev logistic 1.9000 7.9000 1

echo "=============================================="
echo "Total: 23 scenarios × ${N_REPS} reps"
echo "Monitor:   sacct -u \$USER --name=usr_* --format=JobID,State,Elapsed,ExitCode"
echo "Aggregate: Rscript ${SIM_ROOT}/summarize_all_to_csv.R ${RESULTS_DIR} ${SIM_ROOT}/csv_output"
echo "=============================================="
