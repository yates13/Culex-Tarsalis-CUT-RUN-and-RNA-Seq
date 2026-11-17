#!/bin/bash
# config.sh - single source of truth for paths, environments and Slurm settings
# for the Culex tarsalis CUT&RUN + RNA-seq pipeline on Alpine.
#
# Use: add ONE line near the top of every script, after `set -euo pipefail`:
#   source /scratch/alpine/.colostate.edu/c832500103/Cxt_Cut_Run_Pipeline/Culex-Tarsalis-CUT-RUN-and-RNA-Seq/config.sh
# Use the absolute path. Slurm copies the job script elsewhere before running it,
# so a path relative to the script will not resolve.
#
# Rule: no script types a path that is defined here. Change it here, once.
# Lines marked VERIFY came from earlier notes and should be checked with `ls`.

# ---- Roots (the forms WITHOUT "@", which breaks Perl-based tools) ----
export SCRATCH=/scratch/alpine/.colostate.edu/c832500103
export PROJECTS=/projects/.colostate.edu/c832500103

# ---- Conda (in /projects: backed up, not purged after 90 days) ----
export CONDA_BASE="$PROJECTS/miniforge3"        # VERIFY
export CONDA_ENVS="$PROJECTS/conda_envs"        # VERIFY

# One env per step. Activate by full path: load_env "$ENV_FASTP"
export ENV_FASTP="$CONDA_ENVS/fastp_env"
export ENV_MULTIQC="$CONDA_ENVS/rnaPseudo_clean2"
export ENV_HISAT2="$CONDA_ENVS/hisat2_env"
export ENV_MACS2="$CONDA_ENVS/macs2_env"
export ENV_DEEPTOOLS="$CONDA_ENVS/deeptools_kernel_v2"
export ENV_HOMER="$CONDA_ENVS/homer_env"

# ---- Reference: assembly CtarK1 ----
export GENOME_FA="$SCRATCH/RNA_Seq/Genome/CtarK1.fa"
export GTF="$SCRATCH/RNA_Seq/Genome/GTF/Culex-tarsalis_knwr_BASEFEATURES_CtarK1.1_filtered_WNV_rep_fixed.gtf"

# ---- Pipeline layout ----
export PIPE="$SCRATCH/Cxt_Cut_Run_Pipeline"
export REPO="$PIPE/Culex-Tarsalis-CUT-RUN-and-RNA-Seq"   # VERIFY: location of the git clone
export RAW_DIR="$PIPE/Cxt_Cut_Run_Raw_Files"             # VERIFY: parent directory
export GENOME_DIR="$PIPE/00_Cxt_Genome"
export HISAT2_INDEX_DIR="$PIPE/00_Cxt_Indexed_Genome"
export ROUND="$PIPE/Round1_Round2_With_Duplicates"

export FASTP_DIR="$ROUND/01_Cxt_FastP"
export FASTP_OUT="$FASTP_DIR/01.1_Cxt_fastp_out"
export FASTP_MULTIQC="$FASTP_DIR/01.2_FastP_MultiQC/MultiQC_out"
export ALIGN_DIR="$ROUND/02_Cxt_Alignment"
export PEAKS_DIR="$ROUND/03_MACS2_Peaks_New_Peaks"
export DEEPTOOLS_DIR="$ROUND/04_Deeptools"
export HOMER_DIR="$ROUND/05_HomerAnnotation"

# ---- Slurm settings (for reference; #SBATCH lines cannot read variables) ----
# --partition=acpu  --qos=cpu-normal  --mail-type=END,FAIL
# --output=logs/%x.%A_%a.out  --error=logs/%x.%A_%a.err   (mkdir -p logs first)

# ---- Helpers ----

# load_env <env path>: activate a conda env by full path.
load_env() {
    set +u    # conda activation scripts can trip `set -u`
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate "$1"
    set -u
}

# need <path>...: stop with a clear message if any input is missing or empty.
need() {
    local p
    for p in "$@"; do
        if [[ ! -s "$p" ]]; then
            echo "ERROR: missing or empty: $p" >&2
            exit 1
        fi
    done
}

# expect_lines <file> <n>: stop if a file does not have exactly n lines.
# Use after joins and filters so a silent mismatch cannot pass.
expect_lines() {
    local n
    n=$(wc -l < "$1")
    if [[ "$n" -ne "$2" ]]; then
        echo "ERROR: $1 has $n lines, expected $2" >&2
        exit 1
    fi
}

# check_config: run once by hand (`source config.sh; check_config`) to see
# which configured paths exist. Reports; does not exit.
check_config() {
    local v
    for v in SCRATCH PROJECTS CONDA_BASE CONDA_ENVS ENV_FASTP ENV_MULTIQC \
             ENV_HISAT2 ENV_MACS2 ENV_DEEPTOOLS ENV_HOMER GENOME_FA GTF PIPE \
             REPO RAW_DIR GENOME_DIR HISAT2_INDEX_DIR ROUND FASTP_DIR FASTP_OUT \
             ALIGN_DIR PEAKS_DIR DEEPTOOLS_DIR HOMER_DIR; do
        if [[ -e "${!v}" ]]; then
            echo "ok       $v = ${!v}"
        else
            echo "MISSING  $v = ${!v}"
        fi
    done
}
