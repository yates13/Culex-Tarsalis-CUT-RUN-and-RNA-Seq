# Project context: Culex tarsalis CUT&RUN + RNA-seq pipeline

Read this before writing or changing any script in this repo. It holds the
facts an assistant would otherwise have to ask for. Paths live in `config.sh`;
solved bugs live in `GOTCHAS.md`; the current state lives in `HANDOFF.md`.

Items marked **VERIFY** came from earlier notes and may be out of date.

## What this is

CUT&RUN chromatin profiling of *Culex tarsalis* mosquitoes, aligned to genome
assembly **CtarK1**, with matching RNA-seq. The RNA-seq scripts live in their
own repo, `Culex-RNA-Seq`; conda env specs live in `Alpine_Conda_Envs`.
One near-term goal is reproducing
publication-style figures (peaks by genomic feature; peak distance to nearest
gene start) from this data.

## Experimental design

| Factor | Values |
| --- | --- |
| Target | `AC` = H3K27ac (narrow peaks), `ME` = H3K9me3 (broad peaks), `SR` = SREBP (narrow peaks), `IN` = input control |
| Group | `BF`, `MP`, `SF`, `ZH` |
| Day | `D1`, `D3` |
| Replicates | up to 5 per condition |

Sample names: `{group}-{target}-{day}-R{rep}`, for example `BF-AC-D1-R2`.

## Where it runs

- **Cluster:** Alpine, Slurm. Login nodes are for editing and submitting only.
- **Slurm:** `--partition=acpu`, `--qos=cpu-normal`, `--mail-type=END,FAIL`.
- **Logs:** every script writes Slurm logs to `logs/`.
- **Paths:** always the form without `@`:
  - scratch: `/scratch/alpine/.colostate.edu/c832500103`
  - projects: `/projects/.colostate.edu/c832500103`
- **Storage rules:** scratch files are removed 90 days after creation;
  `/projects` is backed up with a 250 GB quota. Conda lives in `/projects`.
  Results worth keeping are backed up off scratch.
- **Conda:** `miniforge3` and `conda_envs` under the projects path. Activate
  by full path with `load_env` from `config.sh`. Env specs (`envs/*.yml`) and
  the build and verify scripts are in the `Alpine_Conda_Envs` repo.

## Pipeline steps

All under `Round1_Round2_With_Duplicates/`.

| Step | Scripts | Env | Notes |
| --- | --- | --- | --- |
| `01_Cxt_FastP` | 01.1 fastp array (96 tasks); rename; MultiQC | `fastp_env`, `rnaPseudo_clean2` | Output in `01.1_Cxt_fastp_out`. Trimmed names are fixed from `_Cxt` to `-Cxt`. MultiQC goes to `01.2_FastP_MultiQC/MultiQC_out`. |
| `02_Cxt_Alignment` | `02.1_runBuildCtark1Genome_alpine.sh`, `02.2_Cxt_AlignArray_Hisat.sh`, `02.2.5` `MultiQC.sh`, `02.3_runBam_array.sh`, `02.4_BamMerge_With_replicate_counts_array.sh` | `hisat2_env`, `rnaPseudo_clean2` | HISAT2 against the CtarK1 index. |
| `03_MACS2_Peaks_New_Peaks` | `03.1_New_flags_Macs2_PeakCalling_with_replicate_count.sh` | `macs2_env` | One job looping over merged BAMs. `BAMPE`; `--broad` for `ME`. Control is `<GROUP>-IN-D1*`. |
| `04_Deeptools` | `00_Create_Bigwigs_Individual_Bams_ARRAY.sh` (array 1-72), `04.1_ComputeMatrix_Individual_array.sh` (array 0-3), `04.1.2_..._RUN_THIS_AFTER.sh` | `deeptools_kernel_v2` | bamCompare input-subtracted bigWigs; computeMatrix per group; then `cbind` of the four group matrices. |
| `05_HomerAnnotation` | HOMER annotation, then R/ggplot2 figures | `homer_env` | |

## Known temporary choices

- **Input pairing for bigWigs is temporary.** Every day uses a D1 input of the
  same group. Switch to replicate-matched inputs once all input replicates exist.
- **`--skipZeros` in 04.1 can make groups differ in region count,** and `cbind`
  matches regions by name against the first matrix. Compare `wc -l region_*.bed`
  across groups; if they differ, rerun 04.1 without `--skipZeros`.

## How I want help delivered

- Recommend only changes that are necessary. Do not restyle working scripts.
- When only a few lines change, say which lines, instead of rewriting the script.
- Single-line commands for the terminal, with no backslash continuations.
- Scripts as copyable code blocks.
- Numbered scripts (`01.1`, `02.3`) in run order; keep existing names.
- Verify each step's output before moving to the next. Tell me what to check.
- Be direct and brief.

## Rules for new or edited scripts

1. `set -euo pipefail`, then `source` the absolute path to `config.sh`.
2. No path typed in a script if `config.sh` defines it.
3. Check inputs with `need` before running; check counts with `expect_lines`
   after joins and filters.
4. Array jobs pick their sample with `sed -n "${SLURM_ARRAY_TASK_ID}p" list.txt`.
5. Write to a temp name and rename on success.
6. Test tools on a compute node with `srun`, not on the login node.
