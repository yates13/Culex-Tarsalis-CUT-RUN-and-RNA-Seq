# Culex tarsalis CUT&RUN and RNA-seq Analysis

This repository holds the analysis pipeline for a study of chromatin regulation in *Culex tarsalis* during infection with Rift Valley Fever Virus (RVFV). We used **CUT&RUN** and **RNA-seq** data to compare infection with two RVFV strains, **ZH501** and **MP12**.

<!-- TODO (Study overview): confirm H3K9me vs H3K9me3, explain the BF and SF sites, confirm strain/site wording -->
## Study overview

- **Virus strains:** RVFV ZH501 and MP12
- **CUT&RUN targets:**
  - H3K27ac (active regulatory regions)
  - H3K9me (repressive chromatin)
  - SREBP (transcription factor)
- **Additional data:** RNA-seq, for comparing chromatin changes with gene expression
- **Timepoints:** Day 1 and Day 3
- **Compute:** CU/CSU Alpine HPC, jobs run with SLURM

### Sample names

Samples are named `<TREATMENT>-<TARGET>-<DAY>-R<n>-Cxt`, for example `ZH-AC-D1-R1-Cxt`.

| Field | Values |
|-------|--------|
| `TREATMENT` | BF, MP, SF, ZH |
| `TARGET` | AC (H3K27ac), ME (H3K9me, broad domains), SR (SREBP), IN (input control) |
| `DAY` | D1, D3 |
| `R<n>` | replicate number, up to 5 |

After replicates are merged (Step 2.4) the name becomes `<TREATMENT>-<TARGET>-<DAY>-rep<digits>`. The digits list which replicates went in: `rep123` = R1, R2 and R3.

## Pipeline

Run the steps in order. Each folder has its own scripts, numbered in run order.

| Step | Folder | What it does | Status |
|------|--------|--------------|--------|
| 1 | [`CUT-RUN/01_Cxt_FastP`](CUT-RUN/01_Cxt_FastP) | Read trimming with fastp, sample renaming, MultiQC summary of the fastp reports | Scripts updated |
| 2 | [`CUT-RUN/02_Cxt_Alignment`](CUT-RUN/02_Cxt_Alignment) | HISAT2 index and alignment, alignment QC, BAM conversion, replicate merging | Scripts updated (`Gvis_Function` and `HOW_TO_MAKE_FASTQTRIMMED_FILE` still to do) |
| 3 | [`CUT-RUN/03_MACS2_Peaks_New_Peaks`](CUT-RUN/03_MACS2_Peaks_New_Peaks) | Peak calling with MACS2 | `03.1` updated |
| 4 | [`CUT-RUN/04_Deeptools`](CUT-RUN/04_Deeptools) | bigWigs, signal matrices, heatmaps and profiles with deepTools | `04.0`, `04.1`, `04.1.2`, `04.2` and `04.3` updated; the representative-replicate scripts still to do |
| 5 | [`CUT-RUN/05_HomerAnnotation`](CUT-RUN/05_HomerAnnotation) | Peak annotation with HOMER | Not yet updated |

To view tracks in IGV, see [`IGV_ChIP-Seq_Instructions.pdf`](CUT-RUN/IGV_ChIP-Seq_Instructions.pdf).

## Setup

### Conda environments

Each step activates its own conda environment.

| Environment | Used in |
|-------------|---------|
| `fastp_env` | Step 1.1, trimming |
| `rnaPseudo_clean` | Step 1.2 (MultiQC) and all of Step 2 (HISAT2, samtools, MultiQC) |
| `macs2_env` | Step 3, peak calling |
| `deeptools_kernel_v2` | Step 4, see `How_to_make_deeptools_kernel_env.md` |
| `homer_env` | Step 5, HOMER annotation |

<!-- TODO: add .yml exports for rnaPseudo_clean and homer_env to CondaEnvs/, and update the conda lines in the scripts if the environments move -->
The folder [`CondaEnvs`](CondaEnvs) holds one `.yml` file per environment. Build the ones that are missing with:

```bash
bash CondaEnvs/make_envs.sh <ENV_DIR>
```

where `<ENV_DIR>` is the folder that will hold the environments (for example a `conda_envs` folder in `/projects`). Environments that already exist are skipped. To save an existing environment so it can be rebuilt later:

```bash
conda env export --from-history -p <ENV_DIR>/<name> | grep -v '^prefix:' > CondaEnvs/<name>.yml
conda env export --no-builds    -p <ENV_DIR>/<name> | grep -v '^prefix:' > CondaEnvs/<name>.full.yml
```

Scripts load conda like this. Set `SCRATCH` to your own base path first:

```bash
SCRATCH=/scratch/alpine/<your_username>
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh
conda activate ${SCRATCH}/conda_envs/<env_name>
```

### Data retention (important)

Files in `/scratch/alpine` are removed 90 days after they are created, and scratch is not backed up. Conda environments and pipeline outputs on scratch can disappear. To protect your work:

- Keep scripts in this repository and environment files in `CondaEnvs/`.
- Build conda environments (and Miniconda) in `/projects`, which is backed up, not in scratch.
- Copy the results you cannot regenerate off scratch: peak calls, matrices, merged BAMs and bigWigs. Keep a second copy of the raw reads somewhere other than scratch.
- Intermediate files (trimmed FASTQs, SAMs, individual BAMs) can be regenerated from the raw reads and these scripts.

## Running

1. Edit the paths at the top of each script to match your directory, and the `--mail-user` address in the `#SBATCH` header.
2. All scripts use `--partition=acpu` and `--qos=cpu-normal`, and write their Slurm logs to a `logs/` folder. Slurm cannot create that folder, and a job whose log folder is missing fails without a log, so run `mkdir -p logs` first.
3. Submit each script from the folder that contains it (see the table below), then check the outputs before moving on.
4. See the `How_to_make_*` files in each folder for building sample and treatment lists.

### Where to run each script

With `P` set to the pipeline output folder (`.../Round1_Round2_With_Duplicates`):

| Scripts | Submit from |
|---------|-------------|
| `01.1_fastp_Cxt.sh`, `01.2_multiqc_fastp.sh` | `$P/01_Cxt_FastP` |
| `02.2`, `02.2.5` (`MultiQC.sh`), `02.3`, `02.4` | `$P/02_Cxt_Alignment` |
| `03.1_New_flags_Macs2_PeakCalling_with_replicate_count.sh` | `$P/03_MACS2_Peaks_New_Peaks` |
| `04.0`, `04.1`, `04.1.2`, `04.2`, `04.3` | `$P/04_Deeptools` (it needs `bw_chunks/` and writes `matrix_*.gz` and `region_*.bed` there) |

## Step details

### Step 1: Trimming and QC (`CUT-RUN/01_Cxt_FastP`)

<!-- TODO: explain "Round1_Round2_With_Duplicates" (the pipeline output folder name) if it should appear here -->
| Script | What it does | Environment |
|--------|--------------|-------------|
| `01.1_fastp_Cxt.sh` | Trims paired-end reads with fastp, one sample per SLURM array task. fastp also writes a QC report per sample. | `fastp_env` |
| `01.1_fix_sample_names.sh` | Renames `*_Cxt_R{1,2}_trimmed.fastq.gz` to `*-Cxt_R{1,2}_trimmed.fastq.gz` so sample names are consistent. | none |
| `01.2_multiqc_fastp.sh` | Combines the per-sample fastp reports into one MultiQC report. | `rnaPseudo_clean` |

**1. Trim reads.** Build the R1 list in the raw FASTQ folder. Its line count is the array size:

```bash
cd <raw_fastq_dir>
ls *R1_001.fastq.gz | sort > r1_file_list.txt
N=$(wc -l < r1_file_list.txt)
```

Edit `SCRATCH`, `INDIR`, `OUTDIR` and `--mail-user` at the top of `01.1_fastp_Cxt.sh`. Then, from the script folder:

```bash
mkdir -p logs
sbatch --array=1-${N} 01.1_fastp_Cxt.sh
```

A sample with no matching R2 file makes its task fail with `No R2 found`, so check `sacct` or the logs afterwards. Outputs per sample (in `OUTDIR`): `<sample>_R1_trimmed.fastq.gz`, `<sample>_R2_trimmed.fastq.gz`, `<sample>_fastp.html`, `<sample>_fastp.json`.

**2. Fix sample names.** Run once on the fastp output folder after all tasks finish:

```bash
bash 01.1_fix_sample_names.sh <path_to_01.1_Cxt_fastp_out>
```

Only the trimmed FASTQ files are renamed. The `.html` and `.json` reports keep their original names.

**3. Summarize QC.**

```bash
sbatch 01.2_multiqc_fastp.sh
```

The MultiQC report is written to `01.2_FastP_MultiQC/MultiQC_out/`.

### Step 2: Index, alignment and BAMs (`CUT-RUN/02_Cxt_Alignment`)

<!-- TODO: 02.1 uses a relative path (cd Cxt_Genome); switch it to absolute paths. Add Gvis_Function and HOW_TO_MAKE_FASTQTRIMMED_FILE. -->
| Script | What it does | Environment |
|--------|--------------|-------------|
| `02.1_runBuildCtark1Genome_apline.sh` | Builds the HISAT2 index for CtarK1 and writes `CtarK1_index.*.ht2` to `Cxt_Cut_Run_Pipeline/00_Cxt_Indexed_Genome/`. Run once. | `rnaPseudo_clean` |
| `02.2_Cxt_AlignArray_Hisat.sh` | Aligns each trimmed read pair with HISAT2, one sample per array task. Writes `<sample>.sam` and `<sample>.hisat2.log` to `02.2_Cxt_HisatAligned_Cleaned/`. | `rnaPseudo_clean` |
| `MultiQC.sh` (Step 2.2.5) | One MultiQC report covering the fastp and HISAT2 outputs. | `rnaPseudo_clean` |
| `02.3_runBam_array.sh` | Converts each SAM to a sorted, indexed BAM with samtools, one sample per array task. Safe to resubmit: finished samples are skipped. | `rnaPseudo_clean` |
| `02.4_BamMerge_With_replicate_counts_array.sh` | Merges the replicate BAMs of each condition, one condition per array task. | `rnaPseudo_clean` |

**2.1 Build the index (once).** Skip this if the eight `CtarK1_index.*.ht2` files already exist in `00_Cxt_Indexed_Genome/`.

```bash
mkdir -p logs
sbatch 02.1_runBuildCtark1Genome_apline.sh
```

**2.2 Align.** Array numbering here starts at 0, so the array is `0-(N-1)`:

```bash
mkdir -p logs
ls <fastp_out_dir> | grep '_R1_' | grep 'trimmed' > FastqTrimmed.txt
N=$(wc -l < FastqTrimmed.txt)
sbatch --array=0-$((N-1))%20 02.2_Cxt_AlignArray_Hisat.sh
```

HISAT2 runs with `--no-spliced-alignment` (CUT&RUN is genomic, not RNA), `--no-mixed` (only pairs where both mates align) and `--mm` (memory-mapped index). An array larger than the list makes the extra tasks quietly do nothing, so a successful job does not prove every sample ran. Check the results:

```bash
ls <OUT_DIR>/*.sam | wc -l                       # should equal N
sacct -j <jobid> --format=JobID,State,ExitCode   # look for FAILED or TIMEOUT
```

A task that hits the time limit leaves a truncated `.sam`. Rerun it with `sbatch --array=<task> 02.2_Cxt_AlignArray_Hisat.sh`.

**2.2.5 Alignment QC.** Run after Step 1 and Step 2.2 have finished for every sample:

```bash
sbatch MultiQC.sh
```

It scans the fastp `.json` reports and the `.hisat2.log` files, and writes `Cxt-CUTRUN-Pipeline-QC-fastp--HISAT2_multiqc_report.html` to `02.2.5_MultiQC_Hisat_aligned/`. Open the General Statistics table: the number of rows should match your sample count.

**2.3 Convert to BAM.** Once every `.sam` is confirmed good, build the list of SAM files (its line count is the array size, counting from 1):

```bash
mkdir -p logs
samdir=<path_to_02.2_Cxt_HisatAligned_Cleaned>
ls ${samdir}/*.sam > ${samdir}/sam_file_list.txt
N=$(wc -l < ${samdir}/sam_file_list.txt)
sbatch --array=1-${N} 02.3_runBam_array.sh
```

Check with `ls <outdir>/*.bam | wc -l` (should equal N) and `sacct`. A task fails if samtools hits a truncated SAM, so rerun Step 2.2 for that sample first. Once every BAM is confirmed good, the SAM files can be deleted to free scratch space.

**2.4 Merge replicates.** The array size is the number of conditions, not samples. Build the conditions list, then submit with that size:

```bash
bamdir=<path>/02.3_Cxt_BamConverted
mergedir=<path>/02.4_Cxt_BamMerged_with_replicate_count
mkdir -p logs $mergedir
ls ${bamdir}/*.bam | xargs -n1 basename | sed 's/-R[0-9]-Cxt\.bam//' | sort -u > ${mergedir}/conditions.txt
sbatch --array=1-$(wc -l < ${mergedir}/conditions.txt) 02.4_BamMerge_With_replicate_counts_array.sh
```

BAM names must follow `<condition>-R<n>-Cxt.bam`. Outputs per condition: `<condition>-rep<digits>-merged-sorted.bam` and its `.bai`. The unsorted `-merged.bam` intermediates can be deleted.

### Step 3: Peak calling (`CUT-RUN/03_MACS2_Peaks_New_Peaks`)

`03.1_New_flags_Macs2_PeakCalling_with_replicate_count.sh` calls peaks with MACS2 on every `*-merged-sorted.bam` from Step 2.4. It is a single job that loops over the samples, not an array. Environment: `macs2_env`.

```bash
mkdir -p logs
sbatch 03.1_New_flags_Macs2_PeakCalling_with_replicate_count.sh
```

- Input controls (names containing `-IN-`) are skipped as samples.
- ME (H3K9me) samples are called with `--broad`. AC and SR samples use MACS2's default narrow peak calling.
- BAMs are treated as paired-end (`-f BAMPE`), duplicates are kept (`--keep-dup all`), and the effective genome size is 800,000,000.
- Each sample is compared against the first `<TREATMENT>-IN-D1*-merged-sorted.bam`. All days (D1 and D3) therefore use a Day 1 input control.
- Output goes to `03.1_MacsPeaks_with_replicate_count/<sample>/`, with the MACS2 messages for all samples appended to `macs2.log`. Rerunning redoes every sample and overwrites its outputs.
- Samples with no input control, or no `rep<digits>` in the name, are reported and skipped. The job stops if MACS2 itself fails.

### Step 4: deepTools (`CUT-RUN/04_Deeptools`)

Environment for all Step 4 scripts: `deeptools_kernel_v2`. Submit from `04_Deeptools`.

**4.0 Input-subtracted bigWigs** (`00_Create_Bigwigs_Individual_Bams_ARRAY.sh`). One BAM per array task; each treatment BAM from Step 2.3 is compared with a Day 1 input BAM of the same treatment using `bamCompare` (subtract, BPM-normalized):

```bash
mkdir -p logs
ls <BAM_INPUT_DIR>/*-Cxt.bam | grep -v -- '-IN-' > treatment_bam_list.txt
sbatch --array=1-$(wc -l < treatment_bam_list.txt) 00_Create_Bigwigs_Individual_Bams_ARRAY.sh
```

Output: `<sample>_inputsubtracted.bw` in `00_BigWig_Files_Individual_Bams/`. Finished samples are skipped on resubmission.

**4.1 Signal matrices** (`04.1_ComputeMatrix_Individual_array.sh`). One group (BF, MP, SF, ZH) per array task, TSS -2 kb to +2 kb. It reads the bigWigs listed in `bw_chunks/chunk_<GROUP>.txt` and the `CtarK1_TSS.bed` regions, and writes `matrix_<GROUP>.gz` and `region_<GROUP>.bed`:

```bash
ls bw_chunks/
sbatch 04.1_ComputeMatrix_Individual_array.sh
```

**4.1.2 Combine the groups** (`04.1.2_ComputeMatrix_Individual_array_RUN_THIS_AFTER.sh`). Run after 4.1 has finished for every group. It uses `computeMatrixOperations cbind` to put the samples of all four groups side by side:

```bash
sbatch 04.1.2_ComputeMatrix_Individual_array_RUN_THIS_AFTER.sh
computeMatrixOperations info -m matrix_combined_gene_only_bed_file.gz
```

`cbind` matches regions by name against the first matrix. Regions missing from a later matrix become NA, and regions found only in a later matrix are dropped, without any message. Every group matrix should therefore be built from the same set of regions, so do not use `--skipZeros` in 4.1. Compare `wc -l region_*.bed` across the four groups to check.

**4.2 Profile plots by group** (`04.2_BF_...`, `04.2_MP_...`, `04.2_SF_...`, `04.2_ZH_Plot_profile_Individual_Bams_Dynamic_sample_naming.sh`). Run after 4.1.2. The four scripts are the same except for the group name, job name, plot title and line colors. Each is a single job (not an array) that draws every replicate of one group as its own line around the TSS (-2 kb to +2 kb), one plot per target (AC, ME, SR) and day (D1, D3):

```bash
mkdir -p logs
sbatch 04.2_BF_Plot_profile_Individual_Bams_Dynamic_sample_naming.sh
```

The samples are found from the bigWigs in `00_BigWig_Files_Individual_Bams/`, so new replicates are picked up automatically, and their names must match the sample labels in the combined matrix. A target/day with no samples for the group is skipped. Each plot can hold 8 replicates, the number of colors in `COLOR_PALETTE`. The y-axis is fixed at -0.3 to 0.8. Run the four group scripts one at a time, because they all write the same `matrix_<TARGET>_<DAY>_ALL.gz` files. Outputs: `Figures/<TARGET>_<DAY>_<GROUP>_profile.png` and `matrix_<TARGET>_<DAY>_<GROUP>.gz`. Check with `ls Figures/*_BF_profile.png` (6 plots) and `grep "samples:" logs/*.log`.

**4.3 Heatmaps** (`04.3_HeatMap_Individual_Bams_By_Treatment_updated.sh`). Run after 4.1.2. One heatmap per target and day, with every replicate of every group as its own panel. A title is required, and the script exits with an error without one. Each heatmap is titled `<Base Title> - <TARGET> <DAY>`:

```bash
mkdir -p logs
sbatch 04.3_HeatMap_Individual_Bams_By_Treatment_updated.sh "Base Title"
```

The sample list s read from the combined matrix, so no other input files are needed. Outputs: `Figures/matrix_<TARGET>_<DAY>_all_pops_heatmap.png` and `matrix_<TARGET>_<DAY>_all_pops.gz`. The log ends with a list of the heatmaps and their sizes.

<!-- TODO: document 04.2_Representative_Replicates and 04.3_Representative_Replicates_heat_maps once they are reviewed -->

### Step 5: HOMER annotation (`CUT-RUN/05_HomerAnnotation`)

<!-- TODO: scripts not yet reviewed -->
Not yet documented.

## Known issues and notes

- **Input pairing is temporary.** Steps 3 and 4.0 use a Day 1 input control for every sample, including Day 3, because not every treatment replicate has a same-numbered input yet. Switch to replicate-matched inputs once all input replicates exist.
- **Timed-out alignments.** A HISAT2 task that reaches its time limit leaves a truncated `.sam`; rerun that task before converting to BAM.
- **Email address.** Every script has a `--mail-user` line. Change it to your own.

## How to use Slurm

Common commands: sbatch, srun, scancel, sacct, squeue, sinfo

```bash
sbatch [options] script [args]                   # submit a job
scancel 12345                                    # cancel job 12345
scancel {31415..31425}                           # cancel a range of sequential job IDs
squeue | grep <your_username>                    # just your jobs
sacct -j <jobid> --format=JobID,State,ExitCode   # how a job finished
```

Put `--array` before the script name (`sbatch --array=1-96 script.sh`). Anything after the script name is passed to the script, not to sbatch.

### Writing a Slurm job script

Start the file with this header (don't put any code above it):

```bash
#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=YourJobNameHere
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
```

Run `mkdir -p logs` before submitting.

<!-- TODO (Reference genome): the annotation listed here is a .gff3 but the working annotation in use is a filtered .gtf; add the assembly citation; rename 02.1_runBuildCtark1Genome_apline.sh ("apline" typo) -->
## Reference genome

- **Assembly:** *Culex tarsalis* CtarK1 [add source/citation]
- **Annotation:** `Culex-tarsalis_knwr_BASEFEATURES_CtarK1.gff3`
- **Derived files:** `genes_only.gff3`, `CtarK1_TSS.bed`, and `Culex-tarsalis_knwr_CtarK1.bed`, built with the steps in `Bed_file_instructions.txt` and `Gene_only_bed_instructions.txt`
- The HISAT2 index is built in Step 2.1 and written to `00_Cxt_Indexed_Genome/`

## Data

Raw sequencing reads are not included in this repository.

## Notes

- Large outputs (BAMs, bigWigs, matrices) and the `logs/` folders are not tracked in this repo.
- Contact: [Shawn / shawn.yates@colostate.edu]

## TODO

- [ ] Finish updating the scripts through Step 4
  - [ ] Step 2: `Gvis_Function`, `HOW_TO_MAKE_FASTQTRIMMED_FILE`
  - [ ] Step 3: any other MACS2 scripts
  - [ ] Step 4: the representative-replicate scripts (`04.2_Representative_Replicates`, `04.3_Representative_Replicates_heat_maps`)
- [ ] Add the R scripts to the repo and document how to run them
- [ ] Update and document Step 5 (HOMER)
- [ ] Rebuild the lost conda environments (`fastp_env`, `macs2_env`, `deeptools_kernel_v2`) in `/projects`, export `.yml` files for all environments, and update the conda lines in the scripts
- [ ] Confirm which environment the Step 2 scripts use (`rnaPseudo_clean`)
- [ ] Back up pipeline results and raw reads off scratch
- [ ] Switch Steps 3 and 4.0 to replicate-matched input controls when all input replicates exist
- [ ] Fix the relative `Cxt_Genome` path in `02.1`, and rename the "apline" file
- [ ] Use one `--mail-user` address across all scripts
- [ ] Fix the Study overview (H3K9me vs H3K9me3, what BF and SF are) and add the genome citation
