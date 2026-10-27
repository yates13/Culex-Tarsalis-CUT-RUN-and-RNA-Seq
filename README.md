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

## Pipeline

Run the steps in order. Each folder has its own scripts, numbered in run order.

| Step | Folder | What it does | Status |
|------|--------|--------------|--------|
| 1 | [`CUT-RUN/01_Cxt_FastP`](CUT-RUN/01_Cxt_FastP) | Read trimming with fastp, sample renaming, MultiQC summary of the fastp reports | Scripts updated |
| 2 | [`CUT-RUN/02_Cxt_Alignment`](CUT-RUN/02_Cxt_Alignment) | HISAT2 index and alignment, BAM conversion, replicate merging, MultiQC | In progress (02.1 and 02.2 updated) |
| 3 | [`CUT-RUN/03_MACS2_Peaks_New_Peaks`](CUT-RUN/03_MACS2_Peaks_New_Peaks) | Peak calling with MACS2 | Not yet updated |
| 4 | [`CUT-RUN/04_Deeptools`](CUT-RUN/04_Deeptools) | bigWigs, signal matrices, heatmaps and profiles with deepTools | Not yet updated |
| 5 | [`CUT-RUN/05_HomerAnnotation`](CUT-RUN/05_HomerAnnotation) | Peak annotation with HOMER | Not yet updated |

To view tracks in IGV, see [`IGV_ChIP-Seq_Instructions.pdf`](CUT-RUN/IGV_ChIP-Seq_Instructions.pdf).

## Setup

Each step activates its own conda environment from `${SCRATCH}/conda_envs/`.

| Environment | Used in |
|-------------|---------|
| `fastp_env` | Step 1, trimming |
| `rnaPseudo_clean` | Step 1 (MultiQC) and Step 2 (HISAT2 index build and alignment) |
| `macs2_env` | Step 3, peak calling |
| `deeptools_kernel_v2` | Step 4, see `How_to_make_deeptools_kernel_env.md` |

<!-- TODO: add the Step 5 environment (homer_env) and confirm the Step 2 scripts after 02.2 use rnaPseudo_clean; confirm .yml files exist in the repo -->
If a `.yml` file exists for an environment (for example in `CondaEnvs/`), create it with:

```bash
conda env create --file YourCondaEnv.yml --prefix ${SCRATCH}/conda_envs/<env_name>
```

Scripts load conda like this. Set `SCRATCH` to your own scratch path first:

```bash
SCRATCH=/scratch/alpine/<your_username>
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh
conda activate ${SCRATCH}/conda_envs/<env_name>
```

## Running

1. Edit the paths at the top of each script to match your directory (and the `--mail-user` address in the `#SBATCH` header).
2. All scripts use `--partition=acpu` and `--qos=cpu-normal`, and write their Slurm logs to a `logs/` folder. Create it in the folder you submit from before submitting (`mkdir -p logs`). Slurm cannot create it, and a job whose log folder is missing fails without a log.
3. Submit each step with `sbatch <script>.sh` and check outputs before moving on.
4. See the `How_to_make_*` files in each folder for building sample and treatment lists.

## Step details

### Step 1: Trimming and QC (`CUT-RUN/01_Cxt_FastP`)

<!-- TODO: explain "Round1_Round2_With_Duplicates" (the pipeline output folder name) if it should appear here -->
Scripts, in run order:

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

### Step 2: Index and alignment (`CUT-RUN/02_Cxt_Alignment`)

Scripts reviewed and updated so far:

| Script | What it does | Environment |
|--------|--------------|-------------|
| `02.1_runBuildCtark1Genome_apline.sh` | Builds the HISAT2 index for CtarK1 from `Cxt_Genome/CtarK1.fa` and writes `CtarK1_index.*.ht2` to `Cxt_Cut_Run_Pipeline/00_Cxt_Indexed_Genome/`. Run once. | `rnaPseudo_clean` |
| `02.2_Cxt_AlignArray_Hisat.sh` | Aligns each trimmed read pair with HISAT2, one sample per SLURM array task. Writes `<sample>.sam` and `<sample>.hisat2.log` to `02.2_Cxt_HisatAligned_Cleaned/`. | `rnaPseudo_clean` |

Still to be reviewed and documented: `02.2.5_MultiQC_Hisat_aligned`, `02.3_runBamm_array.sh`, `02.4_BamMerge_With_replicate_counts_array.sh`, `Gvis_Function`, `HOW_TO_MAKE_FASTQTRIMMED_FILE`.

**1. Build the index (once).** Submit from the folder that contains `Cxt_Genome/`:

```bash
mkdir -p logs
sbatch 02.1_runBuildCtark1Genome_apline.sh
```

**2. Align.** Submit from `02_Cxt_Alignment`. First list the R1 trimmed files. Note that array numbering here starts at 0, so the array is `0-(N-1)`:

```bash
mkdir -p logs
ls <fastp_out_dir> | grep '_R1_' | grep 'trimmed' > FastqTrimmed.txt
N=$(wc -l < FastqTrimmed.txt)
sbatch --array=0-$((N-1))%20 02.2_Cxt_AlignArray_Hisat.sh
```

HISAT2 runs with `--no-spliced-alignment` (CUT&RUN is genomic, not RNA), `--no-mixed` (only pairs where both mates align) and `--mm` (memory-mapped index).

An array larger than the list makes the extra tasks quietly do nothing, so a successful job does not prove every sample ran. Check the results:

```bash
ls <OUT_DIR>/*.sam | wc -l                       # should equal N
sacct -j <jobid> --format=JobID,State,ExitCode   # look for FAILED or TIMEOUT
```

A task that hits the time limit leaves a truncated `.sam`. Rerun it with `sbatch --array=<task> 02.2_Cxt_AlignArray_Hisat.sh`.

## How to use Slurm

Common commands: sbatch, srun, scancel, sacct, squeue, sinfo

### Common Slurm commands

```bash
sbatch [options] script [args]                   # submit a job, e.g. sbatch 04_deeptools.sh path/to/bams/folder
scancel 12345                                    # cancel job 12345
scancel {31415..31425}                           # cancel a range of sequential job IDs
squeue | grep <your_username>                    # just your jobs
sacct -j <jobid> --format=JobID,State,ExitCode   # how a job finished
```

### Writing a Slurm job script

Create the file with `nano YourNewJob.sh`. Start it with this header (don't put any code above it):

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

<!-- TODO (Reference genome): the annotation listed here is a .gff3 but the working annotation in use is a filtered .gtf; add the assembly citation; check the "apline" spelling in the 02.1 script name -->
## Reference genome

- **Assembly:** *Culex tarsalis* CtarK1 [add source/citation]
- **Annotation:** `Culex-tarsalis_knwr_BASEFEATURES_CtarK1.gff3`
- **Derived files:** `genes_only.gff3`, `CtarK1_TSS.bed`, and `Culex-tarsalis_knwr_CtarK1.bed`, built with the steps in `Bed_file_instructions.txt` and `Gene_only_bed_instructions.txt`
- The HISAT2 index is built in Step 2 (`02.1_runBuildCtark1Genome_apline.sh`) and written to `00_Cxt_Indexed_Genome/`

## Data

Raw sequencing reads are not included in this repository.

## Notes

- Large outputs (BAMs, bigWigs, matrices, logs) are not tracked in this repo.
- Contact: [Shawn / shawn.yates@colostate.edu]

## TODO

- [ ] Finish updating the scripts through Step 4
  - [ ] Step 2: `02.2.5_MultiQC_Hisat_aligned`, `02.3_runBamm_array.sh`, `02.4_BamMerge_With_replicate_counts_array.sh`, `Gvis_Function`, `HOW_TO_MAKE_FASTQTRIMMED_FILE`
  - [ ] Step 3: MACS2 peak calling scripts
  - [ ] Step 4: deepTools scripts
- [ ] Add the R scripts to the repo and document how to run them
