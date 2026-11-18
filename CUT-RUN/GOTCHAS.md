# Gotchas: bugs already solved

Check here before debugging. Add an entry whenever a bug costs more than an
hour. Format: symptom, cause, fix.

## Cluster and environment

| Symptom | Cause | Fix |
| --- | --- | --- |
| Perl-based tool fails on a path | The `c832500103@colostate.edu` path contains `@`, which Perl reads as an array | Use the `.colostate.edu/c832500103` form everywhere (`config.sh`). For `hisat2_env_fixed`, also escape the `@` in `Config_heavy.pl` (see `RERUN_2026-11.md`) |
| `hisat2` exits 255 after about 30 s: `Global symbol "@colostate"` in `Config_heavy.pl` | The env's Perl config has the `@` path written into it; the no-`@` `PROJ=` alone is not enough | Escape it once (`RERUN_2026-11.md`); test with `perl -V`, not `perl -MConfig` |
| Job exits 127 at once: `hisat2`, `samtools` or `multiqc` not found | Script activated `rnaPseudo_clean2`, which has none of them | Step 2 uses `hisat2_env_fixed`, MultiQC uses `multiqc_env`; check with `ls <env>/bin/<tool>` before submitting |
| `conda create` hangs for an hour at "Executing transaction" (process in `D` state) | `/projects` stalls on the login node; quota was fine | Build envs with `sbatch --wrap 'conda create -y -p ...'` |
| A multi-line paste runs each line as its own command | The paste broke the line after each trailing backslash | Put long `git add` commands on one line |
| A file pasted into nano has every line doubled | Terminal paste added a blank line after each line | Check `wc -l`; if it is 2N-1, keep the odd lines: `awk 'NR % 2 == 1'` |
| Same file appears under two paths | Alpine mounts scratch under two different path strings | Use only the `config.sh` form; never mix them in one file list |
| `hisat2 --version` segfaults | Login node CPU differs from compute nodes; the env is fine | Test on a compute node with `srun` |
| Env cannot be activated by name | Env was created with `--prefix` | Activate by full path (`load_env`) |
| Env directory exists but tools are missing | Scratch purge removed files after 90 days | Keep envs in `/projects`; rebuild from the `.yml` |
| `bedtools` fails with a libbz2 error | Broken `libbz2.so.1.0` symlink in that env | Rebuild the env, or use an env with the real library |
| Script continues after `command not found` | No `set -e` | `set -euo pipefail` in every script |
| Script dies at `wc -l` | `set -e` plus `wc` on a missing file | Check the file exists first (`need`) |

## Silent wrong results

| Symptom | Cause | Fix |
| --- | --- | --- |
| Zero counts for every gene | Contig names differ between BAM and GTF | Compare names from `samtools idxstats` and the GTF before counting |
| Combined deepTools matrix is huge, with only the first group's labels | `computeMatrixOperations rbind` stacks regions; `cbind` stacks samples | Use `cbind` to combine groups; inspect the matrix JSON header |
| `cbind` output has NA or dropped regions | Regions matched by name to the first matrix; `--skipZeros` changed the lists | Compare region counts per group; rerun without `--skipZeros` |
| Join returns nothing | GMAP IDs end in `.mrna1` or `.path1`; FASTA headers do not | Strip suffixes before joining; assert the joined row count |
| One annotation copied onto many genes | Genes shared one representative transcript | Annotate each gene from its own sequence |
| "0% conservation" result | Missing backslash sent matches to stdout | Single-line commands; check counts |
| Wrong best BLAST hit | `sort -k1,1 -k12,12 -rn` used the wrong column and broke grouping | Single-pass awk or pandas, keyed on the real bitscore column |
| ~1.97M junk rows from a GTF filter | Filtered on feature type (`$3=="gene"`) | Filter on the source column (`maker`, `Geneious`, `VectorBase`) |
| Six mitochondrial genes vanish | `Mt` contig renamed after the real-contigs filter | Rename before filtering |
| Coordinates off by one | BED is 0-based, GTF is 1-based | `gtf_start = bed_start + 1`; comment every conversion |
| Gene IDs not matched by a regex | Pattern written for another species (`gene12345` vs `CPIJ007595-RA`) | Check the ID format in the actual file |

## Files and formats

| Symptom | Cause | Fix |
| --- | --- | --- |
| File corrupted after sorting | In-place `sort -o samefile` | Sort to a temp file, then `mv` |
| Truncated or corrupt SAM | Job hit the time limit mid-write; huge SAM output | Pipe to `samtools sort` and write BAM; run `samtools quickcheck`; size time from a test |
| R script crashes on file paths | Old Windows NAS paths in the metadata CSV | Replace the prefix; `stopifnot(all(file.exists(...)))` |
| Gviz pipeline fails on the metadata CSV | Columns must be `Population`, `Factor`, `Day`, `bamReads`, `Peaks` | Validate column names at the top of the script |
| `makeTxDbFromGFF()` fails | Mitochondrial genes lack `ID=`; function moved to `txdbmaker` | Exclude Mt genes first; install `txdbmaker` |
| GFF3 attributes corrupted | `rtracklayer::export()` | Filter as plain text with `readLines`/`writeLines` |
| Python `re.PatternError: bad escape` | Replacement string contained a backslash sequence | Pass a function (lambda) as the replacement in `re.sub` |
