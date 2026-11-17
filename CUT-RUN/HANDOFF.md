# Session handoff

Paste this at the start of a new chat, with `CONTEXT.md`. Update it at the end
of each session; keep only the current state, not history (git has the history).

Last updated: 2026-11-16 (filled from earlier notes; correct anything stale)

## Current state

- Scripts reviewed and updated: Step 1 (all), 02.1, 02.2, 02.3, 02.4, 03.1,
  Step 4 `00_Create_Bigwigs`, 04.1, 04.1.2.
- All scripts after 02.1 load conda from `/projects` (miniforge3 + conda_envs).
- Pipeline results and raw reads are backed up off scratch.
- Latest pushes to `main`: e0b1e64, 8b9be71.

## Changed this session

- (fill in)

## Open questions

- 02.2.5 `MultiQC.sh` still needs review.
- Do region counts match across BF/MP/SF/ZH after 04.1? If not, rerun without
  `--skipZeros` before the `cbind` in 04.1.2.
- Are `fastp_env`, `macs2_env` and `deeptools_kernel_v2` rebuilt in `/projects`
  with `.yml` files in `CondaEnvs/`?
- `update_readme.sh` is untracked: commit it or add it to `.gitignore`?

## Next steps

1. Finish script review through Step 4, then Step 5 (HOMER).
2. Add the R scripts to the repo.
3. Regenerate the README with `update_readme.sh` once all scripts are done;
   it must mention the temporary D1 input pairing and the name `02.3_runBam_array.sh`.
4. Point scripts at `config.sh`, one step at a time, testing each.

## Do not touch

- `Gvis_Function` (ignored for now).
- README, until all scripts are done.

## Paste when asking for help

- The script, in full.
- `head` of its input list or sample sheet.
- The `.err` log and `sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS`.
