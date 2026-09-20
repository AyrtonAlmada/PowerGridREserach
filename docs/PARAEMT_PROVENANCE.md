# Pinning ParaEMT correctly

The public repository is https://github.com/NatLabRockies/ParaEMT_public.
A link to `tree/main` identifies a project, not the version used by an experiment.
The case files supplied for this release do not contain enough information to
uniquely identify that historical Git revision.

The public history inspected during this task displays `d79d735` as a short tip
identifier, but **that is not asserted to be your run commit** and is not filled
into the publication's version field. Recover the full SHA from the checkout
that actually generated the calibration set.

## In the original local checkout

```sh
git -C /path/to/ParaEMT_public rev-parse HEAD
git -C /path/to/ParaEMT_public status --short
git -C /path/to/ParaEMT_public diff --binary HEAD
```

For your Windows workflow, use the actual folder, for example:

```powershell
conda activate paraemt
python scripts/capture_paraemt_provenance.py `
  "C:\Users\pablo\ParaEMT_public" `
  "provenance\paraemt\calibration-run" `
  --run-label "wecc240-calibration" `
  --extra-file "main_step1_simulation.py"
```

The example path comes from the previously described local workflow; replace it
with the folder that really produced your data. Running this on a new clone only
records the new clone. It does not establish provenance of existing CSVs.

The script records the full SHA, working-tree status and tracked patch, untracked
file names, Python/pip environment, file hashes, and selected case inputs. Inspect
the patch before release. Untracked files are listed, not automatically published;
add your actual fault/switching scripts and any other required files explicitly.
Hash the initialization snapshot with `--extra-file` when that snapshot was used.
Record EMT integration step, output decimation, voltage/phase extraction, clearing
logic and any modified device settings in the run manifest as well.

## Reproduce a pinned checkout

```sh
git clone https://github.com/NatLabRockies/ParaEMT_public.git .external/ParaEMT_public
git -C .external/ParaEMT_public checkout --detach FULL_40_CHARACTER_RUN_SHA
# If the run used local modifications, apply the archived patch from that record.
git -C .external/ParaEMT_public apply /absolute/path/to/tracked_changes.patch
```

Replace the SHA token with the recorded value; do not literally run the placeholder.
Keep a permanent source link of the form:

```text
https://github.com/NatLabRockies/ParaEMT_public/tree/FULL_40_CHARACTER_RUN_SHA
```

A Git submodule is an alternative:

```sh
git submodule add https://github.com/NatLabRockies/ParaEMT_public.git external/ParaEMT_public
git -C external/ParaEMT_public checkout --detach FULL_40_CHARACTER_RUN_SHA
git add .gitmodules external/ParaEMT_public
```

Commit the gitlink and retain any local patch separately. A submodule pinned to an
arbitrary new revision does not document the old run.

## If the original code was downloaded as a ZIP

Do not fabricate a commit. Keep the ZIP/source snapshot, checksums, download record,
and local modifications. Compare to upstream revisions if possible, or rerun the
reference generation with a genuinely pinned checkout. The provenance-capture script
still hashes relevant files without Git, but records the commit as unknown.

## Licenses and case acquisition

Keep the upstream BSD-3-Clause notice from the recorded revision and the notices
inside individual case files. The bundled source case files are private preparation
copies supplied by the author; check their redistribution conditions before upload.
If a case cannot be redistributed, replace its public copy with pinned acquisition
instructions, its expected checksum, and a clear restriction in the DAS. Do not
silently label derived case files as covered by your own code's license.

Files explicitly supplied through `--extra-file` are copied under the captured
record's `extra_files/` directory, as well as hashed. Review those files and the
Git patch for private paths or credentials before release. Automatically discovered
case/source files are hashed; archive their actual bytes with the code/data release.

When a historical checkout has no recoverable commit, do not assign today's
upstream commit. Archive the actual source bundle and record its SHA-256. The
calibration indexer accepts this explicitly incomplete provenance only with
`--allow-unknown-commit`, `paraemt_commit=unknown`, and a
`source_archive_sha256` field. This option records the limitation; it does not
verify that the bundle was the historical generating code.
