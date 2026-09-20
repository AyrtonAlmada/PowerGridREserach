# Preserve Israel results; add WECC240 without rewriting history

The public repository currently contains `0data`, `0examples`, `0src`, `1README`,
`2ExamplesOfSinglePhaseFaults`, and `2ExamplesOfThreePhaseFaults`. Keep the repository
name `PowerGridREserach` unchanged so existing citations continue to resolve.

## 1. Create a recoverable checkpoint

Work on a clone, not directly on your only working directory. Review and commit any
intended local changes before the checkpoint. Do not run a history rewrite, force
push, or `git reset --hard` to perform this migration.

```sh
git status
git switch -c reorganize-wecc240
git tag -a israel-before-wecc -m "Israel study before WECC240 reorganization"
git branch archive/israel-before-wecc
```

Choose a different tag/branch name if it already exists. Publish the checkpoint when
appropriate with `git push origin israel-before-wecc archive/israel-before-wecc`.
A tag and a branch preserve the old layout; moving files does not erase old commits.

## 2. Move the old tree intact

Keep the original numbered subfolders together, rather than mixing old and new
files by their generic names. This preserves many relative links in the notebooks.
Commands below use Git Bash (also available with Git for Windows).

```sh
mkdir -p legacy/israel
git mv 0data legacy/israel/0data
git mv 0examples legacy/israel/0examples
git mv 0src legacy/israel/0src
git mv 1README legacy/israel/1README
git mv 2ExamplesOfSinglePhaseFaults legacy/israel/2ExamplesOfSinglePhaseFaults
git mv 2ExamplesOfThreePhaseFaults legacy/israel/2ExamplesOfThreePhaseFaults
git mv README.md legacy/israel/README.md
```

Check each command against `git ls-tree --name-only HEAD`; do not move an unrelated
new folder into the Israel archive just because it has a similar name. After moving,
run the Israel notebooks with their original working directory under
`legacy/israel/`. Absolute paths or root-relative links still need updating in a new
copy; preserve the original notebook at the checkpoint.

## 3. Copy this release candidate into the repository root

Copy its `legacy/source/` directory alongside the migrated `legacy/israel/`;
do **not** replace the entire `legacy/` folder. The `legacy/source` files are the attachments supplied for
this update, not a substitute for the historical Israel source.

Add the new `src/`, `data/wecc240/`, `config/`, `examples/wecc240/`, `test/`, `docs/`,
`scripts/`, `provenance/`, and `manuscript/` folders. Put the new README and project
files at the root. No remote repository has been changed by preparing this package.

## 4. Verify before committing

```sh
python scripts/audit_inputs.py
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl
git diff --stat
git status
```

Save test outputs and the actual `Manifest.toml`. Do not migrate historical figures
into `results/wecc240` or label new solver outputs as reproductions of older runs.
Use `results/wecc240/runs/<run-id>` for new experiments. Identify legacy visual
examples separately from measurements supporting the current paper.

## 5. Deposit reference trajectories and archive

Keep 1,150 EMT references in `data/wecc240/calibration/emt/` or a linked DOI-bearing
archive. Use `data/wecc240/calibration/surrogate/` for the reduced-model sweep.
Commit manifests and checksums even if large CSVs reside in the data archive.

After all author checks, create a new release tag, e.g. `wecc240-prx-v1.0.0`, and
archive that specific code release. This tag name is a suggestion, not an existing
release. Cite the immutable code release and data version, not only a moving `main`.

Source inspected: https://github.com/AyrtonAlmada/PowerGridREserach
