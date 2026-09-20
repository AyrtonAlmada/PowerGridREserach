# Validation status of this release candidate

This file reports checks actually performed during package preparation. It is not
an assertion that the manuscript was reproduced or that the new solver is calibrated.

## Performed

- Read the supplied Julia and CEM source and compared their stated numerical
  operations with the refactored implementation by inspection.
- Checked all seven original code/input attachments against the archived copies:
  all copied input bytes are identical (SHA-256 records are included).
- Executed the processed-network audit: 243 bus rows, 248 branch rows, 70 connected
  components, effective inertia 1, damping 0.5, and branch limit 4 in supplied units.
  The supplied JSON contains 243 buses, 329 lines, and 122 transformers.
- Compiled the three Python helper scripts to Python bytecode as a syntax check.
- Tested calibration indexing with a complete **temporary test-only** 25 x 46
  fixture, rejection of incomplete grids and existing output files, and explicit
  handling of unknown Git history with a source-bundle hash.
- Tested provenance capture on a temporary local Git repository, including the
  full commit, dirty tracked patch, untracked-file listing, explicit extra-file
  copy, and overwrite protection.
- Independently checked in Python the phase-difference matrix identity G²=0,
  global-phase invariance of that noise, and known importance-weight/SE/ESS values.

The temporary test trajectories were artificial software-test fixtures, were
removed, and are **not** calibration data. No 1,150-reference dataset is supplied.
The machine-readable preparation log is `provenance/preparation_checks.json`.

## Not performed here

**Julia was not installed in the preparation environment.** Therefore the Julia
module was not compiled or executed, `test/runtests.jl` was not run, and no WECC
surrogate/CEM/ParaEMT trajectories were generated. Python algebra checks are not a
substitute for running the actual Julia implementation. No runtime, parameter-fit,
mean-square stability, rare-event precision or reproduction claim is certified.

## Required local checks

From the repository root:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl
julia --project=. examples/wecc240/single_trajectory.jl
```

Save the output and the actual generated `Manifest.toml`. Start with a small,
clearly identified rerun. Check the processed topology and mappings, compare the
deterministic solution against an independent ODE calculation, and test the chosen
stochastic backend and time-step convergence before running the full ensemble.

The pathwise Heun implementation is a newly introduced alternative. It requires
new validation and, where applicable, recalibration. Keeping the old alpha and
sigma numbers in a configuration does not establish that they remain optimal.

A final paper release additionally needs the actual raw-to-processed transformation,
EMT fault/reclosure and projection scripts, calibration fitting routine, figures and
tables generated from saved samples, archived reference trajectories, and confirmed
run provenance. The provided files do not fill missing experimental records.
