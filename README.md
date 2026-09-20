# PowerGridREserach — EMT-informed dynamic contingency screening

Julia tools for reduced swing dynamics, cumulative line-overload indicators, and
cross-entropy (CE) importance sampling. The repository separates the historical
Israel studies from the WECC240 workflow and records the inputs and numerical
choices needed to audit a run.

**Repository:** https://github.com/AyrtonAlmada/PowerGridREserach  
**EMT reference software:** https://github.com/NatLabRockies/ParaEMT_public

## Release status

This update is a **release candidate**, not an archived reproduction of the paper.
The attached source and case files are included unchanged under `legacy/source/`
and `data/wecc240/`. The author reports that the **1,150 EMT calibration
trajectories** are ready for deposition; they are not included in this package.
The exact ParaEMT commit used for those trajectories has not yet been recovered.
No DOI, completed calibration rerun, or measured speedup is claimed here.

The Julia implementation was reviewed, but not executed in the preparation
environment, which did not contain Julia. Run `test/runtests.jl` locally before
launching an ensemble. `docs/VALIDATION_STATUS.md` records the checks actually run.

## What is implemented

- Import of the supplied processed branch/bus CSVs with explicit bus-ID maps.
- The `[omega; theta]` state ordering, nominal and open-phase swing matrices,
  equilibrium recomputation, and a single reclosure switch.
- Deterministic affine propagation and two explicitly distinguished stochastic
  backends, described below.
- Per-line cumulative overload durations, global scores, and phase-separation
  diagnostics calculated on the actual trajectory times.
- The supplied **unweighted-elite CE adaptation**, using a categorical line law
  and one shared exponential duration **rate**.
- Fresh final importance samples using **both** learned proposal components,
  with SE, relative SE, ESS, maximum normalized weight and complete sample logs.

A calibration optimizer, EMT fault-insertion runner, and figure-generation
pipeline were not included in the current source attachments. Add the actual
scripts that generated the paper, rather than describing those components as
already reproduced by this package.

## Repository layout

```text
PowerGridREserach/
├── README.md
├── Project.toml
├── PowerGridsFunctions3.jl          # convenience loader
├── src/PowerGridsFunctions3.jl      # supported simulation + CEM module
├── config/wecc240.toml             # explicit rerun settings
├── examples/wecc240/
├── scripts/
│   ├── run_wecc_cem.jl
│   ├── capture_paraemt_provenance.py
│   ├── audit_inputs.py
│   └── index_calibration.py
├── test/runtests.jl
├── data/wecc240/
│   ├── upstream/                  # original .raw, .xls, and .json
│   ├── processed/                 # branch240E3.csv, bus240E3.csv
│   ├── metadata/                  # checks, ID maps, manifest template
│   └── calibration/
│       ├── emt/                   # 1,150 reference trajectories: pending
│       └── surrogate/             # separate reduced-model outputs
├── results/wecc240/               # new run directories; no archived results here
├── legacy/
│   ├── source/                    # exact supplied .jl and CEM.txt
│   └── israel/                    # move existing Israel tree here; see migration
├── provenance/                    # input hashes, ParaEMT provenance template
├── docs/
└── manuscript/                    # DAS, acknowledgments and references
```

`legacy/israel/` is a migration destination. This package does not contain a copy
of your existing remote Israel results and does not move or overwrite them.
Follow [MIGRATION.md](docs/MIGRATION.md) before merging these files.

## Installation and first check

Use a dedicated Julia project; do not install packages from inside library code.
The declared compatibility is Julia 1.10/1.11, CSV 0.10, and DataFrames 1.6 or later
within major version 1. These are compatibility ranges, not a resolved lockfile.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl
julia --project=. examples/wecc240/single_trajectory.jl
```

On Windows PowerShell the same commands work; use double quotes around the Julia
expression when required by your shell. Commit the **actual generated
`Manifest.toml`** after successful tests and include it in the archived release.
Do not construct a lockfile by hand or run `Pkg.update()` when reproducing a release.

Interactive use:

```julia
using CSV, DataFrames, Random
include("PowerGridsFunctions3.jl")

case = load_case(
    "data/wecc240/processed/branch240E3.csv",
    "data/wecc240/processed/bus240E3.csv"
)
df, df2 = case.df, case.df2

# Deterministic smoke test. The scaled operating state is used consistently.
system = prepare_system(df, df2; injection_scale=1.02)
contingency = prepare_contingency(system, 1; alpha=0.516, sigma=0.0)
DFX = simulate_trajectory(contingency;
    T1=0.0, T2=0.7, T3=2.5, saveat=0.01, include_frequency=true)
Sij = line_overload_indicators(DFX, df)
S = sum(Sij.Sij)
```

Output phase columns are `x1`, ..., `x243` for the supplied inputs. Optional
frequency columns are `omega1`, ..., `omega243`; `Time` contains the actual
sampling instants. An original bus label is **not** its state index; use the
exported `bus_map.csv` and `branch_map.csv`.

## Important input facts

The supplied processed case contains **243 bus rows and 248 branch rows**.
The supplied ParaEMT JSON contains **243 buses, 329 lines, and 122 transformers**.
The processed branch graph has **70 connected components**, including isolated
nodes. This loader does not silently add missing transformers or collapse the
case into a connected network. The processed inputs use `m=1` and `d=0.5` at
every bus; those values are effective model inputs, not imported GENROU inertia
and damping values.

See [DATA_DICTIONARY.md](docs/DATA_DICTIONARY.md) and the machine-readable
`data/wecc240/metadata/input_audit.json`. The `.raw`/`.xls`/`.json` to processed-CSV
transformation script and its choices must be added by the author. A collection
of files alone does not document that transformation.

## Model and backend selection

For state `X = [omega; theta]`, the uploaded matrices correspond to

```math
A = \begin{bmatrix}-M^{-1}D & -M^{-1}L\\ I & 0\end{bmatrix},\qquad
b = \begin{bmatrix}M^{-1}P\\0\end{bmatrix}.
```

`T1` is pole opening, `T2` is reclosure, and `T3` is the screening/evaluation end.
The open-phase duration is `T2-T1`. The altered coupling is active only in that
interval. Draws whose reclosure lies beyond `T3` are observed in the open-phase
configuration throughout the remaining window; the original sampled duration
is retained in the likelihood ratio.

There is no implicit choice of stochastic solver when `sigma > 0`:

| Backend | Meaning | Reproducibility status |
|---|---|---|
| `:deterministic` | Affine matrix-exponential propagation; `sigma=0` | Supported deterministic path |
| `:source_random_map` | Uploaded `make_f1` random-evaluation formula, with explicit RNG | Historical-formula audit, **not an SDE path** |
| `:stratonovich_heun` | New pathwise Stratonovich predictor-corrector integration | Rerun/step-convergence/calibration required |

The uploaded source uses **phase-difference noise**,

```math
G_e=\begin{bmatrix}0&\sigma M^{-1}(e_i-e_j)(e_i-e_j)^T\\0&0\end{bmatrix}.
```

For this matrix `G_e^2 = 0`. In particular, its Stratonovich-to-Itô correction
does **not** renormalize alpha. The optional `:frequency_endpoints` structure is
the different diagonal noise model discussed in the manuscript; it must be chosen
explicitly and is not identified with the uploaded implementation. No claim is
made that `(0.516, 1.5)` remains calibrated after changing the noise structure or
solver. The code exposes the distinction instead of silently resolving it.

Example of a **new** pathwise run, not a reproduction of the old random map:

```julia
c = prepare_contingency(system, 1;
    alpha=0.516, sigma=1.5, noise_structure=:phase_difference)
path = simulate_trajectory(c;
    T1=0.0, T2=0.5, T3=2.5, saveat=0.01, dt=0.001,
    backend=:stratonovich_heun, rng=Xoshiro(2026))
```

`saveat` is the output grid; `dt` is the Heun internal step. Reclosure is included
exactly once. In `:source_random_map`, changing `saveat` changes the sequence of
independent random evaluations, not merely the display resolution.

`moment_rhs` supplies the linear-SDE mean/covariance equations with the declared
calculus convention. It is not an implementation of the manuscript's calibration
optimizer, and those SDE moments must not be attributed to `:source_random_map`.

## Overload definition

The code uses the strict comparison `abs(beta_ij * (theta_i-theta_j)) > F_ij` and
integrates the Boolean indicator with a **left-endpoint rule on actual time
intervals**. It does not count the final endpoint as an additional interval.
The global score is the sum of line durations, not a peak-flow or temperature
measure. The primary output uses nominal monitored-line coefficients as in the
manuscript; it does not implement a time-varying three-phase line-flow model.

Large phase differences invalidate the small-angle interpretation of this flow
proxy. `phase_difference_diagnostics` reports per-path maxima and time fractions;
these time fractions are not fractions of the contingency ensemble. No wrapping
or post hoc shrinking of phase differences is performed.

## Cross-entropy importance sampling

The source algorithm is retained as an **unweighted elite search**: 250 pilot
samples per iteration, 20 iterations, 25 elites, and retention factor `kappa=0.9`.
It is not a likelihood-ratio-weighted rare-event CE fit. It fits categorical
probabilities `phi` and one shared exponential **rate** `r`; it does not fit a
separate rate for every line. The factorized proposal does not learn dependence
between line and duration.

The final law is `phi = res.π`, `r = res.r`. With uniform nominal line selection,

```math
w_k=\frac{\lambda_0}{E\phi_{e_k}r}\exp[(r-\lambda_0)\tau_k].
```

All exponential sampling uses `randexp(rng)/r`. Ordinary IS is not self-normalized;
finite-sample estimates are never clipped to `[0,1]`. Final samples and seeds are
saved. Samples with failed/nonfinite solver scores raise an error rather than being
counted as safe. Zero hits produce an estimate of zero but an **unavailable** error
estimate, not a claimed zero-risk confidence statement.

```sh
julia --project=. scripts/run_wecc_cem.jl config/wecc240.toml results/wecc240/runs/my-new-run
```

Start with a small final batch while checking the configuration. That script
uses distinct adaptation and estimation streams and the same frozen proposal for
all requested thresholds. The adaptation budget is therefore counted **once for
that run**, not once per row of its summary table. Change the configuration and
output directory for a separate experiment.

The integration is serial; do not share its mutable cache or RNG across threads.
See [REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for run records and benchmark rules.
A seed is not a substitute for saved samples and a pinned Julia environment;
Julia's RNG stream can change across versions [Julia RNG documentation].

### Saved outputs

`summary.csv` contains `N_final`, `N_adapt`, `N_total`, `phat`, `se`,
`relative_se`, `ESS`, `max_normalized_weight`, `event_ESS`, `hits`, and stage times.
`final_samples.csv` contains durations, line indices, scores, driver seeds, and
log-weights. Pilot samples, elite flags, every categorical update, and the frozen
proposal are saved separately. `run.toml` records settings, source and input hashes,
Julia/BLAS/thread information, and the declared timing scope. File writing is not
included in the simulation/estimation timer. No old runtimes are reused as new
measurements.

## Calibration archive and version provenance

The base calibration archive is **25 selected faults x 46 operating scales**
(`0.75:0.01:1.20`) at its documented fixed switching/evaluation configuration.
A separate 15-duration sweep contains **17,250** cases. Neither count describes
the number of independent noise realizations unless the manifest says so.

The author must provide the exact 25-line list, train/validation/test membership,
EMT channel projection and units, physical settings, solver commit, and file hashes.
Do not infer these from plot labels. Keep EMT references and surrogate trajectories
in separate directories even when an older filename starts with `ParaEMTDF`.

Populate `calibration_manifest.template.csv` with real records, then run:

```sh
python scripts/index_calibration.py --metadata YOUR_METADATA.csv --data-root . --output data/wecc240/metadata/calibration_manifest.csv --require-1150
```

To pin ParaEMT, run the following in the **actual ParaEMT environment and checkout**:

```sh
python scripts/capture_paraemt_provenance.py /path/to/ParaEMT_public provenance/paraemt/run-record --run-label calibration-release
```

This records the full commit, local patch, file hashes and Python environment.
Do not substitute the current `main` commit for the version used in older runs.
See [PARAEMT_PROVENANCE.md](docs/PARAEMT_PROVENANCE.md).

## Archiving, availability and citation

Keep large reference trajectories in a versioned data archive, with a manifest and
checksums in this repository. Freeze the code with a Git tag and archive that exact
release with a DOI. Link the code release and dataset version in both directions.
GitHub documents [release archiving through Zenodo]; APS requests citations for
public data/software in its [Data Availability guidelines]. A DOI has not been
issued for this package.

The draft and after-release manuscript statements are in `manuscript/`; do not
use the after-release wording before the data and code are actually deposited.

Related papers:
- *Real-Time Stochastic Assessment of Dynamic N-1 Grid Contingencies*,
  arXiv:2510.18007.
- *Real-Time Dynamic N-1 Screening: Identifying High-Risk Lines and Transformers
  After Common Faults*, arXiv:2602.12293.
- M. Xiong et al., *ParaEMT: An Open Source, Parallelizable, and HPC-Compatible EMT
  Simulator for Large-Scale IBR-Rich Power Grids*, IEEE Transactions on Power
  Delivery 39, 911–921 (2024).
- M. Xiong et al., *An Open-Source Parallel EMT Simulation Framework*, Electric
  Power Systems Research 235, 110734 (2024).

Use the exact release/version citations in `manuscript/references.bib` when they
have been completed. Do not invent a DOI or claim the pending calibration archive
is already publicly available.

## Licensing and responsible use

The existing project README declares MIT licensing. Preserve the authors' licensing
choice and add the complete license text with the correct copyright holders before
release. `LICENSE.template` is deliberately not an executed license declaration.
ParaEMT is distributed under BSD-3-Clause terms; its notice is kept with the upstream
inputs. Preserve any additional case-specific notices. Confirm rights before
publicly redistributing third-party case files or Israel data; otherwise provide
pinned acquisition instructions and checksums instead.

This is research software, not a protection relay implementation or an operational
security certification. Classifications and overload thresholds must be explained
in the corresponding study, not treated as universal grid requirements.

[Julia RNG documentation]: https://docs.julialang.org/en/v1/stdlib/Random/#Reproducibility
[release archiving through Zenodo]: https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content
[Data Availability guidelines]: https://journals.aps.org/authors/data-availability-statements
