# PowerGridREserach — EMT-informed dynamic contingency screening

Julia tools for stochastic electromechanical simulation, cumulative line-overload
indicators, and cross-entropy importance sampling. The project uses high-fidelity
EMT trajectories to inform a reduced model for repeated dynamic contingency
assessment.

**Project repository:** [AyrtonAlmada/PowerGridREserach](https://github.com/AyrtonAlmada/PowerGridREserach)  
**EMT reference software:** [NatLabRockies/ParaEMT_public](https://github.com/NatLabRockies/ParaEMT_public)

Historical Israel studies and WECC240 experiments are stored separately. Each
released experiment should identify its input files, model parameters, solver
settings, random seeds, and source-code revision.

## Release status

This README specifies the phase-endpoint noise model below and documents
**`:stratonovich_heun`** as the stochastic trajectory backend.

**Implementation alignment:** this README update does not change the Julia
module. Its noise-construction routine and experiment configuration must be
updated to produce the stated lower-right-block matrices before running this
model. Selecting the integration backend alone does not select or change the
noise matrix. The matrix checks in the example below prevent an incompatible
configuration from running unnoticed.

The 1,150 calibration reference trajectories are prepared for deposition but
are not included in the release-candidate package. The exact ParaEMT revision
used to generate them remains to be recorded. The package's Julia tests have
not been executed in the preparation environment; local test results and the
resolved Julia environment must accompany a reproducible release.

## Workflow

```text
ParaEMT reference trajectories
              ↓
EMT-informed parameter calibration
              ↓
Stratonovich surrogate trajectories
              ↓
Line-specific and global overload durations
              ↓
Importance-weighted probabilities and contingency rankings
```

The surrogate parameters describe the physical reduced model. The CEM proposal
parameters describe how contingencies are sampled. These are separate inputs:
adapting a proposal does not recalibrate the surrogate.

## Repository organization

```text
PowerGridREserach/
├── README.md
├── Project.toml
├── PowerGridsFunctions3.jl          # convenience loader
├── src/PowerGridsFunctions3.jl      # surrogate, overload metrics, and CEM
├── config/wecc240.toml             # experiment settings
├── examples/wecc240/
├── scripts/
│   ├── run_wecc_cem.jl
│   ├── capture_paraemt_provenance.py
│   ├── audit_inputs.py
│   └── index_calibration.py
├── test/runtests.jl
├── data/wecc240/
│   ├── upstream/                   # original ParaEMT case inputs
│   ├── processed/                  # surrogate branch and bus tables
│   ├── metadata/                   # mappings, manifests, and input audit
│   └── calibration/
│       ├── emt/                    # reference trajectories
│       └── surrogate/              # reduced-model trajectories
├── results/wecc240/                # separate directory for each run
├── legacy/
│   ├── israel/                     # preserved historical case studies
│   └── source/                     # original supplied source files
├── provenance/
├── docs/
└── manuscript/
```

Follow [MIGRATION.md](docs/MIGRATION.md) to move the existing Israel folders
without discarding their history or results. The release-candidate package does
not contain a replacement copy of those historical results.

## Network and state variables

The reduced state uses the ordering

```math
X_t = \begin{bmatrix}
\omega_t \\
\theta_t
\end{bmatrix}
\in\mathbb{R}^{2n},
```

where `n` is the number of modeled nodes, `omega` is the frequency-deviation
state, and `theta` contains phase angles in a consistently defined reference
frame. The nominal drift is

```math
A=\begin{bmatrix}
-M^{-1}D & -M^{-1}L\\
I & 0
\end{bmatrix},
\qquad
b(z) = \begin{bmatrix}
M^{-1}p(z) \\
0
\end{bmatrix}.
```

Here, `M` and `D` are diagonal effective inertia and damping matrices, `L` is
the weighted network Laplacian, and `p(z)` is the power-injection vector for
operating state `z`. These are the effective parameters of the supplied reduced
model, not automatically the detailed machine parameters in the EMT model.

The processed inputs contain **243 bus rows and 248 branch rows**. The supplied
ParaEMT JSON contains **243 buses, 329 lines, and 122 transformers**. The processed
branch graph has **70 connected components**; its topology is not silently
completed by adding missing elements. The input audit records these differences
and the supplied effective inertia and damping values.

Original bus IDs and state indices are distinct. Use `bus_map.csv` and
`branch_map.csv` when matching trajectories, fault labels, and original case
files. Input schemas are documented in [DATA_DICTIONARY.md](docs/DATA_DICTIONARY.md).

## Piecewise stochastic surrogate

For a faulted line `e = (i,j)`, let its endpoint set be
$`\partial e=\{i,j\}`$. With pole opening at time zero, the model is


$$
\begin{cases}
dX_t = \bigl[A_e(\alpha)X_t+b(z)\bigr]\,dt + \sum_{r\in\partial e}G_{e,r}(\sigma_{e,z,r})X_t\circ dW_t^r, & 0\leq t<T_{\mathrm{op}} \\
\dot{X}_t = AX_t+b(z), & T_{\mathrm{op}}<t\leq T_{\mathrm{hor}}
\end{cases}
$$


The state is continuous at reclosure. The independent Wiener drivers act only
during the open-phase interval. The affected line has effective coupling
$`\alpha\beta_e`$ during that interval and nominal coupling $`\beta_e`$ afterward.

The code's time arguments are:

| Argument | Meaning |
|---|---|
| `T1` | Pole-opening time |
| `T2` | Reclosure time |
| `T2 - T1` | Open-phase duration |
| `T3` | Final evaluation or screening time |
| `saveat` | Spacing of saved trajectory observations |
| `dt` | Maximum internal integration step during stochastic propagation |

When reclosure occurs after `T3`, the altered network remains active throughout
the observed interval. The original duration draw remains unchanged in the
importance-sampling likelihood ratio.

### Localized phase-endpoint noise

For the common-amplitude model, the endpoint amplitudes satisfy
$`\sigma_{e,z,i}=\sigma_{e,z,j}=\sigma`$. The noise matrix is

```math
G_{e,r}(\sigma_{e,z,r})=
\begin{bmatrix}
0 & 0\\
0 & \sigma M^{-1}\mathbf{e}_r\mathbf{e}_r^{\mathsf T}
\end{bmatrix},
\qquad r\in\partial e,
```

where $`\mathbf{e}_r`$ is the `r`th canonical vector in $`\mathbb{R}^n`$ and
$`M=\mathrm{diag}(m_1,\ldots,m_n)`$, with $`m_r>0`$. For independent endpoint
amplitudes, replace `sigma` by the corresponding $`\sigma_{e,z,r}`$ in each matrix.

With the state ordering above, each matrix is diagonal and rank one for nonzero
amplitude, with its only nonzero entry at

```math
[G_{e,r}]_{n+r,n+r}=\frac{\sigma}{m_r}.
```

**The noise acts on phase, not frequency.** Its contribution to the phase equation
is

```math
d\theta_t=\omega_t\,dt+
\sum_{r\in\partial e}\frac{\sigma}{m_r}
\mathbf{e}_r\theta_r(t)\circ dW_t^r.
```

Consequently, `omega` is the frequency-deviation state supplying the phase drift;
it is not the pathwise derivative of the noisy phase process.

This model depends on the phase reference: multiplying absolute phase coordinates
does not preserve invariance under adding a common constant to every phase.
Calibration data, initial conditions, and surrogate runs must use the same fixed
reference convention. The existing initialization chooses the minimum-norm
pre-fault solution of $`L\theta_0=p(z)`$ when that equation is consistent, with
zero initial frequency deviations. Do not re-center or wrap the saved phase
trajectories to change this convention after simulation.

The diffusion coefficient $`\sigma/m_r`$ must have units
$`\mathrm{s}^{-1/2}`$. If `M` is dimensionless in the adopted normalization, `sigma`
also has those units; otherwise its units must include the inertia scaling.
Record that normalization alongside the fitted amplitudes.

### Stratonovich integration

The stochastic trajectory backend is selected explicitly:

```julia
backend = :stratonovich_heun
```

The predictor and corrector use the **same Wiener increments** within each step.
For an open-phase step of length `h`, the scheme is

```math
\widetilde X=X_n+h f(X_n)+\sum_rG_{e,r}X_n\Delta W_n^r,
```

```math
X_{n+1}=X_n+\frac h2\bigl[f(X_n)+f(\widetilde X)\bigr]
+\frac12\sum_rG_{e,r}(X_n+\widetilde X)\Delta W_n^r,
\qquad \Delta W_n^r\sim\mathcal{N}(0,h),
```

where $`f(X)=A_e(\alpha)X+b(z)`$. The integrator lands on opening and reclosure
times rather than taking a stochastic step across a switching event. After
reclosure, propagation uses the nominal affine dynamics without further noise.

`dt` and `saveat` have different roles. For example, `dt=0.001` and `saveat=0.01`
use internal steps no larger than 0.001 s while saving observations every 0.01 s,
with switching times also retained. Report both settings and check step-size
convergence before interpreting calibrated parameters or overload probabilities.

### Deterministic moment propagation

The equivalent Itô drift during the open-phase interval is

```math
A_{e,\mathrm I}=A_e(\alpha)+\frac12\sum_{r\in\partial e}G_{e,r}^{\,2}.
```

For the specified matrices, this correction belongs to the phase--phase block.
It changes the mean dynamics but is not, in general, a scalar change in the line
coupling parameter `alpha`.

The mean $`\mu=\mathbb{E}[X_t]`$ and covariance $`\Sigma`$ satisfy

```math
\dot\mu=A_{e,\mathrm I}\mu+b(z),
```

```math
\dot\Sigma=A_{e,\mathrm I}\Sigma+\Sigma A_{e,\mathrm I}^{\mathsf T}
+\sum_{r\in\partial e}G_{e,r}
(\Sigma+\mu\mu^{\mathsf T})G_{e,r}^{\mathsf T}.
```

After reclosure, the drift is `A` and the diffusion terms are absent. The function
`moment_rhs(...; convention=:stratonovich)` evaluates these right-hand sides for
the matrices supplied to it. It is not the calibration optimizer itself.

The Itô correction is used when forming Itô moment equations. It must not be added
a second time to the Stratonovich drift supplied to the Heun integrator.

## Installation and model checks

From the repository root:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/runtests.jl
```

Retain the actual generated `Manifest.toml` with the tested release. The existing
tests must also be updated to check the phase-endpoint matrices specified here;
a passing test of a different noise matrix is not validation of this model.

The matrix-construction requirement is straightforward: for each endpoint `r`,
create a `2n`-by-`2n` zero matrix and set its `(n+r,n+r)` entry to `sigma/m_r`.
This construction must be used by both `prepare_contingency` and the contingencies
created within `make_surrogate_score`.

The example below assumes that update has been made. Its explicit matrix check
stops execution if the builder still returns a different noise structure.

```julia
using CSV, DataFrames, Random
include("PowerGridsFunctions3.jl")

case = load_case(
    "data/wecc240/processed/branch240E3.csv",
    "data/wecc240/processed/bus240E3.csv"
)

system = prepare_system(case.df, case.df2; injection_scale=1.02)
line_row = 1
alpha = 0.516
sigma = 1.5

# Requires the phase-endpoint matrix construction documented above.
contingency = prepare_contingency(system, line_row; alpha, sigma)

n = system.n
endpoints = (Int(case.df.From[line_row]), Int(case.df.To[line_row]))
expected_G = [zeros(Float64, 2n, 2n) for _ in endpoints]
for (k, r) in enumerate(endpoints)
    expected_G[k][n+r, n+r] = sigma / system.df2.Inertia[r]
end

@assert length(contingency.G) == length(expected_G)
for k in eachindex(expected_G)
    @assert isapprox(contingency.G[k], expected_G[k]; atol=1e-12, rtol=1e-12) "Noise matrix does not match the README model."
end

DFX = simulate_trajectory(
    contingency;
    T1=0.0,
    T2=0.7,
    T3=2.5,
    saveat=0.01,
    dt=0.001,
    backend=:stratonovich_heun,
    rng=Xoshiro(2026),
    include_frequency=true
)
```

These parameter values specify the example run; changing the model or integration
method does not establish that an earlier calibration remains valid. Match each
reported parameter set to its actual model definition and validation records.

The output contains phase columns `x1`, ..., `xn`, optional frequency columns
`omega1`, ..., `omegan`, and `Time`. Reclosure is saved once. When changing the
operating point, scale the nominal power-injection vector and recompute the
initial equilibrium; use that same scaled operating state in all intervals.

## Dynamic overload indicators

For a monitored line `{i,j}`, the reduced flow proxy is

```math
p_{ij}(t)\approx\beta_{ij}[\theta_i(t)-\theta_j(t)].
```

The line score and global score are

```math
S_{ij}=\int_0^{T_{\mathrm{hor}}}
\mathbb{I}\{|p_{ij}(t)|>\overline p_{ij}\}\,dt,
\qquad S=\sum_{\{i,j\}\in\mathcal{E}_m}S_{ij}.
```

The implementation uses a left-endpoint integration rule on the actual saved time
intervals. `Sij` is cumulative time above a limit; `S` is accumulated line-overload
time across the monitored set. Neither quantity is a peak-flow value or a direct
temperature calculation. The endpoint of the saved grid is not counted as an
additional time interval.

The primary output uses nominal monitored-line coefficients. Large phase
separations limit the validity of the linearized flow proxy.
`phase_difference_diagnostics` records per-path maxima and temporal exceedance
fractions; these are not ensemble probabilities. Severe cases require appropriate
nonlinear or EMT-level assessment.

## Cross-entropy importance sampling

The implementation uses unweighted elite-based adaptation: 250 samples per
iteration, 20 iterations, 25 elites, and retention factor `kappa=0.9`. It learns
categorical fault probabilities `phi` and one shared exponential duration **rate**
`r`. This factorized proposal does not learn line--duration dependence and is not
a likelihood-ratio-weighted CE adaptation.

The final estimation law retains **both learned components**, `res.π` and `res.r`.
With uniform nominal fault selection and nominal duration rate `lambda0`,

```math
p_Z(e,\tau)=\frac1E\lambda_0e^{-\lambda_0\tau},\qquad
q(e,\tau)=\phi_e r e^{-r\tau},
```

```math
w_k=\frac{\lambda_0}{E\phi_{e_k}r}
\exp[(r-\lambda_0)\tau_k].
```

The implementation draws exponential durations using `randexp(rng)/r`. It uses a
fresh final batch with the same conditional stochastic-driver law as the nominal
model, so no additional driver likelihood factor is required. The current
assessment is conditioned on a fault; it does not sample a no-fault category.

All surrogate randomness must use the RNG passed to the scorer. Failed or
nonfinite solves raise an error rather than contributing a safe-event indicator.
Final estimates use ordinary importance weights, not self-normalized weights, and
are not clipped to `[0,1]`. Zero observed events do not establish zero risk.

In the experiment configuration, set:

```toml
backend = "stratonovich_heun"
```

Also align the matrix-construction settings with the phase-endpoint model before
running the driver. The existing configuration is not updated by editing this
README. Once the builder, configuration, and tests agree, the run command is:

```sh
julia --project=. scripts/run_wecc_cem.jl config/wecc240.toml results/wecc240/runs/wecc-run-001
```

Start with a small final batch. Use a new output directory for every run. The
driver uses separate adaptation and estimation streams and can evaluate several
thresholds from the same final batch. In that case, adaptation is counted once,
not once for each threshold row. The scorer and its cached contingency models
are serial; do not share the mutable cache or RNG across worker threads.

## Saved outputs and reproducibility

| File | Contents |
|---|---|
| `summary.csv` | Final sample count, estimate, SE, relative SE, ESS, maximum normalized weight, event ESS, hit count, and stage times |
| `final_samples.csv` | Sampled line, duration, score, driver seed, and log-weight |
| `pilot_samples.csv` | Adaptation samples and elite flags |
| `adaptation_history.csv` | Iteration-level scores, thresholds, and rate history |
| `adaptation_phi.csv` | Categorical proposal history |
| `proposal.csv` | Frozen line probabilities and shared duration rate |
| `run.toml` | Experiment settings, seeds, environment, timing scope, and source/input hashes |
| `bus_map.csv`, `branch_map.csv` | State-index and original-case mappings |
| `checksums.csv` | Checksums of saved run files |
| `Project.toml`, `Manifest.toml` | Project specification and resolved environment when available |

Preserve samples as well as seeds. Record the phase reference, mass normalization,
state ordering, endpoint amplitudes, opening and reclosing times, internal `dt`,
output `saveat`, and event-comparison convention. Keep the exact source revision
that generated each result.

Record adaptation and estimation runtimes separately, including compilation and
initialization policies. The current run driver times computation separately
from file serialization. Read [REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for the
run-record specification and [VALIDATION_STATUS.md](docs/VALIDATION_STATUS.md) for
the checks completed during preparation. Those records must be updated when the
phase-endpoint implementation is tested.

## Calibration archive and ParaEMT provenance

The base calibration set contains **25 selected faulted lines × 46 operating
scales**, with scales `0.75:0.01:1.20`, at the recorded switching and evaluation
settings. A separate sweep over 15 open-phase durations contains **17,250** cases.
Neither count determines the number of independent stochastic realizations.

Keep EMT references and surrogate outputs in different directories. A historical
filename beginning with `ParaEMTDF` does not identify the generating simulator.
Each manifest entry must record that source, the exact line mapping, operating
scale, time settings, projection/units, train/test assignment, and checksum.

After filling the calibration metadata with the actual trajectory records:

```sh
python scripts/index_calibration.py --metadata YOUR_METADATA.csv --data-root . --output data/wecc240/metadata/calibration_manifest.csv --require-1150
```

Capture ParaEMT provenance from the checkout and Python environment that generated
the reference trajectories:

```sh
python scripts/capture_paraemt_provenance.py /path/to/ParaEMT_public provenance/paraemt/calibration-run --run-label wecc240-calibration
```

This records the full commit, local modifications, input/source hashes, and Python
environment. The upstream `main` URL is a project link, not a historical version
identifier. Use the actual recorded commit in the release citation. See
[PARAEMT_PROVENANCE.md](docs/PARAEMT_PROVENANCE.md).

The final reproducibility archive must also include the raw-to-processed case
transformation, EMT switching and projection scripts, calibration optimizer, and
analysis scripts used to generate figures and tables. These were not all present
in the supplied source package.

## Data availability, citations, and licensing

Freeze the tested code as a versioned release and archive the associated reference
data with a manifest and checksums. Link the code release and data archive in both
directions. Fill the release identifiers in `manuscript/references.bib` and use
the after-release Data Availability statement only after the listed materials
are publicly accessible.

Related work and reference software are listed in the manuscript bibliography,
including the N1Plus studies and the ParaEMT publications. Preserve upstream
notices for the original case files and distinguish third-party licensing from
the license governing this project's own code. The repository's license template
must be completed by the rights holders before release.

This software supports research into contingency screening. It is not a protection
relay implementation or an operational security certification. Probability-zone
boundaries and overload-duration thresholds must be identified for each study;
they are not universal protection requirements.
