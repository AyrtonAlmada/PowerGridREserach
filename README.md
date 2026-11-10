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
**`:stratonovich_heun`** as the stochastic trajectory backend. It also records
the corrected line-disjoint calibration/testing partition and the distinct time
settings used for calibration, hierarchy assessment, and temporal diagnostics
in Paper 3. The protocol and fault identifiers below supersede earlier README
descriptions of that partition.

**Implementation alignment:** this README update does not change the Julia
module. Its noise-construction routine and experiment configuration must be
updated to produce the stated lower-right-block matrices before running this
model. Selecting the integration backend alone does not select or change the
noise matrix. The matrix checks in the example below prevent an incompatible
configuration from running unnoticed. The duration sampler, likelihood ratio,
configuration, and saved metadata must also agree on whether a duration parameter
is a mean or a rate; the manuscript convention is specified in
[Cross-entropy importance sampling](#cross-entropy-importance-sampling). This is
a documentation update, not a new execution or verification of the Julia code.

The ParaEMT calibration reference dataframes and the corresponding surrogate
dataframes are now available in the repository as ZIP archives. The reference
archives are stored in
[`data/wecc240/calibration/emt/`](data/wecc240/calibration/emt/), and the
surrogate archives are stored in
[`data/wecc240/calibration/surrogate/`](data/wecc240/calibration/surrogate/).
Each directory contains 14 ZIP archives. The exact ParaEMT revision
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
adapting a proposal does not recalibrate the surrogate. Calibration at reclosure
and subsequent post-reclosure assessment are also separate experiments; their
time settings must not be interchanged.

## WECC240 calibration and evaluation protocol

Paper 3 distinguishes parameter identification at reclosure from assessment of
the subsequent transient. Calibration uses **25 faulted lines and 46 pre-fault
operating states**, with **`T_op = T_ev = 0.8 s`**. The hierarchy comparison in
Figure 4 uses a **different set of 25 faulted lines**, with **`T_op = 1.5 s`**
and **`T_ev = 2.5 s`**, over the same 46 injection scales. Figures 9–10 instead
hold the injections at their nominal values while varying one time parameter.

Figure numbers in this section refer to the Paper 3 manuscript accompanying this
protocol. They do not refer to the earlier Israel or IEEE 118-bus studies.

| Study | Fault set | Operating states | Open-phase duration `T_op` | State-evaluation time `T_ev` or screening horizon `T_h` |
|---|---|---|---|---|
| Calibration/training | 25 training lines listed below | 46 scales, `0.75:0.01:1.20` | `0.8 s` | `T_ev = 0.8 s` |
| Calibration-hierarchy assessment, Figure 4 | 25 testing lines listed below | 46 scales, `0.75:0.01:1.20` | `1.5 s` | `T_ev = 2.5 s` |
| Horizon study, Appendix D2 / Figure 9 | 25 testing lines | Nominal injections; scale `1.0` | Fixed at `1.5 s` | `T_ev` varies |
| Duration study, Appendix D3 / Figure 10 | 25 testing lines | Nominal injections; scale `1.0` | Varies over `0–1.5 s` | Fixed at `T_ev = 1.5 s` |
| Illustrative large-excursion transient, Figures 3 and 11 | Recorded `9#4` case | The operating state of that recorded example | `0.5 s` | Figure 3 stops at reclosure; Figure 11 continues to `10 s` |
| Probabilistic screening | Manuscript catalogue of 329 initiating lines | Operating state specified for the screening run | Sampled from the nominal law or the importance-sampling proposal | `T_h = 10 s` |

### Operating-state family and purpose of the time settings

The operating-state family is generated by uniformly scaling the nominal
active-power injection vector:

```math
p(z_c)=c\,p^{(0)},\qquad c\in\{0.75,0.76,\ldots,1.20\}.
```

Here `c` is the injection scale, denoted by `lambda` in the manuscript's
operating-state formula. It is not the exponential duration mean `lambda0`.
The pre-fault equilibrium is recomputed for each scale. The topology and dynamic
parameters remain fixed, and each EMT–surrogate comparison uses corresponding
initial conditions, the same phase reference, and the same switching sequence.

The calibration target is the electromechanical state at **reclosure**, when the
open-phase interval ends. With `T_op = T_ev = 0.8 s`, no post-reclosure evolution
is included in the fitting interval. The reduced description starts at pole
opening and does not resolve the preceding electrical fault-clearing transient.
The target therefore refers to reclosure, not to a separately simulated breaker
operation or fault-inception event.

The hierarchy assessment is deliberately more demanding. Its open-phase duration
is extended to `1.5 s`, and the state is compared at `2.5 s`, one second after
reclosure, while the trajectories can still be oscillating. It tests transfer to
held-out fault locations and a longer disturbance rather than agreement only at
the calibration time. Both the training and Figure 4 assessment contain
`25 × 46 = 1,150` line/operating-state combinations. These are two distinct
scenario collections, not two subsets of one 1,150-scenario collection.

The illustrated oscillatory response in Figure 11(a) motivates this type of
post-reclosure assessment, but that example has its own `0.5 s` open-phase
duration. It is not the run that defines the Figure 4 timing.

### Training and testing fault identifiers

The partition is **by faulted line**. The two lists are disjoint, while their
46 injection scales are shared. Accordingly, the global-model test examines
held-out line locations within the same prescribed operating-state family; it
is not a test on a disjoint set of injection scales or on a different network.
The listed order is an identifier record, not a risk ranking or a numerical
branch-row index.

```julia
training_lines = [
    "X1#11",  "X2#4",   "X7#18",  "X10#11", "X10#19",
    "X10#80", "X11#14", "X11#19", "X13#14", "X13#19",
    "X13#56", "X13#80", "X18#19", "X18#33", "X22#27",
    "X27#28", "X27#51", "X28#51", "X33#36", "X34#35",
    "X34#37", "X34#49", "X34#56", "X34#80", "X35#49",
]

testing_lines = [
    "X1#10",   "X110#154", "X112#155", "X117#119", "X118#119",
    "X15#228", "X18#30",   "X2#227",   "X2#8",     "X20#24",
    "X22#31",  "X24#31",   "X26#30",   "X34#81",   "X38#39",
    "X38#43",  "X4#233",   "X4#234",   "X40#41",   "X53#55",
    "X53#56",  "X61#71",   "X64#65",   "X67#69",   "X85#86",
]

@assert length(training_lines) == length(unique(training_lines)) == 25
@assert length(testing_lines) == length(unique(testing_lines)) == 25
@assert isempty(intersect(training_lines, testing_lines))

injection_scales = collect(75:120) ./ 100
@assert length(injection_scales) == 46
```

Resolve these identifiers against `branch_map.csv` and the original case inputs
before selecting a line in the solver. Do not infer a processed row number from
the list position or change the spelling of a stored fault identifier.

### Calibration objective and model hierarchy

Parameter identification uses the phase-angle and frequency observations at the
two endpoints of the faulted line (`q = 0`). The first two surrogate moments are
propagated deterministically and compared with the reference information at
`0.8 s`. The covariance target is scatter across the **training operating
states** at fixed training line; it is not an ensemble of repeated random
ParaEMT runs at an identical operating point. Testing observations do not enter
the global fit or its covariance target.

| Model | Parameterization |
|---|---|
| `M0` | Fixed deterministic baseline: `alpha = 2/3`, `sigma = 0` |
| `M1` | One globally fitted deterministic coupling; `sigma = 0` |
| `M2` | One globally fitted coupling and one globally fitted stochastic amplitude |
| `M3` | Separate parameter pairs for individual faulted lines |
| `M4` | Separate parameter pairs for individual line–operating-state combinations |

The reported global pair is **`alpha = 0.516`, `sigma = 1.5`**, in the model's
specified normalization. It is retained in the global-model assessment and the
subsequent screening calculations. The hierarchy comparison uses the same
observation convention and Figure 4 timing for all model levels. The
line-specific and line/state-specific fits require their own fitting records;
their results should not be described as transfer of an unseen line's fitted
parameters without documenting how those parameters were obtained.

### Relative state error and the temporal studies

For a specified fault, operating state, and switching time, the manuscript's
relative state error is

```math
\mathrm{RSE}_{\mathrm{map}}(e,z;T_{\mathrm{ev}},T_{\mathrm{op}})
=
\frac{
\left\|H\mu_{\Theta}(T_{\mathrm{ev}};e,z,T_{\mathrm{op}})
-y^{\mathrm{EMT}}_{e,z}(T_{\mathrm{ev}};T_{\mathrm{op}})\right\|_2
}{
\left\|y^{\mathrm{EMT}}_{e,z}(T_{\mathrm{ev}};T_{\mathrm{op}})\right\|_2
}.
```

Here `H` denotes the common observation convention, `mu` the surrogate mean,
and `y_EMT` the corresponding projected EMT state. This is a finite-time state
comparison, not an error in overload probability or the pathwise discrepancy of
one stochastic realization.

The reported ranges belong to different experiments:

| Result | Timing and operating conditions | Reported relative state error |
|---|---|---|
| Figure 4, global stochastic model `M2` | `T_op = 1.5 s`, `T_ev = 2.5 s`; all 46 injection scales | Approximately `0.0968–0.15` |
| Figure 9 at reclosure | `T_op = T_ev = 1.5 s`; nominal injections | Approximately `0.068–0.086` |
| Figure 9 after reclosure | `T_op = 1.5 s`, `T_ev = 2.5 s`; nominal injections | Approximately `0.11–0.15` |
| Figure 10 at its longest tested duration | `T_op = T_ev = 1.5 s`; nominal injections | Approximately `0.056–0.0815` |

Figures 9–10 contain **no injection-rescaling sweep**. Figure 9 changes the time
at which a fixed switching experiment is observed; Figure 10 changes the
switching time while keeping the observation time fixed. Their reported curves
were generated using different random-seed configurations and are retained as
separate numerical records. The endpoint values are not overwritten to make
the two studies coincide, and the nominal-injection range in Figure 9 is not
substituted for the scaled-injection range in Figure 4.

**Mean versus sampled evaluation:** the exact mean obtained from the
deterministic moment equations is seed-independent. For any figure generated
from seeded surrogate trajectories, its run record must identify the numerical
quantity used: a moment solution, an empirical ensemble mean, or an individual
realization, together with the seeds and any ensemble size. Different seeds
identify different sampled evaluations; they do not change the definition of
`RSE_map`. The existing figures and reported values are retained without
assuming an undocumented ensemble construction.

### Scope of the screening calculation

The 25 training lines and 25 testing lines do not define the full screening
catalogue. Paper 3 samples among **329 initiating lines**, using the same fixed
global parameter pair. This extends the model's application beyond the
line-disjoint assessment set; it is not an EMT validation of all 329 lines.
The screening horizon is `T_h = 10 s`, distinct from the `0.8 s` calibration
time and `2.5 s` hierarchy-assessment time.

The nominal open-phase duration has exponential **mean `0.1 s`**. This
probability law is separate from the fixed durations used for calibration and
stress testing. The conditional risk matrix averages over the nominal duration
and surrogate drivers; the duration-resolved plots condition on duration
instead. Agreement in the state metric does not establish accuracy of every
resulting overload-probability estimate.

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
│       ├── emt/                    # ParaEMT1.zip through ParaEMT14.zip
│       └── surrogate/              # Surrogate1.zip through Surrogate14.zip
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
and the supplied effective inertia and damping values. These processed-file
counts do not establish a one-to-one representation of the manuscript's
329-line screening catalogue. Resolve that mapping before reproducing the
screening experiment; do not replace its catalogue by the training/testing
lists or assume that every source line is a row in the processed branch table.

Original bus IDs and state indices are distinct. Use `bus_map.csv` and
`branch_map.csv` when matching trajectories, fault labels, and original case
files. Input schemas are documented in [DATA_DICTIONARY.md](docs/DATA_DICTIONARY.md).

## Piecewise stochastic surrogate

For a faulted line `e = (i,j)`, let its endpoint set be
$`\partial e=\{i,j\}`$. With pole opening at time zero, the model is

$$
\begin{cases}
dX_t = \bigl[A_e(\alpha)X_t+b(z)\bigr]\,dt
+ \displaystyle\sum_{r\in\partial e}G_{e,r}(\sigma_{e,z,r})X_t\circ dW_t^r,
& 0\leq t<T_{\mathrm{op}}, \\
\dot{X}_t = AX_t+b(z),
& T_{\mathrm{op}}\leq t\leq T_{\mathrm{end}}.
\end{cases}
$$

`T_end` denotes the final time of the particular run: `T_ev` for a finite-time
state comparison or `T_h` for overload screening. It does not identify those
two quantities as the same experimental setting.

The state is continuous at reclosure. The independent Wiener drivers act only
during the open-phase interval. The affected line has effective coupling
$`\alpha\beta_e`$ during that interval and nominal coupling $`\beta_e`$ afterward.

The code's time arguments are:

| Argument | Meaning |
|---|---|
| `T1` | Pole-opening time; `0.0 s` for the manuscript protocols |
| `T2` | Reclosure time; equals `T_op` when `T1 = 0` |
| `T2 - T1` | Open-phase duration |
| `T3` | `T_ev` for state comparisons; `T_h` for screening |
| `saveat` | Spacing of saved trajectory observations |
| `dt` | Maximum internal integration step during stochastic propagation |

For calibration-reference evaluation, use `T1=0.0`, `T2=0.8`, and `T3=0.8`.
For the Figure 4 assessment, use `T1=0.0`, `T2=1.5`, and `T3=2.5`. The state
is not reset at `T2`; when `T2=T3`, the run ends at the reclosure state.

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

using .PowerGridsFunctions3

case = load_case(
    "data/wecc240/processed/branch240E3.csv",
    "data/wecc240/processed/bus240E3.csv"
)

system = prepare_system(case.df, case.df2; injection_scale=1.02)
line_row = 1  # API example only; resolve a manuscript fault using branch_map.csv.
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
    T2=1.5,
    T3=2.5,
    saveat=0.01,
    dt=0.001,
    backend=:stratonovich_heun,
    rng=Xoshiro(2026),
    include_frequency=true
)
```

This single-trajectory example uses the Figure 4 time settings and one of its
injection scales. The example row, seed, internal step, and saved grid are not
a reconstruction of the figure's complete run configuration. To evaluate a
calibration-reference trajectory, set `T2=T3=0.8` and select a training line;
parameter fitting itself uses the deterministic moment equations and the
calibration objective, not this single realization. For Figures 9–10, set
`injection_scale=1.0` and vary only the time specified in the protocol table.
Match each reported parameter set to its actual model definition and numerical
records; changing the model or integration method does not establish that an
earlier calibration remains valid.

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
S_{ij}=\int_0^{T_{\mathrm{h}}}
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

The manuscript uses unweighted elite-based adaptation: 250 samples per
iteration, 20 iterations, 25 elites, and retention factor `kappa=0.9`. It learns
categorical fault probabilities `phi` and one shared exponential duration
**mean** `r`. This factorized proposal does not learn line–duration dependence
and is not a likelihood-ratio-weighted CE adaptation. The initial proposal mean
is `r = 1.0 s`; the nominal mean remains **`lambda0 = 0.1 s`** throughout.

For a catalogue of `E` initiating lines, uniform nominal fault selection and
the shared-mean proposal give

```math
p_Z(e,\tau)=\frac{1}{E\lambda_0}
\exp\!\left(-\frac{\tau}{\lambda_0}\right),\qquad
q_\nu(e,\tau)=\frac{\phi_e}{r}
\exp\!\left(-\frac{\tau}{r}\right).
```

Here `E = 329` for the manuscript screening catalogue; it is not `25`, the
number of training or testing lines. Both `lambda0` and `r` are in seconds.
Their reciprocal rates are `1/lambda0` and `1/r`, respectively.

After the 20 adaptation iterations, freeze the learned line probabilities and
shared duration mean. For final-batch sample `k`, the complete likelihood ratio
is

```math
w_k=\frac{r^{(20)}}{E\phi_{e_k}^{(20)}\lambda_0}
\exp\!\left[
\left(\frac{1}{r^{(20)}}-\frac{1}{\lambda_0}\right)\tau_k
\right].
```

The final estimation law retains **both learned components**. The stored fields
`res.π` and `res.r` must be interpreted according to the generating code's
parameter convention, rather than inferred from their names.

### Mean–rate convention in the sampler

With a mean parameter, Julia's unit-exponential draw is scaled by multiplication:

```julia
using Random
rng = Xoshiro(2026)  # Example stream, not a manuscript run identifier.

lambda0_mean_s = 0.1
proposal_mean_s = 1.0  # Initial value; freeze the learned value for final sampling.

tau_nominal = lambda0_mean_s * randexp(rng)
tau_proposal = proposal_mean_s * randexp(rng)
```

The earlier README described a sampler using `randexp(rng)/r`, while writing
the densities in the mean convention above. That expression treats the code
variable `r` as a **rate**, not a mean. If a code revision retains a rate-based
interface, use explicitly named reciprocal rates and a likelihood consistent
with the distribution actually sampled. For rates `rate0 = 1/lambda0` and
`rateq = 1/r`, the same likelihood is

```math
w_k=\frac{\mathrm{rate0}}{E\phi_{e_k}\,\mathrm{rateq}}
\exp\!\left[(\mathrm{rateq}-\mathrm{rate0})\tau_k\right].
```

This README does not change the Julia sampler. Check the sampler, adaptation
update, configuration, likelihood, and exported metadata together before a new
run. In particular, smoothing means and smoothing rates are different updates;
renaming an existing rate field does not implement the manuscript's mean-based
adaptation. Preserve the convention used in historical outputs instead of
silently relabeling them.

### Probability estimates

Use a fresh final batch with the same conditional stochastic-driver law as the
nominal model, so no additional driver likelihood factor is required. The
assessment is conditioned on a fault; it does not sample a no-fault category.
For a strict global exceedance, the estimate is

```math
\widehat Q(\gamma)=\frac{1}{N}\sum_{k=1}^{N}
w_k\,\mathbb{I}\{S_k>\gamma\}.
```

For the risk-matrix entry conditional on initiating line `e`, use

```math
\widehat Q_{ij}(e;\gamma)=\frac{1}{N\pi_e}
\sum_{k=1}^{N}w_k\,
\mathbb{I}\{e_k=e\}\,
\mathbb{I}\{S_{ij}^{(k)}>\gamma\},
\qquad \pi_e=\frac{1}{E}.
```

The nominal mass `pi_e` is distinct from the learned proposal mass `phi_e`.
The displayed risk matrix uses `gamma = 1 s` and `T_h = 10 s`. Its entries
measure more than one second of accumulated exposure within that horizon,
not the probability that the open-phase duration itself exceeds one second.

All surrogate randomness must use the RNG passed to the scorer. Failed or
nonfinite solves raise an error rather than contributing a safe-event indicator.
Final estimates use ordinary importance weights, not self-normalized weights,
and are not clipped to `[0,1]`. All event comparisons are strict `>`, including
at `gamma = 0`. Zero observed events do not establish zero risk.

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
| `adaptation_history.csv` | Iteration-level scores, thresholds, and duration-parameter history, with the mean/rate convention identified |
| `adaptation_phi.csv` | Categorical proposal history |
| `proposal.csv` | Frozen line probabilities and shared duration parameter; identify whether the stored field is a mean or rate |
| `run.toml` | Experiment settings, seeds, environment, timing scope, and source/input hashes |
| `bus_map.csv`, `branch_map.csv` | State-index and original-case mappings |
| `checksums.csv` | Checksums of saved run files |
| `Project.toml`, `Manifest.toml` | Project specification and resolved environment when available |

Preserve samples as well as seeds. Record the phase reference, mass normalization,
state ordering, endpoint amplitudes, opening and reclosing times, internal `dt`,
output `saveat`, and event-comparison convention. Keep the exact source revision
that generated each result.

For calibration and state-error figures, also record the experiment name,
training/testing membership, injection scale, `T_op`, `T_ev`, observation mask,
model level, and parameter-fit identifier. Distinguish moment-based state
predictions from single trajectories and empirical ensemble means. Include any
trajectory seeds and ensemble sizes without treating a seed as an input to the
exact deterministic moment solution.

Record adaptation and estimation runtimes separately, including compilation and
initialization policies. The current run driver times computation separately
from file serialization. Read [REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) for the
run-record specification and [VALIDATION_STATUS.md](docs/VALIDATION_STATUS.md) for
the checks completed during preparation. Those records must be updated when the
phase-endpoint implementation is tested.

## Calibration archive and ParaEMT provenance

The calibration-reference collection contains **25 training lines × 46
operating scales = 1,150 scenarios**, with `T_op = T_ev = 0.8 s`. The Figure 4
assessment uses **25 disjoint testing lines × the same 46 scales**, with
`T_op = 1.5 s` and `T_ev = 2.5 s`. The exact partition is given under
[Training and testing fault identifiers](#training-and-testing-fault-identifiers).
Archive membership and completeness must be established from scenario metadata,
not from the number of ZIP files.

The previously documented sweep over 15 open-phase durations contains **17,250**
line/scale/duration combinations. It is a separate archive description, not the
definition of Figures 9–10: those figures use nominal injections without an
injection-rescaling sweep. None of these counts determines the number of
independent stochastic realizations.

The calibration dataframes are stored as ZIP archives because of the number
and total size of the CSV files:

| Data source | Repository directory | Available archives |
|---|---|---|
| ParaEMT reference trajectories | [`data/wecc240/calibration/emt/`](data/wecc240/calibration/emt/) | `ParaEMT1.zip` through `ParaEMT14.zip` (14 archives) |
| Surrogate trajectories | [`data/wecc240/calibration/surrogate/`](data/wecc240/calibration/surrogate/) | `Surrogate1.zip` through `Surrogate14.zip` (14 archives) |

Download the required archives from the appropriate directory and extract
each into its own subdirectory, named after the archive, to avoid overwriting
files with matching names. For example, extract `ParaEMT1.zip` into
`data/wecc240/calibration/emt/ParaEMT1/` and `Surrogate1.zip` into
`data/wecc240/calibration/surrogate/Surrogate1/`. Match reference and surrogate
trajectories using their scenario identifiers and metadata, not ZIP numbers.

Keep EMT references and surrogate outputs in different directories. A historical
filename beginning with `ParaEMTDF` does not identify the generating simulator.
Each manifest entry must record that source, the exact line mapping, operating
scale, `T_op`, `T_ev` or `T_h`, projection/units, experiment and train/test
assignment, and checksum. Include the random-seed configuration when the
trajectory generation or numerical assessment uses sampling. Two files with
the same fault label need not belong to the same experiment.

After extraction, populate the calibration metadata with the **training**
trajectory records and paths to the extracted CSV files, then run:

```sh
python scripts/index_calibration.py --metadata YOUR_METADATA.csv --data-root . --output data/wecc240/metadata/calibration_manifest.csv --require-1150
```

The `--require-1150` count applies to the intended 1,150-scenario training
collection. Do not pool training, testing, and multi-duration records merely to
meet that count. Keep separate experiment manifests and verify the actual line
identifiers, the 46 scales, and the recorded times against this README; a count
check alone does not verify the partition.

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

The ParaEMT reference and surrogate calibration dataframes are publicly
accessible in the ZIP archives linked under
[Calibration archive and ParaEMT provenance](#calibration-archive-and-paraemt-provenance).

Freeze the tested code as a versioned release and archive the associated reference
data with a manifest and checksums. Link the code release and data archive in both
directions. Fill the release identifiers in `manuscript/references.bib` and use
the after-release Data Availability statement only after the listed materials
are publicly accessible.

The manuscript cites this repository for the exact fault identifiers and the
training/testing partition documented above. A repository citation is:

> A. Almada, “PowerGridREserach: Power-grid dynamic screening code,” GitHub
> repository. https://github.com/AyrtonAlmada/PowerGridREserach

Associate that citation with the actual code commit and data manifest when
reproducing a reported experiment. A link to the current repository does not
identify the historical revision that generated a result.

Related work and reference software are listed in the manuscript bibliography,
including the N1Plus studies and the ParaEMT publications. Preserve upstream
notices for the original case files and distinguish third-party licensing from
the license governing this project's own code. The repository's license template
must be completed by the rights holders before release.

This software supports research into contingency screening. It is not a protection
relay implementation or an operational security certification. Probability-zone
boundaries and overload-duration thresholds must be identified for each study;
they are not universal protection requirements.
