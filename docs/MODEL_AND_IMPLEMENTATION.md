# Source audit and numerical changes

This document separates facts found in the supplied files from the new release
code. It is intentionally not a claim that the paper's results were reproduced.

## Uploaded source facts

`PowerGridsFunctions3(2).jl` uses `XiMass` with state order `[omega; theta]`,
`A = [-M^-1 D  -M^-1 L; I  0]`, and `b=[M^-1 P;0]`.
`Ysol1` computes the initial phase from `pinv(L_nominal)*P_nominal`.

Its noise matrix is in the **upper-right** block: it multiplies the endpoint phase
difference and enters the frequency equation. Its square is exactly zero. The
frequency-diagonal noise matrices discussed elsewhere are a different model.
The code cannot be cited as showing that a nonzero `G^2/2` renormalizes alpha.

The uploaded `make_f1` evaluates independent Gaussian scalars at every requested
time, with standard deviations `t+1` and `(t+1)^3/3`. It includes nonzero random
perturbations at time zero and does not construct Brownian increments shared across
times. This is not a pathwise solution of the stated multiplicative-noise SDE.
The new `:source_random_map` option preserves those random-evaluation formulas so
the discrepancy can be inspected. It is not described as an exact SDE solution.

`AnalyticalSolution` duplicated the switch time, redrew the state at that time, used
a global `n`, and replaced actual evaluation instants by a linspace. Several helper
functions duplicate definitions or depend on globals. The source overload routine
counts sample indicators with a common interval and converts to Float32.

`CEM.txt` uses a global scorer with physical `alpha=0`, a global `sigma`,
`Exponential(r)` sampling while treating `r` as a rate in the MLE and weights,
unweighted elite selection, one shared duration parameter, and no convergence
criterion beyond the iteration count. Its positional argument list can treat a
boolean as `pi_min`. Those are source observations, not measured consequences.

## Supported refactor

- Positive-inertia all-bus ODEs only; no silent mass-matrix pseudoinverse for loads.
- Component-balanced pre-fault equilibrium; no automatic slack redistribution.
- Original node/branch order and mappings preserved; no deduplication or topology repair.
- Deterministic affine propagation handles zero eigenvalues by value, not position.
- One stored reclosure state, strictly increasing actual output times, no global `n`.
- Float64 left-endpoint overload quadrature on `diff(Time)`; no added final interval.
- Explicit RNG and serial scorer callbacks; every simulation seed is stored.
- Rate-consistent exponential sampling; fresh final learned categorical proposal.
- Fixed-iteration unweighted adaptation remains unweighted, as in the source.
- Full final weights and diagnostics are stored; undefined zero-hit uncertainty is missing.
- An explicit new `:stratonovich_heun` backend provides pathwise SDE integration.
  It is an alternative for new runs, not an editorial correction to old outputs.
- Optional `:frequency_endpoints` noise is a second explicit model; it is not
  substituted for the uploaded upper-right-block noise.

The selected numerical changes can alter trajectories, overload scores, proposal
fits, probabilities, and runtime. Old fitted `(alpha,sigma)` values do not become
validated parameters of a new backend by being copied into a configuration file.

## API compatibility

The supported entry points include `Laplatian`, `XiMass`, `AnalyticalSolution`,
`OverheatingIndicator`, and `OverheatingIndicatorIndv`. The short root loader exports
these names into an interactive session. Nonzero stochastic runs must select a
backend. The plotting/alternative simulation routines from the 3,073-line original
remain under `legacy/source/PowerGridsFunctions3_original.jl`; they are not silently
loaded into the supported module or claimed to work with its reduced dependency set.
Use a separate historical environment when executing that original file.

## Checks still needed for the manuscript

1. Choose and document the actual noise matrix and trajectory generator used for
   each figure and calibration fit.
2. Re-run or revalidate fits after changing stochastic interpretation, backend,
   topology, or overload quadrature. Do not equate zero-mean added waveform noise
   with a change in the exact SDE mean.
3. Document the derivation of the processed network, including missing transformer
   connections, parallel-line handling, effective inertias, and injection construction.
4. Supply the calibration optimizer, empirical covariance construction, H projection,
   masks, operating-state splits, and actual EMT switching scripts. `moment_rhs`
   alone does not reproduce parameter identification.
5. Record numerical convergence and probability uncertainty. The default dt/saveat
   in an example is not a validated accuracy setting or a runtime benchmark.
