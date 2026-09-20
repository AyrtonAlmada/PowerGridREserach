# Reproducibility record

## What a reader needs

A paper release should contain the exact input files, the source/version that ran,
a resolved environment, numerical settings, random streams, generated data, and
analysis scripts that connect those data to the paper's figures and tables.
The manifest should distinguish raw EMT, projected EMT, surrogate trajectories,
and aggregated probability outputs. File names alone do not establish provenance.

## Three records, not one moving branch

1. **Code release:** Git tag + complete commit SHA + source checksum; archived DOI.
2. **ParaEMT environment:** actual checkout SHA, local tracked patch, untracked
   experiment scripts, initialization snapshot hash, Python dependency lock and
   numerical/device settings. Capture from the machine/environment that ran it.
3. **Data release:** immutable trajectories with checksums, scenario metadata, units,
   projection/filtering code and train/test membership; separate dataset DOI if large.

`provenance/input_checksums.csv` records byte identity for the supplied inputs.
`provenance/paraemt_provenance.template.json` has blank historical fields because
those values cannot be inferred from the uploads. `main` is not an immutable version.

## Determinism and samples

The supported sampler accepts an explicit RNG. Every score call gets its own driver
seed, stored as hexadecimal text. Adaptation and estimation use distinct streams.
Record the Julia version and generated `Manifest.toml`; do not promise bitwise
agreement across Julia/BLAS/platform versions. Retaining final samples/logweights
allows reanalysis without relying exclusively on seed replay.

The fresh final batch must use `res.π` and `res.r` in BOTH sampling and likelihood
weighting. All exponential parameters are rates. Record any deviation from the
uniform nominal fault law, the conditional no-fault convention, and the treatment
of durations longer than the observation window. The provided implementation
samples faults only; it does not infer or estimate fault-occurrence frequencies.

## Timing

`run_cem` measures adaptation and estimation/diagnostic calculations separately.
CSV serialization is outside those timers. The driver additionally records input
loading and nominal-model preparation time. Lazy per-line model preparation is
included at its first use during a timed stage. No warmup is performed implicitly;
record compilation/caching policy and use the same scope when comparing methods.

The summary may contain several gamma rows derived from a single final batch.
They share one cost and one proposal; do not add that adaptation time once per row.
Different thresholds have different event-specific uncertainties, but all-weight
ESS and maximum normalized weight are identical when the final batch is shared.

## Diagnostics

Report the exact final batch size, event count, probability estimate, absolute and
relative SE, weight ESS and maximum normalized weight. Retain the samples behind
those values. For zero hits, the code records SE as missing rather than using a zero
empirical variance to claim certainty. Ordinary finite-sample IS estimates may
exceed one; they are not clipped or silently self-normalized.

For exponential nominal/proposal rates lambda/r, r < 2lambda is sufficient for
finite second moments of all weights when all categorical probabilities are
positive. An observed finite SE is not proof of that property. A high weight ESS is
not itself evidence of accurate rare-event probabilities.

## Calibration archive

`index_calibration.py` accepts author-supplied metadata, checks real files, hashes
them and, with `--require-1150`, verifies 25 lines x 46 lambda values at one pair of
switching/evaluation times. It does not invent the selected line list, a train/test
split, the generating solver, or a historical commit.

Archive the actual EMT reference trajectories, not only surrogate replicas. For
large CSV collections, a compressed archive with an index and checksums is adequate;
keep each file's source designation and measurement/processing history. Include
both raw and projected outputs when projection/filtering changes interpretation.
A separate 15-duration surrogate sweep is not the original 1,150-reference set.
When a historical Git commit genuinely cannot be recovered, use the explicit
unknown-commit workflow in PARAEMT_PROVENANCE.md and archive the actual source
bundle. Do not substitute an unrelated current commit.

## Release gates

- Run the Julia tests and save their output; freeze the actual environment.
- Resolve the noise/backend and processed-topology differences in MODEL_AND_IMPLEMENTATION.md.
- Add the omitted raw-to-processed conversion, calibration fitting and plot scripts.
- Capture the actual ParaEMT revision and experiment-specific local modifications.
- Deposit and index the 1,150 EMT references; check scenario/channel metadata.
- Verify generated figures/tables against stored final data; do not substitute
  illustrative formatting examples for measured outputs.
- Confirm data redistribution permissions and licenses.
- Fill release/DOI fields in citations; only then use the public-release DAS.

APS policy: https://journals.aps.org/authors/data-availability-statements
Julia RNG guidance: https://docs.julialang.org/en/v1/stdlib/Random/#Reproducibility
