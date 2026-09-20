# WECC240 input and output dictionary

## Input byte identity

The seven user-supplied input/source files are preserved in the package. Original
upload names with `(1)`, `(2)`, and `(3)` are normalized only at the filesystem level;
the contents are unchanged. SHA-256 hashes are in `provenance/input_checksums.csv`.
The upstream filenames are not proof of the exact Git commit that supplied them.

## Upstream files

- `240bus.xls`: original dynamic-parameter workbook, retained in its original format.
- `WECC240.raw`: original power-flow case, including its research-only/header notices.
- `pfd_240_1_1.json`: saved ParaEMT power-flow data; this supplied file contains 243
  buses, 329 line entries, and 122 transformer entries. The system base field is 100 MVA.

The current loader consumes the processed CSVs, not these upstream formats. It does
not claim to reconstruct the exact conversion from the three upstream files. Add
the original conversion script and a record of parameter overrides to close that gap.

## Processed branches: `branch240E3.csv`

| Column | Loader interpretation |
|---|---|
| `index` | original branch identifier, preserved in `branch_map.csv` |
| `f_bus`, `t_bus` | external bus IDs; mapped using the order of the bus CSV |
| `br_x` | positive reactance-like input; coupling is computed as `1/br_x` |
| `rate_a` | nonnegative flow limit used directly as `F` |

The CSV does not encode units/base conversion. The loader uses the numbers as
supplied; the author must document the shared base for injections, coupling, and
limits. No voltage magnitudes, tap ratios, phase shifts, or transformer branches
are added. The 248 branch entries are not the complete 329-line/122-transformer
network recorded in the JSON. No parallel circuits are merged by this loader.

## Processed buses: `bus240E3.csv`

| Column | Loader interpretation |
|---|---|
| `index` | external bus ID; 243 unique IDs |
| `p` | nominal reduced-model active-power injection |
| `m` | effective positive inertia; all supplied values equal 1 |
| `d` | effective nonnegative damping; all supplied values equal 0.5 |
| `va` | retained as `Theta0` metadata, not used as the solved initial state |
| `longitude`, `latitude` | ignored by the dynamic solver; not certified GIS coordinates |

The initial phase is recomputed from `L*theta=P` with a pseudoinverse gauge. The
processed graph has 70 components; component balance is checked rather than merely
testing the total injection sum. Uniform operating-state scaling acts on `p` before
that equilibrium is recomputed, and the scaled injections are used in both stages.

## Trajectories and manifests

A surrogate CSV contains `x1`, ..., `xn`, optional `omega1`, ..., `omegan`, and
`Time`. State indices are distinct from external bus labels. Do not silently mix
radians, degrees, per-unit angular speed, and rad/s. The new state equation uses
`omega = dtheta/dt`; physical units/bases must be recorded when preparing the inputs.

Record the initial-state convention, alpha, sigma, noise structure, backend, time
step, output interval, T1, T2, T3, operating scale and seed for each trajectory.
The old `ParaEMTDFX...` prefix does not establish whether the generator was ParaEMT;
use the explicit `source` field. EMT output channel mapping requires its own H or
projection script; generic `x1` columns do not establish channel equivalence.

For the base 1,150-reference set, metadata must contain 25 selected fault rows and
all integer lambda codes 75 through 120 at one documented switching/evaluation
configuration. Empty templates are not uploaded data. A 15-duration sweep contains
17,250 cases and belongs in a separate collection.
