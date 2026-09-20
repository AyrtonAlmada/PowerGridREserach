# EMT reference trajectories — deposition pending

The author reports that the 1,150 calibration trajectories are ready for upload.
They were not included in the attachments used to prepare this release candidate.
No trajectory has been created or labeled as an EMT result by this package.

Deposit the actual 25-line x 46-operating-state set here, or in a DOI-bearing data
archive. Maintain an explicit manifest with source, line IDs, lambda, opening and
reclosure times, evaluation horizon, channel order/units, integration step, output
sampling step, version, and checksum. `scripts/index_calibration.py` verifies the
base Cartesian grid after metadata is populated. The actual 25 line rows must be
provided by the author; they are not inferred from figure labels or filenames.

The 15-duration surrogate sweep contains 17,250 trajectories, not 1,150 EMT
references. Store that sweep separately under `../surrogate/`.
