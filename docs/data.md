# Nuclear data

nuclide bundles a fixed **NUBASE2020** snapshot. Queries work offline and return
the same values across installations of a package version. The snapshot was
downloaded on 2026-11-16; this date does not make it a 2026 evaluation. NUBASE2020
considers experimental literature available through 2020-10-30.

The source is the [AMDC machine-readable NUBASE2020 table](https://nds.iaea.org/amdc/ame2020/nubase_4.mas20.txt).
The evaluation, its policies, and its state conventions are described in
[Kondev et al., Chinese Physics C 45, 030001 (2021)](https://doi.org/10.1088/1674-1137/abddae).
The [data attribution notice](../src/nuclide/data/NOTICE.md) contains the license,
source hash, and transformation summary. Code and data have separate licenses:
the Python implementation is MIT, and the NUBASE2020 content is CC BY 3.0.

## Coverage

The source contains 5,843 rows. This release includes **5,375 states**:

| State | Records |
| --- | ---: |
| Ground state | 3,555 |
| First isomer, `m` | 1,366 |
| Second isomer, `n` | 454 |

These are records in the evaluation, including estimated properties; the count
does not mean that every included nuclide has been experimentally observed.
The neutron is included as `n-1`. Higher isomers, other excited levels,
resonances, and isobaric analogue states are outside this release's scope.
The source uses some of the same labels for higher isomers and other levels,
so they are not automatically mapped to extra isomer numbers.

Among the ground/`m`/`n` candidates, 11 entries explicitly marked `RN`
(reported nonexistent) and 13 marked `EU` (existence uncertain) are omitted.
All 468 excluded rows, including unsupported states, are listed by name and
reason in the generated metadata. The unmodified compressed source retains
their complete information.

## Half-lives and qualifiers

`half_life_status` distinguishes the meaning of a value:

| Status | Meaning | Records |
| --- | --- | ---: |
| `measured` | Evaluated central value without a special qualifier | 4,507 |
| `estimated` | Value marked `#`, estimated from systematic trends | 413 |
| `approximate` | Value marked `~` | 28 |
| `stable` | Source reports `stbl` | 254 |
| `lower_limit` | Only a lower limit is given | 58 |
| `upper_limit` | Only an upper limit is given | 52 |
| `unknown` | No numerical half-life or recognized limit is given | 60 |
| `unstable` | Particle-unstable, with no numerical half-life | 3 |

Here `measured` describes the absence of a source qualifier; it does not imply
an exact value. Evaluated central values still have measurement uncertainties.
In NUBASE, `stbl` also covers a nuclide for which no finite half-life has been
established. It is not a proof of infinite physical lifetime.

For example, Rn-222 has `3.8215 d`, Tc-99m has `6.0066 h`, and B-18 has only
`<26 ns`. The importer preserves the last value as an upper limit, including
when the source places it in the uncertainty column. A bound is never promoted
to a measured half-life. Numerical decay calculations require a central value
or a source-stable state; missing values and bounds cannot determine a unique
decay constant.

The evaluation uses a **tropical year of 31,556,926 seconds**, as specified in
section 2.5 of the paper. This differs from a Julian year of 365.25 days.
The imported second values use the NUBASE convention. Converting a resulting
quantity into another library's year unit uses that unit's own definition.

`half_life_uncertainty_raw` preserves the source's uncertainty field, including
asymmetric errors and supplementary limits. `raw_line` preserves every original
fixed-width row. The current numerical API uses central values and does not
propagate uncertainties. Mass-excess estimates retain their `#` flag via
`mass_excess_estimated`; mass-excess uncertainty and excitation energy are also
stored in the data.

## Decay modes are descriptive data

The source's decay-mode strings are preserved. They are **not a complete set
of mutually exclusive transitions between specific nuclear states**.
For example, `B-=100;B-n=16` includes a delayed-neutron channel within the beta
decay total. Adding these values would double-count part of the decay.
An internal transition can also feed another isomer. Therefore this release
does not infer daughter nuclides, complete decay chains, or daughter buildup
from these strings. Parent survival calculations describe disappearance of
the initial nuclide only.

## Reproducing the bundle

The original text is included as a deterministic gzip file. From the repository
root, run:

```sh
python tools/import_nubase.py
```

Or provide the uncompressed snapshot explicitly:

```sh
python tools/import_nubase.py /path/to/nubase_4.mas20.txt
```

The importer checks the source SHA-256 before conversion and generates
`src/nuclide/data/nubase2020.json` deterministically. It rejects unreviewed source
changes rather than silently replacing the data. Importing a different
evaluation requires reviewing its format and scientific changes, updating the
importer and provenance, and rerunning the data tests.
