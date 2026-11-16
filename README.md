# nuclide

Nuclear data and radioactive decay calculations that read like physics.

```python
from nuclide import Nuclide
from unitx.si import si

rn = Nuclide("Rn222")
rn.Z, rn.N, rn.A  # (86, 136, 222)
rn.half_life.to(si.d)  # 3.8215 d (NUBASE2020)
rn.survival_fraction(rn.half_life)  # 0.5
rn.activity(mass=1 * si.ug).to(si.Bq)  # activity of pure Rn-222
```

## Install

```sh
pip install nuclide
```

Python 3.12 or newer. The only runtime dependency is
[unitx](https://github.com/DaitTHU/unitx). Nuclear data are bundled: queries and
calculations work offline, with no NumPy or SciPy required.

For a local checkout before publication: `python -m pip install -e .`.

## A nuclide is an identity

```python
assert Nuclide("Rn222") == Nuclide("Rn-222") == Nuclide("222Rn")
assert Nuclide("Tc99m") != Nuclide("Tc99")

rn = Nuclide("Rn222")
print(rn)  # Rn-222
repr(rn)  # "Nuclide('Rn-222')"
rn.element  # 'Rn'
rn.state  # '' (ground state); 'm' or 'n' for the first two isomers
```

Objects are immutable and hashable, so they work as dictionary keys. The spellings
`Tc-99m`, `Tc99m1`, and `99mTc` identify the same first isomer; `m2` is an alias for
the source's second-isomer suffix `n`. Element symbols are case insensitive.
`Nuclide("n1")` denotes the free neutron.

```python
from nuclide import nuclides

radon_states = nuclides("Rn")
stable_carbon = [n for n in nuclides("C") if n.is_stable]
```

Unknown names and states outside the bundled evaluation raise `ValueError`.

## Physical properties carry units

```python
rn.half_life.to(si.h)
rn.decay_constant.to("1/s")
rn.mean_lifetime.to(si.d)
rn.atomic_mass.to(si.u)
rn.mass_excess.to(si.keV)
rn.molar_mass.to(si.g / si.mol)
rn.specific_activity.to(si.Bq / si.g)

tc = Nuclide("Tc99m")
tc.excitation_energy.to(si.keV)
tc.spin_parity
tc.decay_modes  # source descriptions, in percent
```

Use `.to(unit)` to keep a quantity, or `.value_in(unit)` to get a scalar for
numerical code. `atomic_mass` is the **neutral atomic mass**, including electrons,
computed as `A * u + mass_excess / c**2`. It is not the bare nuclear mass.

## Count, mass, activity, and time

Each conversion takes exactly one named input. A nuclide carries no mutable
sample amount; you specify the amount for each calculation.

```python
rn = Nuclide("Rn222")

count = rn.atoms(mass=1 * si.ug)
activity = rn.activity(atoms=count)
mass = rn.mass(activity=activity).to(si.ug)

remaining = count * rn.survival_fraction(7 * si.d)
remaining_activity = rn.activity(atoms=remaining)
```

These are **single-parent, pure-isotope** calculations. Survival follows
`N(t) / N(0) = exp(-ln(2) * t / half_life)`, with no production, transport, or
daughter growth. Inputs must be finite, nonnegative scalars; times, masses and
activities must carry compatible units. Counts may be nonintegral expected
populations.

For stable nuclides the half-life is infinite, the decay constant and activity
are zero, and the survival fraction is one. Activity cannot determine a stable
nuclide's mass or atom count, so those inverse operations raise `ValueError`.

## Know what the evaluation says

The data are a fixed, attributed **NUBASE2020** snapshot, including ground states
and the first two isomers. Evaluated values, estimates, bounds, and missing values
remain distinct:

```python
rn.half_life_status  # 'measured'
rn.half_life_text  # source value and unit
rn.half_life_uncertainty_text  # original uncertainty notation
rn.mass_is_estimated  # source '#' flag
rn.mass_excess_uncertainty  # tabulated energy uncertainty
```

Estimated or approximate central values can be used in calculations; check
`half_life_status` and `mass_is_estimated` when that distinction matters. A
half-life given only as a limit, unknown, or particle-unstable without a lifetime
raises `DataUnavailableError` on numerical access. Uncertainties are not
propagated. Ground/isomer ordering qualifications remain in the bundled source
records; the package follows NUBASE labels.

`decay_modes` preserves source descriptions such as `A=100`. These descriptions
may overlap: beta-delayed particle emission can be part of a reported total beta
rate. They are **not** mutually exclusive, daughter-state-resolved branching
fractions. Version 0.1 does not compute full decay chains, daughter populations,
emission spectra, or dose/shielding quantities.

See [data provenance and conversion rules](docs/data.md). In particular, NUBASE
uses a tropical year of 31,556,926 seconds; input `si.yr` follows unitx's 365-day
definition. The package converts source lifetimes to seconds before calculation.

## Origin and license

This is a new Python implementation inspired by
[DaitTHU/Radiation](https://github.com/DaitTHU/Radiation), an undergraduate MATLAB
project. It keeps the idea of a natural nuclear-physics interface while using
evaluated data with explicit nuclear-state identities.

Python code: [MIT](LICENSE). NUBASE2020 data: separately attributed under
[CC BY 3.0](https://creativecommons.org/licenses/by/3.0/), with the full source and
transformation notice in [NOTICE.md](src/nuclide/data/NOTICE.md).

F. G. Kondev, M. Wang, W. J. Huang, S. Naimi and G. Audi,
*The NUBASE2020 evaluation of nuclear physics properties*, Chinese Physics C
**45**, 030001 (2021), [doi:10.1088/1674-1137/abddae](https://doi.org/10.1088/1674-1137/abddae).

Development commands and publishing setup are in [CONTRIBUTING.md](CONTRIBUTING.md).
Notes from using unitx for these workflows are in [unitx feedback](docs/unitx-feedback.md).
