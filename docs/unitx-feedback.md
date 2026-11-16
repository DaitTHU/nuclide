# unitx feedback from nuclear calculations

The initial nuclide implementation uses the published `unitx` 0.17.13 API.
Scalar calculations work without NumPy. The observations below come from
running a consumer against the installed wheel and checking its annotations
with Pyright; they distinguish useful existing behavior from possible
improvements. No changes to unitx are required for these calculations.

## Quantities compose naturally

Multiplication by a unit, ordinary arithmetic, and `.to()` cover the basic
nuclear calculations. `si` contains units; `SI` contains physical constants.
For example, with an illustrative half-life and atomic mass:

```python
from math import log

from unitx import math as umath
from unitx.si import SI, si

half_life = 4.0 * si.d
atomic_mass = 222.0 * si.Da
decay_constant = log(2) / half_life
specific_activity = (decay_constant / atomic_mass).to(si.Bq / si.g)
molar_mass = (atomic_mass * SI.NA).to(si.g / si.mol)

atoms = ((1 * si.ug) / atomic_mass).n
activity = (atoms * decay_constant).to(si.Bq)
mass = (activity / decay_constant * atomic_mass).to(si.ug)
remaining_fraction = umath.exp(-decay_constant * half_life).n
assert remaining_fraction == 0.5
assert mass.isclose(1 * si.ug)
```

The numerical values here demonstrate the interface; they are not nuclear
reference data. In nuclide, the properties come from the bundled dataset.

`.n` explicitly extracts a coherent dimensionless number. It fits both a
number of atoms and a survival fraction. `.value_in(unit)` is the corresponding
boundary when passing a physical quantity to a numerical routine.

The example retains `Quantity[float]` for physical results and `float` for the
dimensionless results, without casts or forced numerical conversions.
A parameter annotated as `Quantity[float]` also accepts the natural integer
construction `1 * si.d`.

## Inverse-unit conversion is easy to mistype

The natural-looking expression below raises `TypeError`:

```python
decay_constant.to(1 / si.s)
```

`1 / si.s` constructs a quantity under unitx's number–unit arithmetic rules,
whereas `.to()` takes a unit or unit spelling. The supported forms are:

```python
decay_constant.to(si.s**-1)
decay_constant.to("1/s")
activity.to(si.Bq)
```

This behavior is consistent with the existing object model. A targeted error
message suggesting `si.s**-1` or `'1/s'` would make this common mistake easier
to correct.

## Curie is missing from the standard catalogue

`si.Bq` and its common SI prefixes are available. `si.Ci` and `si.uCi` raise
`AttributeError` in 0.17.13. Curie and its prefixes would be useful additions
for nuclear calculations. nuclide uses becquerels and does not modify unitx's
global registry on import.

## A year needs an explicit convention

In unitx, `si.yr` is exactly 365 days, while `si.a_j` is the Julian year of
365.25 days. A nuclear dataset's stated year convention must determine the
conversion; its `y` label cannot automatically be mapped to `si.yr`.
Both units are documented, but the difference is easy to miss in data imports.

NUBASE2020 specifically adopts a tropical year, giving `1 y = 31,556,926 s`.
Neither `si.yr` nor `si.a_j` matches that convention. nuclide converts the
dataset's year-based values to seconds using its stated factor. See
[NUBASE2020, section 2.3](https://cpc.ihep.ac.cn/en/article/doi/10.1088/1674-1137/abddae).

## Atomic mass and constants already fit

`si.u` and `si.Da` refer to the same atomic mass unit and display as `Da`.
`SI.c` is the speed-of-light quantity; `SI.NA` is the Avogadro constant.
An atomic mass excess can therefore be used directly:

```python
atomic_mass = (mass_number * si.u + mass_excess / SI.c**2).to(si.u)
```

Multiplying an atomic mass by `SI.NA` derives its molar mass while retaining
the current constants. These are atomic masses, including the atomic
electrons; a nuclear mass requires a separate physical model.
