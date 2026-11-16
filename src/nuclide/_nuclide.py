"""Nuclear identities, evaluated properties, and single-parent calculations."""

import re
from dataclasses import dataclass
from math import inf, isfinite, log
from numbers import Real

from unitx import Quantity, Unit
from unitx import math as umath
from unitx.si import SI, si

from ._data import Record, records


class DataUnavailableError(ValueError):
    """The evaluation does not provide a usable value for this calculation."""


def _canonical_name(name: str) -> str:
    if not isinstance(name, str):
        raise TypeError("Nuclide requires a name such as 'Rn222' or 'Tc99m'.")
    text = name.strip()
    if match := re.fullmatch(r"([A-Za-z]{1,3})-?([1-9][0-9]*)(m[12]?|n)?", text):
        element, mass, state = match.groups()
    elif match := re.fullmatch(r"([1-9][0-9]*)(m[12]?|n)?-?([A-Za-z]{1,3})", text):
        mass, state, element = match.groups()
    else:
        raise ValueError(
            f"Invalid nuclide name {name!r}; use 'Rn222', 'Rn-222', or '222Rn'."
        )
    # The source uses m/n for the first/second isomer, respectively.
    state = {None: "", "m1": "m", "m2": "n"}.get(state, state)
    symbol = "n" if element.lower() == "n" and mass == "1" else element.capitalize()
    return f"{symbol}-{mass}{state}"


def _nonnegative(value: float, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{label} must be a real scalar.")
    if not isfinite(value) or value < 0:
        raise ValueError(f"{label} must be finite and nonnegative.")
    return value


def _quantity_value(value: Quantity[float], unit: Unit, label: str) -> float:
    if not isinstance(value, Quantity):
        raise TypeError(f"{label} must carry units; for example, 1 * si.{unit}.")
    return _nonnegative(value.value_in(unit), label)


@dataclass(frozen=True, slots=True, repr=False)
class Nuclide:
    """An immutable nuclear identity with NUBASE2020 evaluated properties.

    Names such as ``Rn222``, ``Rn-222`` and ``222Rn`` describe the same
    ground state. ``Tc99m`` is distinct from ``Tc99``. Numeric calculations
    use central values; source qualifiers remain available as properties.
    """

    name: str

    def __post_init__(self) -> None:
        name = _canonical_name(self.name)
        if name not in records():
            raise ValueError(
                f"Nuclide {self.name!r} is not in the bundled NUBASE2020 data. "
                "Supported states are ground states and the first two isomers (m/n)."
            )
        object.__setattr__(self, "name", name)

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return f"Nuclide({self.name!r})"

    @property
    def _record(self) -> Record:
        return records()[self.name]

    @property
    def Z(self) -> int:
        """Number of protons."""
        return self._record["Z"]

    @property
    def A(self) -> int:
        """Mass number (protons + neutrons)."""
        return self._record["A"]

    @property
    def N(self) -> int:
        """Number of neutrons."""
        return self.A - self.Z

    @property
    def element(self) -> str:
        """Element symbol, or ``n`` for the free neutron."""
        return self.name.split("-")[0]

    @property
    def state(self) -> str:
        """Empty string for the ground state; ``m`` or ``n`` for an isomer."""
        return self._record["state"]

    @property
    def is_stable(self) -> bool:
        """Whether NUBASE2020 explicitly classifies this state as stable."""
        return self.half_life_status == "stable"

    @property
    def half_life_status(self) -> str:
        """Source classification: measured, estimated, approximate, stable, etc."""
        return self._record["half_life_status"]

    @property
    def half_life_text(self) -> str:
        """Original half-life value and unit, preserving limits and qualifiers."""
        return self._record["half_life_raw"]

    @property
    def half_life_uncertainty_text(self) -> str:
        """Original uncertainty field; not a propagated numerical uncertainty."""
        return self._record["half_life_uncertainty_raw"]

    @property
    def half_life(self) -> Quantity[float]:
        """Half-life as a time quantity; infinity for an explicitly stable state.

        Limits, unknown half-lives and particle-unstable states without a
        lifetime raise DataUnavailableError. Estimated and approximate central
        values are usable, with their status retained in half_life_status.
        """
        if self.is_stable:
            return inf * si.s
        seconds = self._record["half_life_seconds"]
        if seconds is None or self.half_life_status not in {
            "measured",
            "estimated",
            "approximate",
        }:
            raise DataUnavailableError(
                f"{self.name} has no usable half-life: {self.half_life_text!r} "
                f"({self.half_life_status}). A limit is not an exact lifetime."
            )
        return seconds * si.s

    @property
    def decay_constant(self) -> Quantity[float]:
        """Parent loss rate ln(2) / half_life, in inverse seconds."""
        return (log(2) / self.half_life).to(si.s**-1)

    @property
    def mean_lifetime(self) -> Quantity[float]:
        """Mean lifetime half_life / ln(2)."""
        return self.half_life / log(2)

    @property
    def mass_excess(self) -> Quantity[float]:
        """Atomic mass excess (M - A u)c²; this can be negative."""
        value = self._record["mass_excess_kev"]
        if value is None:
            raise DataUnavailableError(f"{self.name} has no evaluated mass excess.")
        return value * si.keV

    @property
    def mass_excess_uncertainty(self) -> Quantity[float]:
        """Tabulated mass-excess uncertainty, without propagation."""
        value = self._record["mass_excess_uncertainty_kev"]
        if value is None:
            raise DataUnavailableError(f"{self.name} has no mass-excess uncertainty.")
        return value * si.keV

    @property
    def mass_is_estimated(self) -> bool:
        """Whether the mass was estimated from systematics (source '#')."""
        return self._record["mass_excess_estimated"]

    @property
    def atomic_mass(self) -> Quantity[float]:
        """Neutral atomic mass A u + mass_excess/c² (free mass for n-1)."""
        return (self.A * si.u + self.mass_excess / SI.c**2).to(si.u)

    @property
    def molar_mass(self) -> Quantity[float]:
        """Atomic mass times the Avogadro constant."""
        return (self.atomic_mass * SI.NA).to(si.g / si.mol)

    @property
    def excitation_energy(self) -> Quantity[float]:
        """Excitation energy of the nuclear state; zero for ground states."""
        if not self.state:
            return 0.0 * si.keV
        value = self._record["excitation_energy_kev"]
        if value is None:
            raise DataUnavailableError(f"{self.name} has no excitation energy value.")
        return value * si.keV

    @property
    def spin_parity(self) -> str:
        """Evaluated spin/parity text, retaining source qualifiers."""
        return self._record["spin_parity"]

    @property
    def decay_modes(self) -> tuple[str, ...]:
        """Source decay descriptions (percent), excluding isotopic abundance.

        These may overlap and are not daughter-state-resolved branching
        probabilities. For example, beta-delayed neutron decay is included
        in the total beta decay rate.
        """
        return tuple(
            part.strip()
            for part in self._record["decay_modes_raw"].split(";")
            if part.strip() and not part.strip().startswith("IS")
        )

    @property
    def specific_activity(self) -> Quantity[float]:
        """Activity per mass of this pure isotope, in Bq/g."""
        return (self.decay_constant / self.atomic_mass).to(si.Bq / si.g)

    def survival_fraction(self, time: Quantity[float]) -> float:
        """Fraction of the original parent nuclei remaining after time.

        Assumes an isolated population with no production or removal other
        than radioactive decay. Does not include daughter growth.
        """
        seconds = _quantity_value(time, si.s, "time")
        return umath.exp(-self.decay_constant * (seconds * si.s)).n

    def activity(
        self, *, atoms: float | None = None, mass: Quantity[float] | None = None
    ) -> Quantity[float]:
        """Activity of a pure sample, specifying either atom count or mass."""
        if (atoms is None) == (mass is None):
            raise TypeError("Specify exactly one of atoms= or mass=.")
        if mass is not None:
            count = self.atoms(mass=mass)
        else:
            assert atoms is not None
            count = _nonnegative(atoms, "atoms")
        return (count * self.decay_constant).to(si.Bq)

    def atoms(
        self,
        *,
        mass: Quantity[float] | None = None,
        activity: Quantity[float] | None = None,
    ) -> float:
        """Atom count from a pure sample's mass or activity.

        Counts may be fractional expected populations. Stable nuclides have
        zero activity at every mass, so activity cannot determine their count.
        """
        if (mass is None) == (activity is None):
            raise TypeError("Specify exactly one of mass= or activity=.")
        if mass is not None:
            grams = _quantity_value(mass, si.g, "mass")
            return ((grams * si.g) / self.atomic_mass).n
        assert activity is not None
        becquerels = _quantity_value(activity, si.Bq, "activity")
        if self.is_stable:
            raise ValueError(
                f"Activity cannot determine the amount of stable {self.name}."
            )
        return ((becquerels * si.Bq) / self.decay_constant).n

    def mass(
        self, *, atoms: float | None = None, activity: Quantity[float] | None = None
    ) -> Quantity[float]:
        """Mass of a pure sample, specifying either atom count or activity."""
        if (atoms is None) == (activity is None):
            raise TypeError("Specify exactly one of atoms= or activity=.")
        if activity is not None:
            count = self.atoms(activity=activity)
        else:
            assert atoms is not None
            count = _nonnegative(atoms, "atoms")
        return (count * self.atomic_mass).to(si.g)


def nuclides(element: str | None = None) -> tuple[Nuclide, ...]:
    """List supported nuclear states, optionally for an element symbol.

    Results follow the evaluation's mass-number/proton-number/state order.
    Lowercase ``n`` selects the free neutron; ``N`` selects nitrogen.
    """
    if element is not None:
        if not isinstance(element, str):
            raise TypeError("element must be a symbol such as 'Rn'.")
        element = element.strip()
        element = element if element == "n" else element.capitalize()
        if element not in {name.split("-")[0] for name in records()}:
            raise ValueError(f"Unknown element symbol {element!r}.")
    return tuple(
        Nuclide(name)
        for name in records()
        if element is None or name.split("-")[0] == element
    )
