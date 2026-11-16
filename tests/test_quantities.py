"""Unit-aware calculations, physical limits, and scalar input boundaries."""

import math

import pytest
from unitx import Quantity
from unitx.exceptions import DimensionError
from unitx.si import si

from nuclide import DataUnavailableError, Nuclide


def test_half_life_preserves_nubase_time_conventions():
    assert Nuclide("Rn222").half_life.value_in(si.d) == pytest.approx(3.8215)
    assert Nuclide("Co60").half_life.value_in(si.s) == pytest.approx(
        5.2714 * 31_556_926
    )
    assert Nuclide("He5").half_life.value_in(si.s) == pytest.approx(
        602e-24, rel=1e-12, abs=0
    )


def test_nuclide_quantities_retain_units():
    rn = Nuclide("Rn222")
    assert isinstance(rn.half_life, Quantity)
    assert isinstance(rn.decay_constant, Quantity)
    assert isinstance(rn.atomic_mass, Quantity)
    assert rn.mass_excess.value_in(si.keV) == pytest.approx(16372.0)
    assert rn.atomic_mass.value_in(si.u) == pytest.approx(222.017576, rel=0, abs=1e-6)
    assert rn.decay_constant.value_in(si.s**-1) == pytest.approx(
        math.log(2) / (3.8215 * 86_400)
    )
    assert rn.mean_lifetime.value_in(si.s) == pytest.approx(
        3.8215 * 86_400 / math.log(2)
    )


def test_atomic_mass_carbon_12_and_isomer_energy():
    carbon = Nuclide("C12")
    assert carbon.atomic_mass.value_in(si.u) == 12.0
    assert carbon.mass_excess.value_in(si.keV) == 0.0
    assert carbon.excitation_energy.value_in(si.keV) == 0.0
    assert Nuclide("Tc99m").excitation_energy.value_in(si.keV) == pytest.approx(
        142.6836
    )
    assert Nuclide("Tc99m").atomic_mass > Nuclide("Tc99").atomic_mass


def test_survival_fraction_follows_half_lives_and_unit_conversion():
    rn = Nuclide("Rn222")
    assert rn.survival_fraction(0 * si.s) == 1.0
    assert rn.survival_fraction(rn.half_life) == pytest.approx(0.5)
    assert rn.survival_fraction(2 * rn.half_life) == pytest.approx(0.25)
    assert rn.survival_fraction(1 * si.d) == pytest.approx(
        rn.survival_fraction(24 * si.h)
    )
    assert rn.survival_fraction(10_000 * rn.half_life) == 0.0


def test_atom_mass_activity_conversions_agree():
    rn = Nuclide("Rn222")
    sample_mass = 1 * si.ug
    atom_count = rn.atoms(mass=sample_mass)
    activity = rn.activity(mass=sample_mass)
    assert atom_count > 0
    assert rn.mass(atoms=atom_count).isclose(sample_mass)
    assert rn.mass(activity=activity).isclose(sample_mass)
    assert rn.activity(atoms=atom_count).isclose(activity)
    assert rn.atoms(activity=activity) == pytest.approx(atom_count)
    assert rn.activity(mass=1000 * si.ng).isclose(activity)
    assert (rn.specific_activity * sample_mass).isclose(activity)


def test_molar_mass_reproduces_a_mole_of_atoms():
    carbon = Nuclide("C12")
    mass = carbon.molar_mass * (1 * si.mol)
    assert carbon.atoms(mass=mass) == pytest.approx(6.02214076e23)
    assert carbon.molar_mass.value_in(si.g / si.mol) == pytest.approx(12.0)


def test_zero_sample_is_valid():
    rn = Nuclide("Rn222")
    assert rn.activity(atoms=0).value_in(si.Bq) == 0.0
    assert rn.activity(mass=0 * si.g).value_in(si.Bq) == 0.0
    assert rn.atoms(activity=0 * si.Bq) == 0.0
    assert rn.atoms(mass=0 * si.g) == 0.0
    assert rn.mass(atoms=0).value_in(si.g) == 0.0
    assert rn.mass(activity=0 * si.Bq).value_in(si.g) == 0.0


def test_stable_nuclide_has_zero_activity():
    carbon = Nuclide("C12")
    assert math.isinf(carbon.half_life.value_in(si.s))
    assert math.isinf(carbon.mean_lifetime.value_in(si.s))
    assert carbon.decay_constant.value_in(si.s**-1) == 0.0
    assert carbon.specific_activity.value_in(si.Bq / si.g) == 0.0
    assert carbon.activity(atoms=1e20).value_in(si.Bq) == 0.0
    assert carbon.activity(mass=1 * si.g).value_in(si.Bq) == 0.0
    assert carbon.survival_fraction(1e30 * si.s) == 1.0
    assert carbon.mass(atoms=carbon.atoms(mass=1 * si.g)).isclose(1 * si.g)


@pytest.mark.parametrize("activity", [0 * si.Bq, 1 * si.Bq])
def test_stable_nuclide_cannot_infer_amount_from_activity(activity):
    carbon = Nuclide("C12")
    with pytest.raises(ValueError):
        carbon.atoms(activity=activity)
    with pytest.raises(ValueError):
        carbon.mass(activity=activity)


@pytest.mark.parametrize("value", [-1.0, math.inf, -math.inf, math.nan])
def test_nonphysical_scalar_inputs_are_rejected(value):
    rn = Nuclide("Rn222")
    with pytest.raises(ValueError):
        rn.survival_fraction(value * si.s)
    with pytest.raises(ValueError):
        rn.activity(atoms=value)
    with pytest.raises(ValueError):
        rn.mass(atoms=value)
    with pytest.raises(ValueError):
        rn.activity(mass=value * si.g)
    with pytest.raises(ValueError):
        rn.atoms(mass=value * si.g)
    with pytest.raises(ValueError):
        rn.atoms(activity=value * si.Bq)
    with pytest.raises(ValueError):
        rn.mass(activity=value * si.Bq)


def test_dimension_errors_identify_wrong_physical_inputs():
    rn = Nuclide("Rn222")
    with pytest.raises(DimensionError):
        rn.survival_fraction(1 * si.m)
    with pytest.raises(DimensionError):
        rn.activity(mass=1 * si.s)
    with pytest.raises(DimensionError):
        rn.atoms(mass=1 * si.s)
    with pytest.raises(DimensionError):
        rn.atoms(activity=1 * si.g)
    with pytest.raises(DimensionError):
        rn.mass(activity=1 * si.g)


def test_unknown_half_life_is_not_silently_stable():
    lithium = Nuclide("Li12")
    with pytest.raises(DataUnavailableError):
        _ = lithium.half_life
    with pytest.raises(DataUnavailableError):
        _ = lithium.decay_constant
    with pytest.raises(DataUnavailableError):
        lithium.survival_fraction(1 * si.s)
    assert lithium.atomic_mass.value_in(si.u) > 12.0


@pytest.mark.parametrize("name", ["Ac218m", "Ag95m", "B6"])
def test_lifetime_limit_is_not_an_exact_lifetime(name):
    nuclide = Nuclide(name)
    with pytest.raises(DataUnavailableError):
        _ = nuclide.half_life
    with pytest.raises(DataUnavailableError):
        nuclide.activity(atoms=1.0)


def test_estimated_half_life_can_be_used():
    assert Nuclide("Ac237").half_life_status == "estimated"
    assert Nuclide("Ac237").half_life.value_in(si.s) == 23.0
    assert Nuclide("Ac237").survival_fraction(23 * si.s) == pytest.approx(0.5)


def test_amount_conversion_requires_exactly_one_keyword():
    rn = Nuclide("Rn222")
    with pytest.raises(TypeError):
        rn.activity()
    with pytest.raises(TypeError):
        rn.activity(atoms=1, mass=1 * si.g)
    with pytest.raises(TypeError):
        rn.atoms()
    with pytest.raises(TypeError):
        rn.atoms(mass=1 * si.g, activity=1 * si.Bq)
    with pytest.raises(TypeError):
        rn.mass()
    with pytest.raises(TypeError):
        rn.mass(atoms=1, activity=1 * si.Bq)


def test_time_mass_and_activity_require_units():
    rn = Nuclide("Rn222")
    with pytest.raises(TypeError):
        rn.survival_fraction(1)
    with pytest.raises(TypeError):
        rn.activity(mass=1)
    with pytest.raises(TypeError):
        rn.atoms(activity=1)


@pytest.mark.parametrize("atoms", [True, 1j, [1], "1"])
def test_atom_count_requires_a_real_scalar(atoms):
    rn = Nuclide("Rn222")
    with pytest.raises(TypeError):
        rn.activity(atoms=atoms)
    with pytest.raises(TypeError):
        rn.mass(atoms=atoms)
