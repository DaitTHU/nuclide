"""Public nuclear identities and accepted name forms."""

import pickle
from dataclasses import FrozenInstanceError

import pytest

from nuclide import Nuclide, nuclides


@pytest.mark.parametrize(
    "alias",
    ["Rn222", "Rn-222", "222Rn", "222-Rn", "rn222", "RN-222", "  Rn222  "],
)
def test_ground_state_aliases_share_identity(alias):
    reference = Nuclide("Rn-222")
    other = Nuclide(alias)
    assert other == reference
    assert hash(other) == hash(reference)
    assert {reference: "radon"}[other] == "radon"
    assert str(other) == "Rn-222"
    assert repr(other) == "Nuclide('Rn-222')"
    assert (other.Z, other.A, other.N, other.element, other.state) == (
        86,
        222,
        136,
        "Rn",
        "",
    )


@pytest.mark.parametrize("alias", ["Tc99m", "Tc-99m1", "99mTc", "99m1-Tc"])
def test_first_isomer_is_a_distinct_identity(alias):
    isomer = Nuclide(alias)
    assert isomer == Nuclide("Tc-99m")
    assert isomer != Nuclide("Tc99")
    assert isomer.state == "m"
    assert len({isomer, Nuclide("Tc99")}) == 2


@pytest.mark.parametrize("alias", ["Hf178n", "Hf-178m2", "178nHf", "178m2-Hf"])
def test_second_isomer_aliases(alias):
    isomer = Nuclide(alias)
    assert isomer == Nuclide("Hf-178n")
    assert isomer.state == "n"
    assert len({isomer, Nuclide("Hf178m"), Nuclide("Hf178")}) == 3


@pytest.mark.parametrize("alias", ["n1", "n-1", "1n", "1-n"])
def test_free_neutron_identity(alias):
    neutron = Nuclide(alias)
    assert neutron == Nuclide("n-1")
    assert (neutron.Z, neutron.A, neutron.N, neutron.element) == (0, 1, 1, "n")


def test_identity_is_immutable_and_survives_pickle():
    nuclide = Nuclide("Tc99m")
    with pytest.raises(FrozenInstanceError):
        nuclide.name = "Tc-99"
    assert pickle.loads(pickle.dumps(nuclide)) == nuclide


@pytest.mark.parametrize(
    "name", ["", "Rn", "222", "Rn0", "Rn-0222", "Rn222!", "Tc99m3"]
)
def test_malformed_name_is_rejected(name):
    with pytest.raises(ValueError, match="Invalid nuclide name"):
        Nuclide(name)


@pytest.mark.parametrize("name", ["Xx222", "Rn999", "C1"])
def test_unlisted_nuclear_identity_is_rejected(name):
    with pytest.raises(ValueError, match="not in the bundled NUBASE2020 data"):
        Nuclide(name)


@pytest.mark.parametrize("name", [None, 222, (86, 222)])
def test_nonstring_name_is_rejected(name):
    with pytest.raises(TypeError, match="requires a name"):
        Nuclide(name)


def test_catalogue_identities_round_trip_without_duplicates():
    catalogue = nuclides()
    assert catalogue
    assert len(catalogue) == len(set(catalogue))
    assert all(Nuclide(str(item)) == item for item in catalogue)
    assert all(item.A == item.Z + item.N for item in catalogue)


def test_element_filter_uses_normalized_symbols():
    radon = nuclides("Rn")
    assert radon == nuclides(" rn ")
    assert Nuclide("Rn222") in radon
    assert all(item.element == "Rn" for item in radon)


def test_free_neutron_filter_is_distinct_from_nitrogen():
    assert nuclides("n") == (Nuclide("n1"),)
    nitrogen = nuclides("N")
    assert Nuclide("N14") in nitrogen
    assert all(item.element == "N" for item in nitrogen)


@pytest.mark.parametrize("symbol", ["", "radon", "Xx", "Rn222"])
def test_unknown_element_filter_is_rejected(symbol):
    with pytest.raises(ValueError, match="Unknown element symbol"):
        nuclides(symbol)


def test_nonstring_element_filter_is_rejected():
    with pytest.raises(TypeError, match="element must be a symbol"):
        nuclides(86)
