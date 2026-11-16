"""Check the evaluated data boundary, including meaningful source exceptions."""

import gzip
import hashlib
import importlib.util
import json
import math
from collections import Counter
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "src" / "nuclide" / "data"


@pytest.fixture(scope="module")
def dataset():
    return json.loads((DATA / "nubase2020.json").read_text(encoding="utf-8"))


def test_snapshot_checksum_and_reproducibility(dataset):
    source = gzip.decompress((DATA / "nubase_4.mas20.txt.gz").read_bytes())
    assert hashlib.sha256(source).hexdigest() == dataset["metadata"]["source_sha256"]
    spec = importlib.util.spec_from_file_location(
        "import_nubase", ROOT / "tools" / "import_nubase.py"
    )
    assert spec is not None and spec.loader is not None
    importer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(importer)
    assert importer.convert(source) == dataset
    with pytest.raises(ValueError, match="differs from the reviewed"):
        importer.convert(source + b"\n")


def test_dataset_coverage_and_numeric_integrity(dataset):
    records = dataset["nuclides"]
    assert len(records) == dataset["metadata"]["record_count"] == 5375
    assert dataset["metadata"]["source_row_count"] == 5843
    assert Counter(r["state"] or "ground" for r in records.values()) == {
        "ground": 3555,
        "m": 1366,
        "n": 454,
    }
    for name, record in records.items():
        assert name == record["name"]
        assert record["A"] >= record["Z"] >= 0
        assert record["raw_line"]
        seconds = record["half_life_seconds"]
        if record["half_life_status"] in {"stable", "unknown", "unstable"}:
            assert seconds is None
        else:
            assert seconds > 0 and math.isfinite(seconds)
        for field in (
            "mass_excess_kev",
            "mass_excess_uncertainty_kev",
            "excitation_energy_kev",
        ):
            value = record[field]
            assert value is None or math.isfinite(value)


def test_reference_nuclides_and_source_year(dataset):
    records = dataset["nuclides"]
    assert records["Rn-222"]["half_life_seconds"] == pytest.approx(3.8215 * 86400)
    assert records["Tc-99m"]["half_life_seconds"] == pytest.approx(6.0066 * 3600)
    assert records["Tc-99m"]["excitation_energy_kev"] == 142.6836
    assert records["Tc-99m"]["state"] == "m"
    assert records["Tc-99"]["state"] == ""
    assert records["H-1"]["half_life_status"] == "stable"
    assert records["H-1"]["mass_excess_kev"] == 7288.971064
    assert records["C-12"]["mass_excess_kev"] == 0
    assert records["H-3"]["half_life_seconds"] == pytest.approx(12.32 * 31_556_926)
    assert dataset["metadata"]["year_seconds"] == 31_556_926


def test_limits_approximations_and_missing_values_remain_distinct(dataset):
    records = dataset["nuclides"]
    assert records["B-18"]["half_life_raw"] == "<26 ns"
    assert records["B-18"]["half_life_status"] == "upper_limit"
    assert records["B-18"]["half_life_seconds"] == pytest.approx(26e-9)
    assert records["Ac-218m"]["half_life_status"] == "lower_limit"
    assert records["Ac-216n"]["half_life_status"] == "approximate"
    assert records["Ac-237"]["half_life_status"] == "estimated"
    assert records["Ag-133"]["half_life_status"] == "unknown"
    assert records["Li-3"]["half_life_status"] == "unstable"
    assert records["Li-3"]["mass_excess_estimated"] is True


def test_nonexistent_and_disputed_states_excluded(dataset):
    excluded = {row["name"]: row["reason"] for row in dataset["metadata"]["excluded"]}
    assert excluded["Cu-76m"] == "reported_nonexistent"
    assert excluded["Pm-138m"] == "existence_uncertain"
    assert "Cu-76m" not in dataset["nuclides"]
    assert "Pm-138m" not in dataset["nuclides"]
    assert Counter(excluded.values()) == {
        "reported_nonexistent": 11,
        "existence_uncertain": 13,
        "state_outside_ground_m_n": 444,
    }


def test_decay_modes_preserve_inclusive_source_notation(dataset):
    assert dataset["nuclides"]["He-8"]["decay_modes_raw"] == "B-=100;B-n=16 1;B-t=0.9 1"
