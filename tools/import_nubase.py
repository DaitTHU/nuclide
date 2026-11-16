#!/usr/bin/env python3
"""Rebuild the bundled NUBASE2020 JSON without third-party dependencies."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "src" / "nuclide" / "data"
SOURCE_URL = "https://nds.iaea.org/amdc/ame2020/nubase_4.mas20.txt"
SOURCE_SHA256 = "1585a5eea86c5e17e90307c7e6e786d060049c4039e392a261ff6db977df9859"
# NUBASE2020 section 2.5 specifies the tropical year, rounded to whole seconds.
YEAR_SECONDS = 31_556_926
TIME_UNITS = {
    "ys": 1e-24,
    "zs": 1e-21,
    "as": 1e-18,
    "fs": 1e-15,
    "ps": 1e-12,
    "ns": 1e-9,
    "us": 1e-6,
    "ms": 1e-3,
    "s": 1.0,
    "m": 60.0,
    "h": 3600.0,
    "d": 86400.0,
    "y": YEAR_SECONDS,
    "ky": 1e3 * YEAR_SECONDS,
    "My": 1e6 * YEAR_SECONDS,
    "Gy": 1e9 * YEAR_SECONDS,
    "Ty": 1e12 * YEAR_SECONDS,
    "Py": 1e15 * YEAR_SECONDS,
    "Ey": 1e18 * YEAR_SECONDS,
    "Zy": 1e21 * YEAR_SECONDS,
    "Yy": 1e24 * YEAR_SECONDS,
}


def number(field: str) -> float | None:
    """Read a central value, leaving its source qualifier in the raw row."""
    value = field.strip()
    if not value:
        return None
    if not re.fullmatch(r"-?\d+(?:\.\d*)?#?", value):
        raise ValueError(f"Unrecognized numeric field: {field!r}")
    return float(value.removesuffix("#"))


def half_life(row: str) -> tuple[float | None, str, str]:
    value = row[69:78].strip()
    unit = row[78:80].strip()
    raw = f"{value} {unit}".strip()
    if value == "stbl":
        return None, "stable", raw
    if value == "p-unst":
        return None, "unstable", raw
    if not value:
        # Some entries have no central value, only a bound in the dT column.
        alternative = row[81:88].strip()
        match = re.fullmatch(r"([<>~]\d+(?:\.\d*)?)\s*([A-Za-z]+)", alternative)
        if match is None:
            if alternative:
                raise ValueError(f"Unrecognized half-life field: {alternative!r}")
            return None, "unknown", raw
        value, unit = match.groups()
        raw = alternative
    match = re.fullmatch(r"([<>~]?)(\d+(?:\.\d*)?)(#?)", value)
    if match is None:
        raise ValueError(f"Unrecognized half-life value: {value!r}")
    qualifier, central, estimated = match.groups()
    status = {
        "<": "upper_limit",
        ">": "lower_limit",
        "~": "approximate",
        "": "estimated" if estimated else "measured",
    }[qualifier]
    seconds = float(central) * TIME_UNITS[unit]
    if seconds <= 0:
        raise ValueError(f"Nonpositive half-life: {raw!r}")
    return seconds, status, raw


def convert(source: bytes) -> dict:
    """Parse exactly the documented fixed-width ground/m/n state records."""
    digest = hashlib.sha256(source).hexdigest()
    if digest != SOURCE_SHA256:
        raise ValueError(
            "Source differs from the reviewed NUBASE2020 snapshot. "
            "Review the new data and update SOURCE_SHA256 before importing."
        )
    records = {}
    excluded = []
    source_rows = 0
    for line in source.decode("utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        source_rows += 1
        A, Z, index = int(line[0:3]), int(line[4:7]), int(line[7])
        state = line[16:17].strip()
        symbol_match = re.fullmatch(r"(\d+)([A-Za-z]+)", line[11:16].strip())
        if symbol_match is None or int(symbol_match[1]) != A:
            raise ValueError(f"Invalid isotope label: {line!r}")
        symbol = symbol_match[2]
        name = f"{symbol}-{A}{state}"
        if index not in {0, 1, 2}:
            excluded.append({"name": name, "reason": "state_outside_ground_m_n"})
            continue
        if state != {0: "", 1: "m", 2: "n"}[index]:
            raise ValueError(f"Inconsistent state label: {line!r}")
        origin = line[65:67]
        if origin in {"RN", "EU"}:
            excluded.append(
                {
                    "name": name,
                    "reason": {
                        "RN": "reported_nonexistent",
                        "EU": "existence_uncertain",
                    }[origin],
                }
            )
            continue
        seconds, status, raw = half_life(line)
        record = {
            "name": name,
            "Z": Z,
            "A": A,
            "state": state,
            "mass_excess_kev": number(line[18:31]),
            "mass_excess_uncertainty_kev": number(line[31:42]),
            "mass_excess_estimated": "#" in line[18:42],
            "half_life_seconds": seconds,
            "half_life_status": status,
            "half_life_raw": raw,
            "half_life_uncertainty_raw": line[81:88].strip(),
            "decay_modes_raw": line[119:209].strip(),
            "spin_parity": line[88:102].strip(),
            "excitation_energy_kev": number(line[42:54]),
            "raw_line": line,
        }
        if name in records:
            raise ValueError(f"Duplicate nuclide: {name}")
        records[name] = record
    counts = Counter(record["half_life_status"] for record in records.values())
    return {
        "metadata": {
            "dataset": "NUBASE2020",
            "schema_version": 1,
            "source_url": SOURCE_URL,
            "source_sha256": digest,
            "retrieved_on": "2026-11-16",
            "citation": (
                "F.G. Kondev, M. Wang, W.J. Huang, S. Naimi, and G. Audi, "
                "The NUBASE2020 evaluation of nuclear physics properties, "
                "Chinese Physics C 45, 030001 (2021)."
            ),
            "doi": "10.1088/1674-1137/abddae",
            "license": "CC-BY-3.0",
            "year_seconds": YEAR_SECONDS,
            "source_row_count": source_rows,
            "record_count": len(records),
            "state_counts": dict(
                Counter(record["state"] or "ground" for record in records.values())
            ),
            "half_life_status_counts": dict(counts),
            "excluded": excluded,
        },
        "nuclides": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "source", nargs="?", type=Path, default=DATA / "nubase_4.mas20.txt.gz"
    )
    parser.add_argument("--output", type=Path, default=DATA / "nubase2020.json")
    args = parser.parse_args()
    source = args.source.read_bytes()
    if args.source.suffix == ".gz":
        source = gzip.decompress(source)
    payload = convert(source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            payload, ensure_ascii=False, allow_nan=False, indent=2, sort_keys=True
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {payload['metadata']['record_count']} records to {args.output}")


if __name__ == "__main__":
    main()
