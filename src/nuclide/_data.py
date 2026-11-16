"""Load the bundled evaluation once, without network access."""

import json
from functools import cache
from importlib.resources import files
from typing import TypedDict


class Record(TypedDict):
    name: str
    Z: int
    A: int
    state: str
    mass_excess_kev: float | None
    mass_excess_uncertainty_kev: float | None
    mass_excess_estimated: bool
    half_life_seconds: float | None
    half_life_status: str
    half_life_raw: str
    half_life_uncertainty_raw: str
    decay_modes_raw: str
    spin_parity: str
    excitation_energy_kev: float | None
    raw_line: str


@cache
def records() -> dict[str, Record]:
    text = files("nuclide").joinpath("data/nubase2020.json").read_text(encoding="utf-8")
    data: dict[str, Record] = json.loads(text)["nuclides"]
    return dict(
        sorted(
            data.items(),
            key=lambda item: (item[1]["A"], item[1]["Z"], item[1]["state"]),
        )
    )
