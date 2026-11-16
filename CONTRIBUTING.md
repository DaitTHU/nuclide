# Contributing

Use Python 3.12 or newer. From the repository root:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[test]"
python -m pytest
python -m ruff check .
python -m ruff format --check .
pyright
python -m build
python -m twine check --strict dist/*
```

On Windows, activate the environment with `.venv\Scripts\activate`.

Keep physical quantities explicit at public boundaries and make missing nuclear
data distinguishable from stability. New data must include a traceable source,
its version or retrieval date, and the meaning of any uncertainties or omitted
decay channels. Changes to decay calculations should include tests against
analytic results and physical conservation laws within the model's stated scope.

CI runs the test suite on Python 3.12–3.14, checks formatting and typing, and
builds both source and wheel distributions. It also imports the installed wheel
outside the source checkout to check that its nuclear data are included.

The Python implementation is MIT-licensed. The bundled NUBASE2020 evaluation
is distributed under CC BY 3.0; keep its attribution in
[`src/nuclide/data/NOTICE.md`](src/nuclide/data/NOTICE.md) with any redistributed
data. Source distributions include the original compressed evaluation and
`tools/import_nubase.py` so that the packaged JSON can be reproduced.

## Publishing

The repository uses [PyPI Trusted Publishing](https://docs.pypi.org/trusted-publishers/).
Before the first release, configure a pending publisher for project `nuclide` in
your PyPI account with these exact values:

| Setting | Value |
| --- | --- |
| GitHub owner | `DaitTHU` |
| Repository | `nuclide` |
| Workflow filename | `release.yml` |
| Environment | `pypi` |

Create the GitHub environment `pypi` to match that publisher. No PyPI API token
is needed in the repository.

For a release, update the version in `pyproject.toml`, run the checks above,
and commit the release. Push a tag matching that version, such as `v0.1.0`.
The release workflow checks the tag, runs the full CI workflow, and publishes
the distributions built by that run. Verify the resulting PyPI project and
install the released version in a fresh environment before announcing it.
