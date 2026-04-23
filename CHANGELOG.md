# Changelog

## [0.3.0] — 2026-04-23

### Added

- **`primehelix.analysis` module** — `scan_range()`, `compare_summaries()`, and `build_time_series()` are now importable as a library without the CLI. Results are typed dataclasses.
- **`primehelix.schema` module** — `ScanResult`, `CompareRow`, `WindowSummary`, `TimeSeriesResult` with serialization helpers (`to_json_dict()`, `to_dict()`).
- **`--export-json`** on all five commands — write the full JSON payload to a file alongside or instead of stdout.
- **`--export-csv`** on `structure-scan`, `compare-ranges`, and `structure-time-series`.
- **`--fast` flag** on `structure-scan`, `compare-ranges`, and `structure-time-series` — skips geometry, returns classification-only labels (`detail="classification"` in the API). ~10% faster on unfiltered scans; much faster with `--only-classification`.
- **`--only-classification` validation** — unknown values are rejected immediately with a clear error listing valid options (`prime`, `semiprime`, `composite`, `invalid`). Uses `click.Choice` in the CLI; raises `ValueError` from the library API.
- **Shannon entropy** in `structure-scan` and `compare-ranges` JSON output (`entropy`, `entropy_delta`).
- **`--profile` flag** on `structure-scan` — adds factorization method distribution to output.
- **Golden fixture regression tests** locking exact label counts, entropy, and compare deltas for `[1, 10000)`.
- **10M scaling data** in README Findings — lopsided share grows from 73.2% at 1M to 78.5% at 10M.
- **Python API section** in README with importable usage examples.

### Changed

- Version is now a single source of truth in `primehelix/__init__.py`; `cli.py` reads `__version__` rather than hardcoding it.
- `cli.py` payload assembly uses schema serialization helpers — command bodies are thinner.
- Lazy geometry in `scan_range()` — residue profiling and coil computation are skipped for integers that fail the `--only-classification` filter.
- `detail="fast"` renamed to `detail="classification"` in the API (CLI flag `--fast` unchanged).
- README restructured to lead with empirical findings (result → tool → engine).

### Fixed

- `classify --json` now always includes the balance tier in the `structure` label (was missing without `--coil`).
- `factor --json` no longer includes `classification: null`; `steps` is `[]` without `--verbose`.
- Float and non-integer input to `classify` now produces a clean error instead of a raw traceback.
- `primehelix/__init__.py` version was `0.1.0` while `pyproject.toml` and CLI reported `0.2.0` — fixed.

## [0.2.0] — 2026-04-21

Initial public release. Five commands: `classify`, `factor`, `structure-scan`, `compare-ranges`, `structure-time-series`. Baillie-PSW primality, full factoring pipeline (trial → Pollard p−1 → Williams p+1 → Pollard Rho → ECM → Quadratic Sieve), structure labels, ASCII helix visualization, JSON output, matplotlib plots.
