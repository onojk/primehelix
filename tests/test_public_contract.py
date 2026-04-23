"""
Interface stability tests — verify the public contract documented in README.

These tests only check outward promises (CLI flags exist, JSON keys are present,
output is parseable). They do not test internals or specific algorithm choices.
"""
import json
import subprocess
import sys

import pytest


def cli(*args):
    r = subprocess.run(
        [sys.executable, "-m", "primehelix.cli", *args],
        capture_output=True, text=True,
    )
    return r


def cli_json(*args):
    r = cli(*args)
    assert r.returncode == 0, f"non-zero exit: {r.stderr}"
    return json.loads(r.stdout)


# ---------------------------------------------------------------------------
# classify
# ---------------------------------------------------------------------------

class TestClassifyContract:
    def test_flags_exist(self):
        r = cli("classify", "--help")
        for flag in ("--coil", "--helix", "--residue", "--json"):
            assert flag in r.stdout, f"missing flag {flag}"

    def test_json_required_keys(self):
        d = cli_json("classify", "91", "--json")
        for key in ("command", "n", "classification", "factors",
                    "factorization", "method", "complete", "structure"):
            assert key in d, f"missing key {key}"

    def test_json_command_field(self):
        d = cli_json("classify", "91", "--json")
        assert d["command"] == "classify"

    def test_json_prime(self):
        d = cli_json("classify", "97", "--json")
        assert d["classification"] == "prime"
        assert d["complete"] is True

    def test_json_semiprime(self):
        d = cli_json("classify", "91", "--json")
        assert d["classification"] == "semiprime"
        assert d["factors"]["7"] == 1
        assert d["factors"]["13"] == 1

    def test_json_composite(self):
        d = cli_json("classify", "210", "--json")
        assert d["classification"] == "composite"

    def test_json_invalid_zero(self):
        d = cli_json("classify", "0", "--json")
        assert d["classification"] == "invalid"

    def test_json_invalid_one(self):
        d = cli_json("classify", "1", "--json")
        assert d["classification"] == "invalid"

    def test_json_coil_contains_insight(self):
        d = cli_json("classify", "1300039", "--coil", "--json")
        assert "coil" in d
        assert "insight" in d["coil"]
        assert d["coil"]["bit_gap"] == 13

    def test_json_residue_present(self):
        d = cli_json("classify", "91", "--json")
        assert "residue" in d

    def test_structure_label_format(self):
        d = cli_json("classify", "91", "--json")
        assert "|" in d["structure"] or d["structure"] in ("prime", "composite", "invalid")

    def test_helix_flag_accepted(self):
        r = cli("classify", "91", "--helix")
        assert r.returncode == 0

    def test_float_input_clean_error(self):
        r = cli("classify", "3.7", "--json")
        assert r.returncode != 0
        assert "integer" in r.stderr.lower()

    def test_non_numeric_input_clean_error(self):
        r = cli("classify", "abc", "--json")
        assert r.returncode != 0


# ---------------------------------------------------------------------------
# factor
# ---------------------------------------------------------------------------

class TestFactorContract:
    def test_flags_exist(self):
        r = cli("factor", "--help")
        for flag in ("--verbose", "--budget", "--json"):
            assert flag in r.stdout, f"missing flag {flag}"

    def test_json_required_keys(self):
        d = cli_json("factor", "91", "--json")
        for key in ("command", "n", "factors", "factorization", "method", "complete"):
            assert key in d, f"missing key {key}"

    def test_json_command_field(self):
        d = cli_json("factor", "91", "--json")
        assert d["command"] == "factor"

    def test_json_complete_true(self):
        d = cli_json("factor", "2147483646", "--json")
        assert d["complete"] is True
        assert d["factors"]["331"] == 1
        assert d["factors"]["151"] == 1

    def test_verbose_shows_steps(self):
        d = cli_json("factor", "91", "--json", "--verbose")
        assert isinstance(d["steps"], list)
        assert len(d["steps"]) > 0

    def test_without_verbose_steps_empty(self):
        d = cli_json("factor", "91", "--json")
        assert d.get("steps", []) == [] or isinstance(d["steps"], list)

    def test_plain_output_is_not_repr(self):
        r = cli("factor", "91")
        assert r.returncode == 0
        assert "FactorResult(" not in r.stdout


# ---------------------------------------------------------------------------
# structure-scan
# ---------------------------------------------------------------------------

class TestStructureScanContract:
    def test_flags_exist(self):
        r = cli("structure-scan", "--help")
        for flag in ("--start", "--stop", "--json",
                     "--only-classification", "--only-structure", "--fast"):
            assert flag in r.stdout, f"missing flag {flag}"

    def test_json_required_keys(self):
        d = cli_json("structure-scan", "--start", "1", "--stop", "100", "--json")
        for key in ("command", "start", "stop", "total", "entropy", "counts"):
            assert key in d, f"missing key {key}"

    def test_json_entropy_is_nonneg_float(self):
        d = cli_json("structure-scan", "--start", "1", "--stop", "100", "--json")
        assert isinstance(d["entropy"], float)
        assert d["entropy"] >= 0.0

    def test_json_counts_nonempty(self):
        d = cli_json("structure-scan", "--start", "1", "--stop", "100", "--json")
        assert len(d["counts"]) > 0

    def test_only_classification_filter(self):
        d = cli_json("structure-scan", "--start", "1", "--stop", "200",
                     "--only-classification", "prime", "--json")
        for label in d["counts"]:
            assert label.startswith("prime"), f"non-prime label: {label}"

    def test_invalid_classification_rejected(self):
        r = cli("structure-scan", "--start", "1", "--stop", "100",
                "--only-classification", "banana", "--json")
        assert r.returncode != 0
        assert "banana" in r.stderr or "invalid" in r.stderr.lower()

    def test_only_classification_case_insensitive(self):
        d = cli_json("structure-scan", "--start", "1", "--stop", "200",
                     "--only-classification", "Prime", "--json")
        for label in d["counts"]:
            assert label.startswith("prime"), f"non-prime label: {label}"

    def test_profile_flag_adds_methods(self):
        d = cli_json("structure-scan", "--start", "1", "--stop", "200",
                     "--profile", "--json")
        assert "methods" in d
        assert isinstance(d["methods"], dict)
        assert sum(d["methods"].values()) > 0

    def test_without_profile_no_methods(self):
        d = cli_json("structure-scan", "--start", "1", "--stop", "200", "--json")
        assert "methods" not in d


# ---------------------------------------------------------------------------
# compare-ranges
# ---------------------------------------------------------------------------

class TestCompareRangesContract:
    def test_flags_exist(self):
        r = cli("compare-ranges", "--help")
        for flag in ("--a-start", "--a-stop", "--b-start", "--b-stop",
                     "--top-delta", "--json", "--only-classification"):
            assert flag in r.stdout, f"missing flag {flag}"

    def test_json_required_keys(self):
        d = cli_json("compare-ranges",
                     "--a-start", "1", "--a-stop", "100",
                     "--b-start", "100", "--b-stop", "200", "--json")
        for key in ("command", "rows", "entropy_delta"):
            assert key in d, f"missing key {key}"
        assert "entropy" in d["a"] and "entropy" in d["b"]

    def test_json_entropy_delta(self):
        d = cli_json("compare-ranges",
                     "--a-start", "1", "--a-stop", "100",
                     "--b-start", "100", "--b-stop", "200", "--json")
        assert isinstance(d["entropy_delta"], float)

    def test_json_rows_have_expected_keys(self):
        d = cli_json("compare-ranges",
                     "--a-start", "1", "--a-stop", "100",
                     "--b-start", "100", "--b-stop", "200", "--json")
        for row in d["rows"]:
            for key in ("structure", "a_count", "b_count", "delta"):
                assert key in row, f"row missing key {key}"


# ---------------------------------------------------------------------------
# structure-time-series
# ---------------------------------------------------------------------------

class TestStructureTimeSeriesContract:
    def test_flags_exist(self):
        r = cli("structure-time-series", "--help")
        for flag in ("--start", "--stop", "--window", "--step",
                     "--plot", "--json", "--only-classification"):
            assert flag in r.stdout, f"missing flag {flag}"

    def test_plot_is_optional(self):
        r = cli("structure-time-series",
                "--start", "1", "--stop", "500",
                "--window", "250", "--step", "250",
                "--only-classification", "prime")
        assert r.returncode == 0

    def test_json_required_keys(self):
        d = cli_json("structure-time-series",
                     "--start", "1", "--stop", "500",
                     "--window", "250", "--step", "250",
                     "--only-classification", "prime", "--json")
        for key in ("command", "start", "stop", "labels", "windows"):
            assert key in d, f"missing key {key}"
