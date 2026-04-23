"""
Structure label contract tests.

Verifies the documented label format from structure_summary:
  Part 1: classification
  Part 2: balance tier (semiprimes with 2 distinct primes)
  Part 3: residue family
"""
import json
import subprocess
import sys

import pytest

from primehelix.display.json_output import structure_summary


def cli_json(*args):
    r = subprocess.run(
        [sys.executable, "-m", "primehelix.cli", *args],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, f"non-zero exit: {r.stderr}"
    return json.loads(r.stdout)


# ---------------------------------------------------------------------------
# structure_summary unit tests
# ---------------------------------------------------------------------------

class FakeCoil:
    def __init__(self, bit_gap, balance):
        self.bit_gap = bit_gap
        self.balance = balance


class TestLabelFormat:
    def test_prime_2_is_single_part(self):
        label = structure_summary("prime", coil=None, residue={"prime_family": None})
        assert label == "prime"

    def test_prime_gaussian(self):
        label = structure_summary("prime", coil=None, residue={"prime_family": "gaussian"})
        assert label == "prime | gaussian"

    def test_prime_pythagorean(self):
        label = structure_summary("prime", coil=None, residue={"prime_family": "pythagorean"})
        assert label == "prime | pythagorean"

    def test_composite_is_single_part(self):
        label = structure_summary("composite")
        assert label == "composite"

    def test_invalid_is_single_part(self):
        label = structure_summary("invalid")
        assert label == "invalid"

    def test_semiprime_balanced(self):
        coil = FakeCoil(bit_gap=0, balance=0.05)
        label = structure_summary("semiprime", coil=coil, residue={"semiprime_mod4_pair": "1x1"})
        assert label == "semiprime | balanced | mod4_1x1"

    def test_semiprime_moderate(self):
        coil = FakeCoil(bit_gap=4, balance=2.0)
        label = structure_summary("semiprime", coil=coil, residue={"semiprime_mod4_pair": "1x3"})
        assert label == "semiprime | moderate | mod4_1x3"

    def test_semiprime_lopsided(self):
        coil = FakeCoil(bit_gap=13, balance=87.0)
        label = structure_summary("semiprime", coil=coil, residue={"semiprime_mod4_pair": "1x3"})
        assert label == "semiprime | lopsided | mod4_1x3"

    def test_semiprime_prime_square_no_balance(self):
        label = structure_summary("semiprime", coil=None, residue={"semiprime_mod4_pair": "3x3"})
        assert label == "semiprime | mod4_3x3"

    def test_label_parts_order(self):
        coil = FakeCoil(bit_gap=5, balance=3.0)
        label = structure_summary("semiprime", coil=coil, residue={"semiprime_mod4_pair": "3x3"})
        parts = [p.strip() for p in label.split("|")]
        assert parts[0] == "semiprime"
        assert parts[1] in ("balanced", "moderate", "lopsided")
        assert parts[2].startswith("mod4_")

    def test_all_balance_tiers_covered(self):
        assert "balanced" in structure_summary("semiprime", FakeCoil(0, 0.01), {"semiprime_mod4_pair": "1x1"})
        assert "moderate" in structure_summary("semiprime", FakeCoil(3, 1.0), {"semiprime_mod4_pair": "1x3"})
        assert "lopsided" in structure_summary("semiprime", FakeCoil(15, 200.0), {"semiprime_mod4_pair": "3x3"})


# ---------------------------------------------------------------------------
# End-to-end label consistency: classify and scan agree
# ---------------------------------------------------------------------------

class TestLabelConsistency:
    def test_classify_91_label_has_balance(self):
        d = cli_json("classify", "91", "--json")
        parts = [p.strip() for p in d["structure"].split("|")]
        assert parts[0] == "semiprime"
        assert parts[1] in ("balanced", "moderate", "lopsided")
        assert parts[2].startswith("mod4_")

    def test_classify_1300039_is_lopsided(self):
        d = cli_json("classify", "1300039", "--json")
        assert d["structure"] == "semiprime | lopsided | mod4_1x3"

    def test_scan_labels_match_classify_labels(self):
        scan = cli_json("structure-scan", "--start", "88", "--stop", "95", "--json")
        for n in [91]:
            classify = cli_json("classify", str(n), "--json")
            label = classify["structure"]
            assert label in scan["counts"], (
                f"classify gives '{label}' for n={n} but not in scan counts"
            )


# ---------------------------------------------------------------------------
# Plotting: import works (matplotlib installed in dev)
# ---------------------------------------------------------------------------

class TestPlottingImport:
    def test_matplotlib_importable(self):
        pytest.importorskip("matplotlib", reason="matplotlib not installed")
        from primehelix.display.plots import save_structure_time_series_plot
        assert callable(save_structure_time_series_plot)

    def test_plot_import_error_message(self, monkeypatch):
        import builtins
        real_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if name == "matplotlib":
                raise ImportError("No module named 'matplotlib'")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", mock_import)

        # Re-import to force the lazy loader to run
        import importlib
        import primehelix.display.plots as plots_mod
        importlib.reload(plots_mod)

        with pytest.raises(ImportError, match="pip install"):
            plots_mod._require_matplotlib()
