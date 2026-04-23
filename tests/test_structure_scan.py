import json
import subprocess
import sys


def run_cli(*args):
    cmd = [sys.executable, "-m", "primehelix.cli", *args]
    completed = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return completed.stdout


def test_structure_scan_json_counts_exist():
    out = run_cli("structure-scan", "--start", "10", "--stop", "20", "--json")
    data = json.loads(out)

    assert data["command"] == "structure-scan"
    assert data["start"] == 10
    assert data["stop"] == 20
    assert data["total"] == 10
    assert isinstance(data["counts"], dict)
    assert len(data["counts"]) > 0


def test_structure_scan_contains_prime_families():
    out = run_cli("structure-scan", "--start", "10", "--stop", "20", "--json")
    data = json.loads(out)

    labels = list(data["counts"].keys())
    assert any("prime | pythagorean" in label for label in labels)
    assert any("prime | gaussian" in label for label in labels)
