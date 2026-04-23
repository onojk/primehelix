import json
import subprocess
import sys


def run_cli(*args):
    cmd = [sys.executable, "-m", "primehelix.cli", *args]
    completed = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return completed.stdout


def test_classify_json_schema():
    out = run_cli("classify", "2147483646", "--json")
    data = json.loads(out)

    assert data["command"] == "classify"
    assert data["n"] == 2147483646
    assert data["classification"] == "composite"
    assert data["method"] == "rho"
    assert data["complete"] is True
    assert data["factors"]["331"] == 1
    assert data["prime_factors"] == [2, 3, 3, 7, 11, 31, 151, 331]
    assert data["factorization"] == "2 * 3^2 * 7 * 11 * 31 * 151 * 331"
    assert any(s.startswith("rho:") for s in data["steps"])


def test_factor_json_schema():
    out = run_cli("factor", "2147483646", "--json")
    data = json.loads(out)

    assert data["command"] == "factor"
    assert data["n"] == 2147483646
    assert data["classification"] is None
    assert data["method"] == "rho"
    assert data["complete"] is True
    assert data["factors"]["151"] == 1
    assert data["factors"]["331"] == 1
    assert data["factorization"] == "2 * 3^2 * 7 * 11 * 31 * 151 * 331"
    assert "trial: 31" in data["steps"]


def test_classify_coil_json_contains_insight():
    out = run_cli("classify", "1300039", "--coil", "--json")
    data = json.loads(out)

    assert data["command"] == "classify"
    assert data["classification"] == "semiprime"
    assert data["coil"]["p"] == 13
    assert data["coil"]["q"] == 100003
    assert data["coil"]["bit_gap"] == 13
    assert "insight" in data["coil"]
    assert isinstance(data["coil"]["insight"], str)
    assert len(data["coil"]["insight"]) > 0
