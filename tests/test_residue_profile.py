import json
import subprocess
import sys

from primehelix.geometry.residue import residue_profile


def run_cli(*args):
    cmd = [sys.executable, "-m", "primehelix.cli", *args]
    completed = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return completed.stdout


def test_prime_family_labels():
    p13 = residue_profile(13, {13: 1}, classification="prime")
    p7 = residue_profile(7, {7: 1}, classification="prime")

    assert p13["prime_family"] == "pythagorean"
    assert p7["prime_family"] == "gaussian"


def test_semiprime_mod4_pair():
    prof = residue_profile(1300039, {13: 1, 100003: 1}, classification="semiprime")
    assert prof["semiprime_mod4_pair"] == "1x3"


def test_classify_json_contains_residue():
    out = run_cli("classify", "1300039", "--coil", "--json")
    data = json.loads(out)

    assert "residue" in data
    assert data["residue"]["n_mod_4"] == 3
    assert data["residue"]["semiprime_mod4_pair"] == "1x3"
    assert data["structure"] == "semiprime | lopsided | mod4_1x3"


def test_prime_json_structure_summary():
    out = run_cli("classify", "13", "--json")
    data = json.loads(out)

    assert data["classification"] == "prime"
    assert "prime" in data["structure"]
    assert "pythagorean" in data["structure"]
