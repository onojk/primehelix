import json
import subprocess
import sys


def run_cli(*args):
    cmd = [sys.executable, "-m", "primehelix.cli", *args]
    completed = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return completed.stdout


def test_compare_ranges_json_shape():
    out = run_cli(
        "compare-ranges",
        "--a-start", "10",
        "--a-stop", "20",
        "--b-start", "20",
        "--b-stop", "30",
        "--json",
    )
    data = json.loads(out)

    assert data["command"] == "compare-ranges"

    assert data["a"]["start"] == 10
    assert data["a"]["stop"] == 20
    assert data["a"]["total"] == 10

    assert data["b"]["start"] == 20
    assert data["b"]["stop"] == 30
    assert data["b"]["total"] == 10

    assert isinstance(data["a"]["counts"], dict)
    assert isinstance(data["b"]["counts"], dict)
    assert isinstance(data["rows"], list)
    assert len(data["rows"]) > 0


def test_compare_ranges_json_rows_have_expected_keys():
    out = run_cli(
        "compare-ranges",
        "--a-start", "10",
        "--a-stop", "20",
        "--b-start", "20",
        "--b-stop", "30",
        "--json",
    )
    data = json.loads(out)

    row = data["rows"][0]
    assert "structure" in row
    assert "a_count" in row
    assert "a_percent" in row
    assert "b_count" in row
    assert "b_percent" in row
    assert "delta" in row
    assert "ratio" in row


def test_compare_ranges_json_contains_composite_row():
    out = run_cli(
        "compare-ranges",
        "--a-start", "10",
        "--a-stop", "20",
        "--b-start", "20",
        "--b-stop", "30",
        "--json",
    )
    data = json.loads(out)

    composite = next(row for row in data["rows"] if row["structure"] == "composite")
    assert composite["a_count"] == 3
    assert composite["b_count"] == 4
    assert composite["a_percent"] == 30.0
    assert composite["b_percent"] == 40.0
    assert composite["delta"] == 10.0
    assert composite["ratio"] == "1.33x"


def test_compare_ranges_top_delta_json_trimmed():
    out = run_cli(
        "compare-ranges",
        "--a-start", "10",
        "--a-stop", "20",
        "--b-start", "20",
        "--b-stop", "30",
        "--top-delta", "5",
        "--json",
    )
    data = json.loads(out)

    assert data["top_delta"] == 5
    assert len(data["rows"]) == 5


def test_compare_ranges_top_delta_includes_new_or_zero_ratio():
    out = run_cli(
        "compare-ranges",
        "--a-start", "10",
        "--a-stop", "20",
        "--b-start", "20",
        "--b-stop", "30",
        "--top-delta", "5",
        "--json",
    )
    data = json.loads(out)

    ratios = {row["ratio"] for row in data["rows"]}
    assert ("new" in ratios) or ("0.00x" in ratios)


def test_compare_ranges_only_classification_prime_json():
    out = run_cli(
        "compare-ranges",
        "--a-start", "10",
        "--a-stop", "20",
        "--b-start", "20",
        "--b-stop", "30",
        "--only-classification", "prime",
        "--json",
    )
    data = json.loads(out)

    assert data["only_classification"] == "prime"
    assert all(row["structure"].startswith("prime") for row in data["rows"])


def test_compare_ranges_only_structure_gaussian_json():
    out = run_cli(
        "compare-ranges",
        "--a-start", "10",
        "--a-stop", "20",
        "--b-start", "20",
        "--b-stop", "30",
        "--only-structure", "gaussian",
        "--json",
    )
    data = json.loads(out)

    assert data["only_structure"] == "gaussian"
    assert len(data["rows"]) >= 1
    assert all("gaussian" in row["structure"].lower() for row in data["rows"])
