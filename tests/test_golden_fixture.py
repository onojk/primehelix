"""
Golden-fixture regression tests for [1, 10000).

These lock exact label counts, entropy, and compare-delta directions so that
any change to factorization, labeling, or scan logic is caught immediately.
They call analysis.scan_range / analysis.compare_summaries directly —
no subprocess, no CLI layer.

If a finding changes intentionally (e.g., a bug fix in the label grammar),
update the fixture values here and document why in the commit message.
"""
import pytest
from primehelix.analysis import scan_range, compare_summaries
from primehelix.display.json_output import label_entropy
from primehelix.schema import VALID_CLASSIFICATIONS


FIXTURE_RANGE = (1, 10000)

# Exact label counts for [1, 10000) — locked 2026-04-23
GOLDEN_COUNTS = {
    "composite": 6144,
    "invalid": 1,
    "prime": 1,
    "prime | gaussian": 619,
    "prime | pythagorean": 609,
    "semiprime | balanced | mod4_1x1": 5,
    "semiprime | balanced | mod4_1x3": 16,
    "semiprime | balanced | mod4_3x3": 5,
    "semiprime | lopsided | mod4_1x1": 102,
    "semiprime | lopsided | mod4_1x3": 355,
    "semiprime | lopsided | mod4_2x1": 308,
    "semiprime | lopsided | mod4_2x3": 315,
    "semiprime | lopsided | mod4_3x3": 255,
    "semiprime | mod4_1x1": 11,
    "semiprime | mod4_2x2": 1,
    "semiprime | mod4_3x3": 13,
    "semiprime | moderate | mod4_1x1": 246,
    "semiprime | moderate | mod4_1x3": 600,
    "semiprime | moderate | mod4_2x1": 21,
    "semiprime | moderate | mod4_2x3": 24,
    "semiprime | moderate | mod4_3x3": 348,
}

GOLDEN_TOTAL = 9999
GOLDEN_ENTROPY = 2.2477


@pytest.fixture(scope="module")
def full_scan():
    start, stop = FIXTURE_RANGE
    return scan_range(start, stop)


@pytest.fixture(scope="module")
def half_scans():
    s1 = scan_range(1, 5000)
    s2 = scan_range(5000, 10000)
    return s1, s2


class TestGoldenCounts:
    def test_total(self, full_scan):
        assert full_scan.total == GOLDEN_TOTAL

    def test_exact_label_counts(self, full_scan):
        assert dict(full_scan.counts) == GOLDEN_COUNTS

    def test_no_extra_labels(self, full_scan):
        extra = set(full_scan.counts.keys()) - set(GOLDEN_COUNTS.keys())
        assert not extra, f"unexpected labels appeared: {extra}"

    def test_no_missing_labels(self, full_scan):
        missing = set(GOLDEN_COUNTS.keys()) - set(full_scan.counts.keys())
        assert not missing, f"expected labels disappeared: {missing}"

    def test_entropy(self, full_scan):
        entropy = label_entropy(full_scan.counts, full_scan.total)
        assert entropy == GOLDEN_ENTROPY

    def test_prime_count(self, full_scan):
        # 1229 primes in [1, 10000) — consistent with prime counting function
        gaussian = full_scan.counts.get("prime | gaussian", 0)
        pythagorean = full_scan.counts.get("prime | pythagorean", 0)
        prime_2 = full_scan.counts.get("prime", 0)
        assert gaussian + pythagorean + prime_2 == 1229

    def test_invalid_count(self, full_scan):
        # Only n=1 is invalid in [1, 10000)
        assert full_scan.counts.get("invalid", 0) == 1

    def test_prime_residue_symmetry(self, full_scan):
        gaussian = full_scan.counts["prime | gaussian"]
        pythagorean = full_scan.counts["prime | pythagorean"]
        total_odd = gaussian + pythagorean
        # Neither family dominates by more than 2%
        assert abs(gaussian - pythagorean) / total_odd < 0.02


class TestGoldenCompare:
    def test_lopsided_gains_in_upper_half(self, half_scans):
        s1, s2 = half_scans
        rows = compare_summaries(s1, s2)
        by_label = {r.structure: r.delta for r in rows}
        # Lopsided semiprimes gain share as range grows
        assert by_label["semiprime | lopsided | mod4_1x3"] > 0

    def test_moderate_shrinks_in_upper_half(self, half_scans):
        s1, s2 = half_scans
        rows = compare_summaries(s1, s2)
        by_label = {r.structure: r.delta for r in rows}
        # Moderate semiprimes lose share as range grows
        assert by_label["semiprime | moderate | mod4_1x3"] < 0

    def test_entropy_decreases_in_upper_half(self, half_scans):
        s1, s2 = half_scans
        ea = label_entropy(s1.counts, s1.total)
        eb = label_entropy(s2.counts, s2.total)
        # Upper half is less diverse — lopsided dominance concentrates the distribution
        assert eb < ea

    def test_entropy_values(self, half_scans):
        s1, s2 = half_scans
        assert label_entropy(s1.counts, s1.total) == 2.3308
        assert label_entropy(s2.counts, s2.total) == 2.1401


class TestFastMode:
    def test_fast_total_matches_full(self):
        fast = scan_range(1, 1000, detail="classification")
        full = scan_range(1, 1000)
        assert fast.total == full.total

    def test_fast_labels_are_classification_only(self):
        fast = scan_range(1, 1000, detail="classification")
        allowed = {"prime", "semiprime", "composite", "invalid"}
        assert set(fast.counts.keys()) <= allowed

    def test_fast_counts_sum_to_total(self):
        fast = scan_range(1, 1000, detail="classification")
        assert sum(fast.counts.values()) == fast.total

    def test_fast_with_only_classification(self):
        fast = scan_range(1, 1000, only_classification="prime", detail="classification")
        assert set(fast.counts.keys()) == {"prime"}

    def test_fast_prime_count_matches_full(self):
        # Total primes found must be the same regardless of detail level
        fast = scan_range(1, 1000, only_classification="prime", detail="classification")
        full = scan_range(1, 1000, only_classification="prime")
        assert fast.total == full.total


class TestClassificationValidation:
    def test_valid_classifications_are_known(self):
        assert VALID_CLASSIFICATIONS == {"prime", "semiprime", "composite", "invalid"}

    def test_invalid_raises_value_error(self):
        with pytest.raises(ValueError, match="unknown classification"):
            scan_range(1, 100, only_classification="banana")

    def test_error_message_names_the_bad_value(self):
        with pytest.raises(ValueError, match="banana"):
            scan_range(1, 100, only_classification="banana")

    def test_all_valid_classifications_accepted(self):
        for cls in VALID_CLASSIFICATIONS:
            result = scan_range(1, 50, only_classification=cls)
            assert isinstance(result.total, int)

    def test_case_insensitive_valid_input(self):
        # Library accepts mixed case — same as the CLI
        result = scan_range(1, 50, only_classification="Prime")
        assert isinstance(result.total, int)
