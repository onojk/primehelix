"""Core factoring and primality tests."""
import pytest
from primehelix.core.primes import is_prime
from primehelix.core.rho import pollard_rho
from primehelix.core.pm1 import pollard_pm1, williams_pp1
from primehelix.core.ecm import ecm
from primehelix.core.qs import quadratic_sieve
from primehelix.core.factor import factor, classify


# --- Primality ---

def test_small_primes():
    assert all(is_prime(p) for p in [2, 3, 5, 7, 11, 13, 97, 997])

def test_composites_not_prime():
    assert not any(is_prime(n) for n in [1, 4, 9, 15, 100, 561])

def test_large_prime():
    assert is_prime(2**31 - 1)  # Mersenne prime

def test_large_composite():
    assert not is_prime(2**31 - 1 + 2)


# --- Rho ---

def test_rho_simple():
    assert pollard_rho(15) in (3, 5)

def test_rho_semiprime():
    n = 10000019 * 11000027
    f = pollard_rho(n, timeout_ms=5000)
    assert f is not None and n % f == 0 and 1 < f < n


# --- p-1 ---

def test_pm1_smooth():
    # 2^31 - 1 = 2147483647 (prime), but (2^31-1)*2 has factor 2
    assert pollard_pm1(2 * 97) in (2, 97)

def test_pm1_finds_smooth_factor():
    # p=17, p-1=16=2^4 (smooth to B1=20)
    # q=1013, q-1=1012=4*11*23 (23>20, not smooth — so gcd gives 17 not n)
    n = 17 * 1013
    f = pollard_pm1(n, B1=20)
    assert f is not None and n % f == 0 and 1 < f < n


# --- ECM ---

def test_ecm_tiny():
    f = ecm(8051, B1=1000, curves=10)
    assert f is not None and 8051 % f == 0

def test_ecm_medium():
    n = 104729314187  # = 10007 * 10466341
    f = ecm(n, B1=10000, curves=50, timeout_ms=10000)
    assert f is not None and n % f == 0 and 1 < f < n


# --- QS ---

def test_qs_small():
    # 8633 = 89 * 97
    f = quadratic_sieve(8633)
    assert f in (89, 97)

def test_qs_medium():
    # Two 4-digit primes, both > trial division limit
    n = 1009 * 1013
    f = quadratic_sieve(n)
    assert f is not None and n % f == 0 and 1 < f < n


# --- Full pipeline ---

def test_factor_prime():
    r = factor(97)
    assert r.factors == {97: 1}

def test_factor_semiprime():
    n = 91  # 7 * 13
    r = factor(n)
    assert r.factors.get(7) == 1 and r.factors.get(13) == 1

def test_factor_power():
    r = factor(2**10)
    assert r.factors == {2: 10}

def test_classify_prime():
    cls, r = classify(2147483647)
    assert cls == "prime"

def test_classify_semiprime():
    cls, r = classify(91)
    assert cls == "semiprime"

def test_classify_composite():
    cls, r = classify(2 * 3 * 5 * 7)
    assert cls == "composite"
