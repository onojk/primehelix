from lopsided import is_lopsided

def test_examples():
    assert not is_lopsided(91, 7, 13)   # 7^3=343 > 13 → not lopsided
    assert is_lopsided(202, 2, 101)     # 2^3=8 <= 101 → lopsided
