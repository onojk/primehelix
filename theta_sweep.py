from lopsided import summarize_lopsided

# Small sanity dataset for now
semiprimes = [
    (91, 7, 13),
    (202, 2, 101),
    (15, 3, 5),
    (77, 7, 11),
    (26, 2, 13),
    (221, 13, 17),
    (403, 13, 31),
    (899, 29, 31),
    (106, 2, 53),
    (1189, 29, 41),
]

thetas = [0.20, 0.25, 0.30, 0.35, 0.40]

print("\n=== THETA SWEEP ===")
print(f"{'theta':>8} {'total':>8} {'lopsided':>10} {'balanced':>10} {'ratio':>10}")
print("-" * 52)

for theta in thetas:
    s = summarize_lopsided(semiprimes, theta=theta)
    print(
        f"{s['theta']:>8.2f} "
        f"{s['total']:>8} "
        f"{s['lopsided']:>10} "
        f"{s['balanced']:>10} "
        f"{s['lopsided_ratio']:>10.4f}"
    )
