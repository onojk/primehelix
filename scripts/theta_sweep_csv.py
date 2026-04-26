import csv
import sys
from pathlib import Path

from lopsided import summarize_lopsided


THETAS = [0.20, 0.25, 0.30, 0.35, 0.40]


def find_col(fieldnames, options):
    lower_map = {name.lower(): name for name in fieldnames}
    for opt in options:
        if opt.lower() in lower_map:
            return lower_map[opt.lower()]
    return None


def load_semiprimes(csv_path):
    rows = []

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames or []

        n_col = find_col(fields, ["n", "number", "value"])
        p_col = find_col(fields, ["p", "factor_p", "small_factor", "min_factor"])
        q_col = find_col(fields, ["q", "factor_q", "large_factor", "max_factor"])

        if not n_col or not p_col or not q_col:
            raise ValueError(
                f"Could not find n/p/q columns in {csv_path}\n"
                f"Found columns: {fields}"
            )

        for row in reader:
            n = int(row[n_col])
            p = int(row[p_col])
            q = int(row[q_col])
            rows.append((n, p, q))

    return rows


def main():
    if len(sys.argv) < 2:
        print("Usage:")
        print("  python theta_sweep_csv.py path/to/semiprime_labels.csv")
        sys.exit(1)

    csv_path = Path(sys.argv[1])
    semiprimes = load_semiprimes(csv_path)

    print(f"\n=== THETA SWEEP: {csv_path.name} ===")
    print(f"{'theta':>8} {'total':>10} {'lopsided':>12} {'balanced':>12} {'ratio':>10}")
    print("-" * 60)

    for theta in THETAS:
        s = summarize_lopsided(semiprimes, theta=theta)
        print(
            f"{s['theta']:>8.2f} "
            f"{s['total']:>10} "
            f"{s['lopsided']:>12} "
            f"{s['balanced']:>12} "
            f"{s['lopsided_ratio']:>10.4f}"
        )


if __name__ == "__main__":
    main()
