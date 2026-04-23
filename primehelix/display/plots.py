from __future__ import annotations

from collections import Counter
import math

def _require_matplotlib():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except ImportError:
        raise ImportError(
            "matplotlib is required for plotting. "
            "Install it with: pip install 'primehelix[plot]'"
        )


def _classification_prefix(label: str) -> str:
    if not label:
        return "unknown"
    return label.split("|", 1)[0].strip().lower()


def _classification_rank(label: str) -> int:
    order = {
        "composite": 0,
        "prime": 1,
        "semiprime": 2,
        "invalid": 3,
    }
    return order.get(_classification_prefix(label), 99)


def _sorted_items(counts: Counter, limit: int):
    items = list(counts.items())

    if not items:
        return []

    items.sort(
        key=lambda kv: (
            _classification_rank(kv[0]),
            -kv[1],
            kv[0].lower(),
        )
    )
    return items[:limit]


def save_structure_plot(
    counts: Counter,
    output_path: str,
    title: str = "Structure Distribution",
    limit: int = 20,
):
    plt = _require_matplotlib()
    items = _sorted_items(counts, limit)

    if not items:
        raise ValueError("No structure counts available to plot.")

    labels = [label for label, _ in items]
    values = [count for _, count in items]

    fig_height = max(4, min(14, 0.6 * len(labels) + 1.5))
    fig_width = 12

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    y = list(range(len(labels)))
    ax.barh(y, values)

    ax.set_title(title)
    ax.set_xlabel("Count")
    ax.set_ylabel("Structure")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()

    max_val = max(values) if values else 1
    pad = max(0.02 * max_val, 0.1)
    for yi, val in zip(y, values):
        ax.text(val + pad, yi, str(val), va="center")

    ax.set_xlim(0, max_val * 1.15 if max_val > 0 else 1)

    fig.tight_layout()
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def save_compare_plot(
    counts_a: Counter,
    counts_b: Counter,
    label_a: str,
    label_b: str,
    output_path: str,
    title: str = "Range Comparison",
    limit: int = 20,
    top_delta: int | None = None,
):
    plt = _require_matplotlib()
    labels = set(counts_a.keys()) | set(counts_b.keys())
    items = [(label, counts_a.get(label, 0), counts_b.get(label, 0)) for label in labels]

    if not items:
        raise ValueError("No structure counts available to compare.")

    def row_delta(row):
        _, a, b = row
        return abs(b - a)

    if top_delta:
        items.sort(
            key=lambda row: (
                -row_delta(row),
                row[0].lower(),
            )
        )
        items = items[:top_delta]
    else:
        items.sort(
            key=lambda row: (
                _classification_rank(row[0]),
                -(row[1] + row[2]),
                row[0].lower(),
            )
        )
        items = items[:limit]

    plot_labels = [row[0] for row in items]
    values_a = [row[1] for row in items]
    values_b = [row[2] for row in items]

    fig_height = max(4, min(16, 0.65 * len(plot_labels) + 1.5))
    fig_width = 13

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    y = list(range(len(plot_labels)))
    bar_h = 0.38

    ax.barh([yy - bar_h / 2 for yy in y], values_a, height=bar_h, label=label_a)
    ax.barh([yy + bar_h / 2 for yy in y], values_b, height=bar_h, label=label_b)

    ax.set_title(title)
    ax.set_xlabel("Count")
    ax.set_ylabel("Structure")
    ax.set_yticks(y)
    ax.set_yticklabels(plot_labels)
    ax.invert_yaxis()
    ax.legend()

    max_val = max(values_a + values_b) if (values_a or values_b) else 1
    pad = max(0.02 * max_val, 0.1)

    for yi, val in zip([yy - bar_h / 2 for yy in y], values_a):
        ax.text(val + pad, yi, str(val), va="center")
    for yi, val in zip([yy + bar_h / 2 for yy in y], values_b):
        ax.text(val + pad, yi, str(val), va="center")

    ax.set_xlim(0, max_val * 1.2 if max_val > 0 else 1)

    fig.tight_layout()
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def save_structure_time_series_plot(
    series_map: dict[str, list[float]],
    window_labels: list[str],
    output_path: str,
    title: str = "Structure Time Series",
    ylabel: str = "Percent",
):
    plt = _require_matplotlib()
    if not series_map:
        raise ValueError("No time series data available to plot.")
    if not window_labels:
        raise ValueError("No window labels available to plot.")

    num_points = len(window_labels)
    fig_width = max(10, min(18, 8 + num_points * 0.25))
    fig_height = 7

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    x = list(range(num_points))

    for label, values in series_map.items():
        ax.plot(x, values, marker="o", linewidth=2, markersize=4, label=label)

    ax.set_title(title)
    ax.set_xlabel("Window")
    ax.set_ylabel(ylabel)

    # Reduce x-label clutter for many windows
    if num_points <= 12:
        tick_positions = x
    else:
        step = max(1, math.ceil(num_points / 12))
        tick_positions = x[::step]
        if tick_positions[-1] != x[-1]:
            tick_positions.append(x[-1])

    tick_labels = [window_labels[i] for i in tick_positions]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=30, ha="right")

    ax.legend(loc="best")
    ax.grid(True, axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(fig)
