#!/usr/bin/env python3
"""
factor_ratio_spectral.py

Tests whether the factor ratio p/q in lopsided semiprimes carries spectral
structure beyond what random prime pairs produce.

Three spectral functions:
  C_ratio(t)       = |Σ exp(i·t·log(p/q))|  / N   [actual semiprimes]
  C_random_ratio(t)= |Σ exp(i·t·log(r1/r2))| / N   [shuffled prime pairs]
  C_theta(t)       = |Σ exp(i·t·θ)|          / N   [structural coord θ=log(p)/log(n)]

Also: autocorrelation and power spectrum of the log(p/q) sequence.
"""

import array
import math
import random
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── Config ────────────────────────────────────────────────────────────────────
N_MAX      = 500_000
N_T        = 5000
T_MAX      = 200.0
CHUNK      = 50              # t-values per memory batch
SIG_RATIO  = 1.5             # flag when C_ratio/C_rand exceeds this

DARK_BG   = "#0a0a12"
PANEL_BG  = "#12121e"
SPINE_COL = "#2a2a40"
TEXT_COL  = "#ccccee"

# ── SPF sieve ─────────────────────────────────────────────────────────────────

def spf_sieve(limit: int):
    spf = array.array("i", range(limit + 1))
    i = 2
    while i * i <= limit:
        if spf[i] == i:
            for j in range(i * i, limit + 1, i):
                if spf[j] == j:
                    spf[j] = i
        i += 1
    return spf

# ── Spectral function ─────────────────────────────────────────────────────────

def compute_C(vals: np.ndarray, t_arr: np.ndarray) -> np.ndarray:
    """C(t) = |Σ exp(i·t·v)| / N, computed in CHUNK-sized t batches."""
    N = len(vals)
    C = np.empty(len(t_arr))
    for i in range(0, len(t_arr), CHUNK):
        t_c = t_arr[i : i + CHUNK]            # (chunk,)
        ph  = np.outer(t_c, vals)             # (chunk, N)
        cs  = np.cos(ph).sum(axis=1)
        ss  = np.sin(ph).sum(axis=1)
        C[i : i + CHUNK] = np.sqrt(cs * cs + ss * ss)
    return C / N

# ── Signal analysis ───────────────────────────────────────────────────────────

def power_spectrum(x: np.ndarray):
    """One-sided periodogram of real sequence x."""
    xc    = x - x.mean()
    fft_x = np.fft.rfft(xc)
    psd   = (np.abs(fft_x) ** 2) / len(x)
    freqs = np.fft.rfftfreq(len(x))
    return freqs, psd


def top_peaks(freqs: np.ndarray, psd: np.ndarray,
              n_peaks: int = 10, min_sep: int = 100) -> list:
    """Return (freq, psd) for top n_peaks local maxima, min_sep bins apart."""
    candidates = [
        k for k in range(1, len(psd) - 1)
        if psd[k] > psd[k - 1] and psd[k] > psd[k + 1]
    ]
    candidates.sort(key=lambda k: -psd[k])
    selected = []
    for c in candidates:
        if all(abs(c - s) >= min_sep for s in selected):
            selected.append(c)
        if len(selected) == n_peaks:
            break
    return [(float(freqs[k]), float(psd[k])) for k in selected]


def autocorrelation(x: np.ndarray, max_lag: int) -> np.ndarray:
    """Normalised linear autocorrelation via zero-padded FFT. R(0)=1."""
    xc    = x - x.mean()
    N     = len(xc)
    fft_x = np.fft.rfft(xc, n=2 * N)
    acorr = np.fft.irfft(np.abs(fft_x) ** 2)[:max_lag]
    return acorr / max(acorr[0], 1e-30)

# ── Styling ───────────────────────────────────────────────────────────────────

def _style(ax, title: str = "", xlabel: str = "", ylabel: str = ""):
    ax.set_facecolor(PANEL_BG)
    ax.tick_params(colors=TEXT_COL, labelsize=8)
    for s in ax.spines.values():
        s.set_color(SPINE_COL)
    if title:   ax.set_title(title,  color=TEXT_COL, fontsize=9, pad=6)
    if xlabel:  ax.set_xlabel(xlabel, color=TEXT_COL, fontsize=9)
    if ylabel:  ax.set_ylabel(ylabel, color=TEXT_COL, fontsize=9)
    ax.grid(color=SPINE_COL, lw=0.4, alpha=0.55)

# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 65)
    print("factor_ratio_spectral.py")
    print("=" * 65)

    # ── Build sieve and collect lopsided semiprimes ───────────────────────────
    print("\nBuilding SPF sieve …", flush=True)
    spf = spf_sieve(N_MAX)

    print("Collecting lopsided semiprimes (bit_gap > 8) …", flush=True)
    ns_l, ps_l, qs_l = [], [], []
    for n in range(4, N_MAX + 1):
        p = spf[n]
        if p == n:
            continue
        q = n // p
        if spf[q] != q:
            continue
        if q.bit_length() - p.bit_length() > 8:
            ns_l.append(n)
            ps_l.append(p)
            qs_l.append(q)

    N = len(ns_l)
    ps      = np.array(ps_l, dtype=np.float64)
    qs      = np.array(qs_l, dtype=np.float64)
    ns_arr  = np.array(ns_l, dtype=np.float64)

    log_pq   = np.log(ps / qs)                  # all negative; range ~ [-11.7, -5.5]
    theta    = np.log(ps) / np.log(ns_arr)       # ∈ (0, 0.5)

    from collections import Counter
    p_freq = Counter(ps_l)

    print(f"  {N:,} lopsided semiprimes")
    print(f"  Distinct p values: {len(p_freq)}")
    print(f"  log(p/q) range:  [{log_pq.min():.4f}, {log_pq.max():.4f}]")
    print(f"  theta range:     [{theta.min():.5f}, {theta.max():.5f}]")

    # ── Random control: shuffle p and q independently ─────────────────────────
    rng      = np.random.default_rng(42)
    ps_shuf  = rng.permutation(ps)
    qs_shuf  = rng.permutation(qs)
    log_rand = np.log(ps_shuf / qs_shuf)

    # ── t array ───────────────────────────────────────────────────────────────
    t_arr = np.linspace(0.0, T_MAX, N_T)

    # ── Spectral functions ────────────────────────────────────────────────────
    t0 = time.time()
    print("\nComputing spectral functions …")

    print("  C_ratio       …", end=" ", flush=True); t1 = time.time()
    C_ratio = compute_C(log_pq,   t_arr);  print(f"{time.time()-t1:.1f}s")

    print("  C_random_ratio…", end=" ", flush=True); t1 = time.time()
    C_rand  = compute_C(log_rand, t_arr);  print(f"{time.time()-t1:.1f}s")

    print("  C_theta       …", end=" ", flush=True); t1 = time.time()
    C_theta = compute_C(theta,    t_arr);  print(f"{time.time()-t1:.1f}s")

    print(f"  Total: {time.time()-t0:.1f}s")

    # ── Ratio C_ratio / C_rand ────────────────────────────────────────────────
    noise_floor = 3.0 / math.sqrt(N)
    with np.errstate(divide="ignore", invalid="ignore"):
        quotient = np.where(C_rand > noise_floor, C_ratio / C_rand, np.nan)

    signal_mask = (quotient > SIG_RATIO) & ~np.isnan(quotient)
    sig_t  = t_arr[signal_mask]
    sig_q  = quotient[signal_mask]

    # ── Autocorrelation & power spectrum ──────────────────────────────────────
    max_lag      = min(N // 2, 30_000)
    acorr        = autocorrelation(log_pq, max_lag)
    freqs, psd   = power_spectrum(log_pq)
    peaks        = top_peaks(freqs, psd, n_peaks=10, min_sep=100)

    # ── Print: signal t-values ────────────────────────────────────────────────
    print(f"\nNoise floor (3/√N): {noise_floor:.5f}")
    print(f"Signal threshold:   C_ratio/C_rand > {SIG_RATIO:.1f}")
    valid_q   = quotient[~np.isnan(quotient)]
    mean_q_gt10 = float(np.nanmean(quotient[t_arr > 10]))

    print(f"Mean ratio (t>10):  {mean_q_gt10:.4f}")
    print(f"Signal t-values:    {len(sig_t)}")

    if len(sig_t) == 0:
        print("\n  → None. C_ratio never significantly exceeds C_random_ratio.")
    else:
        top_order = np.argsort(-sig_q)[:20]
        print(f"\n  Top flagged t-values (up to 20):")
        print(f"  {'t':>10}  {'ratio':>8}  {'C_ratio':>10}  {'C_rand':>10}")
        print("  " + "─" * 44)
        for i in top_order:
            ti = sig_t[i]
            qi = sig_q[i]
            cr = C_ratio[signal_mask][i]
            cd = C_rand[signal_mask][i]
            print(f"  {ti:>10.4f}  {qi:>8.3f}  {cr:>10.6f}  {cd:>10.6f}")

    # ── Print: power spectrum peaks ───────────────────────────────────────────
    print(f"\nTop-10 power spectrum peaks of log(p/q) sequence:")
    print(f"  {'freq (cyc/sp)':>14}  {'period (sp)':>13}  {'PSD':>12}")
    print("  " + "─" * 44)
    for f_pk, psd_pk in peaks:
        period = 1.0 / f_pk if f_pk > 0 else float("inf")
        print(f"  {f_pk:>14.7f}  {period:>13.1f}  {psd_pk:>12.5f}")

    # ── Print: autocorrelation decay ──────────────────────────────────────────
    # Find lag where |acorr| first drops below 1.96/sqrt(N) (95% CI)
    ci95 = 1.96 / math.sqrt(N)
    decay_lag = next((i for i, v in enumerate(acorr) if abs(v) < ci95), len(acorr))
    print(f"\nAutocorrelation:")
    print(f"  R(1)  = {acorr[1]:.6f}  (lag-1 correlation)")
    print(f"  R(10) = {acorr[10]:.6f}")
    print(f"  R(100)= {acorr[100]:.6f}")
    print(f"  Drops below 95% CI ({ci95:.5f}) at lag {decay_lag}")

    # ── Verdict ───────────────────────────────────────────────────────────────
    print()
    if len(sig_t) > 20 and mean_q_gt10 > 1.15:
        verdict = "STRUCTURE DETECTED — factor ratios show sustained coherence above random pairs"
    elif len(sig_t) > 0:
        verdict = "MARGINAL — isolated t-values exceed threshold; no sustained excess over random"
    else:
        verdict = "NO STRUCTURE — factor ratios are spectrally indistinguishable from random prime pairs"
    print(f"KEY QUESTION: {verdict}")

    # ═════════════════════════════════════════════════════════════════════════
    # PLOT 1 — spectral functions
    # ═════════════════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(
        2, 1, figsize=(14, 9),
        gridspec_kw={"height_ratios": [2.2, 1]},
    )
    fig.patch.set_facecolor(DARK_BG)

    # Top panel: three C(t) curves
    ax = axes[0]
    _style(ax,
           title=(r"$C(t) = |\Sigma\,e^{i\,t\,x}|\,/\,N$  for three spectral domains  "
                  f"(N={N:,} lopsided semiprimes, bit_gap>8)"),
           ylabel="C(t) / N")
    ax.plot(t_arr, C_ratio, color="#3ecf96", lw=0.9, alpha=0.95,
            label=r"$C_{\rm ratio}$  — log(p/q)  [semiprimes]")
    ax.plot(t_arr, C_rand,  color="#f07050", lw=0.9, alpha=0.85,
            label=r"$C_{\rm rand}$   — log(r₁/r₂) [shuffled pairs]")
    ax.plot(t_arr, C_theta, color="#9f8ff8", lw=0.9, alpha=0.85,
            label=r"$C_{\theta}$    — θ = log(p)/log(n)")
    if len(sig_t) > 0:
        ax.scatter(sig_t, C_ratio[signal_mask], s=12, color="#ffd080",
                   zorder=5, label=f"Signal (ratio>{SIG_RATIO}×, n={len(sig_t)})")
    ax.legend(facecolor=PANEL_BG, edgecolor=SPINE_COL, labelcolor=TEXT_COL,
              fontsize=9, loc="upper right")
    ax.set_xlim(0, T_MAX)
    ax.set_ylim(bottom=0)

    # Bottom panel: ratio quotient
    ax2 = axes[1]
    _style(ax2,
           title=r"$C_{\rm ratio}(t)\,/\,C_{\rm rand}(t)$  "
                 f"(masked where C_rand < noise floor {noise_floor:.4f})",
           xlabel="t",
           ylabel=r"$C_{\rm ratio}\,/\,C_{\rm rand}$")
    valid_t    = t_arr[~np.isnan(quotient)]
    valid_quot = quotient[~np.isnan(quotient)]
    ax2.plot(valid_t, valid_quot, color="#ffd080", lw=0.8, alpha=0.85)
    ax2.axhline(SIG_RATIO, color="#f07050", lw=1.1, ls="--",
                label=f"threshold {SIG_RATIO}×")
    ax2.axhline(1.0, color="#aaaaaa", lw=0.7, ls=":")
    if len(sig_t) > 0:
        ax2.scatter(sig_t, sig_q, s=10, color="#f07050", zorder=5)
    ax2.set_xlim(0, T_MAX)
    y_hi = min(float(np.nanmax(quotient)) * 1.05, 6.0)
    ax2.set_ylim(0, y_hi)
    ax2.legend(facecolor=PANEL_BG, edgecolor=SPINE_COL, labelcolor=TEXT_COL, fontsize=8)

    plt.tight_layout(pad=0.5)
    f1 = "factor_ratio_spectral.png"
    plt.savefig(f1, dpi=130, facecolor=DARK_BG)
    plt.close()
    print(f"\nSaved: {f1}")

    # ═════════════════════════════════════════════════════════════════════════
    # PLOT 2 — autocorrelation + power spectrum
    # ═════════════════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.patch.set_facecolor(DARK_BG)

    # Left: autocorrelation
    ax = axes[0]
    show_lags = min(2000, len(acorr))
    _style(ax,
           title="Autocorrelation of log(p/q)  (ordered by n)\n"
                 "Red dotted = 95% white-noise CI",
           xlabel="Lag  (semiprimes ordered by n)",
           ylabel="R(lag)")
    ax.plot(np.arange(show_lags), acorr[:show_lags],
            color="#3ecf96", lw=0.8, alpha=0.9)
    ax.axhline(0,     color="#888888", lw=0.6, ls="--")
    ax.axhline( ci95, color="#f07050", lw=0.7, ls=":", alpha=0.8)
    ax.axhline(-ci95, color="#f07050", lw=0.7, ls=":", alpha=0.8)
    ax.set_xlim(0, show_lags)

    # Right: power spectrum (low-frequency focus)
    ax = axes[1]
    freq_hi  = 0.05    # show up to 5% Nyquist — long-period structure
    f_mask   = (freqs > 0) & (freqs <= freq_hi)
    _style(ax,
           title="Power spectrum of log(p/q)  (low-frequency)\n"
                 "Yellow dashed = top-10 peak frequencies",
           xlabel="Frequency  (cycles per semiprime)",
           ylabel="PSD  (log scale)")
    ax.semilogy(freqs[f_mask], psd[f_mask], color="#9f8ff8", lw=0.8)
    for f_pk, _ in peaks:
        if 0 < f_pk <= freq_hi:
            ax.axvline(f_pk, color="#ffd080", lw=0.9, ls="--", alpha=0.75)
    ax.set_xlim(0, freq_hi)

    plt.tight_layout(pad=0.5)
    f2 = "factor_ratio_autocorr.png"
    plt.savefig(f2, dpi=130, facecolor=DARK_BG)
    plt.close()
    print(f"Saved: {f2}")

    # ── Summary table ─────────────────────────────────────────────────────────
    print()
    print("─" * 65)
    print("SUMMARY TABLE")
    print("─" * 65)
    print(f"  N lopsided semiprimes (bit_gap>8) : {N:,}")
    print(f"  t range                           : [0, {T_MAX}]  ({N_T} points)")
    print(f"  Noise floor (3/√N)                : {noise_floor:.5f}")
    print(f"  Mean C_ratio/C_rand (t>10)        : {mean_q_gt10:.4f}")
    print(f"  Signal t-values (>{SIG_RATIO}×)         : {len(sig_t)}")
    print(f"  Autocorr lag-1                    : {acorr[1]:.6f}")
    print(f"  Autocorr 95%CI decay lag          : {decay_lag}")
    print(f"  Top PSD peak frequency            : {peaks[0][0]:.7f} cyc/sp  "
          f"(period={1/peaks[0][0]:.1f} sp)" if peaks else "  No peaks found")
    print()
    print(f"  VERDICT: {verdict}")
    print("─" * 65)


if __name__ == "__main__":
    main()
