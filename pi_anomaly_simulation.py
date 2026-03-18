"""
================================================================================
PI ANOMALY IN OFF-DIAGONAL AUBRY-ANDRÉ SYSTEMS
Simulation Code Package v1.0

Author: Bradley Joseph Downey Jr, Independent Researcher (2026)
Preprint: [arXiv link once posted]
Patent: USPTO Provisional Application 63/XXX,XXX

Description:
    Compares five irrational modulation parameters (pi, phi, sqrt2, e, generic)
    in the off-diagonal Aubry-André model across system sizes N=200 to N=3000.
    Measures localization length, level-spacing r-statistic, and spectral
    fractal dimension at the localization transition lambda/t = 2.0.

Requirements:
    Python 3.8+
    numpy
    scipy
    matplotlib

Install dependencies:
    pip install numpy scipy matplotlib

Usage:
    python pi_anomaly_simulation.py

    This will run the full simulation and produce:
    - pi_anomaly_results.png  (main figure)
    - pi_anomaly_data.csv     (raw numerical data)
    - pi_anomaly_summary.txt  (printed summary of key results)

Runtime estimate:
    N=200-1000: ~5 minutes on a standard laptop
    N=2000:     ~15 minutes
    N=3000:     ~25 minutes
    Full run:   ~50 minutes

To reproduce specific system sizes only, edit the N_SIZES list below.

Key result to reproduce:
    At lambda/t = 2.0, the pi/phi localization length ratio should be:
    N=200:  ~1.16
    N=500:  ~1.67
    N=1000: ~2.59
    N=2000: ~2.83
    N=3000: ~3.23
    with no plateau observed.
================================================================================
"""

import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import csv
import time
import sys

# ── Reproducibility ────────────────────────────────────────────────────────
RANDOM_SEED = 42
rng = np.random.default_rng(RANDOM_SEED)

# ── Irrational modulation parameters ──────────────────────────────────────
PI    = np.pi                        # [3; 7, 15, 1, 292, 1, 1, 1, 2, ...]
PHI   = (1 + np.sqrt(5)) / 2        # [1; 1, 1, 1, 1, ...] Golden Ratio
SQRT2 = np.sqrt(2)                   # [1; 2, 2, 2, 2, ...]
E     = np.e                         # [2; 1, 2, 1, 1, 4, 1, 1, ...]
GEN   = 1.2020569 / 10 * 3          # Generic irrational, no special structure

ALPHAS = {
    "pi":      PI,
    "phi":     PHI,
    "sqrt2":   SQRT2,
    "e":       E,
    "generic": GEN,
}

ALPHA_LABELS = {
    "pi":      "π",
    "phi":     "φ (Golden Ratio)",
    "sqrt2":   "√2",
    "e":       "e",
    "generic": "generic",
}

COLORS = {
    "pi":      "#00d4ff",
    "phi":     "#ff6b35",
    "sqrt2":   "#a8ff78",
    "e":       "#f7971e",
    "generic": "#cc44ff",
}

# ── Simulation parameters ──────────────────────────────────────────────────
# Edit N_SIZES to run only specific system sizes
N_SIZES  = [200, 500, 1000, 2000, 3000]

# Number of random phase realisations per (N, alpha, lambda) point
# Higher N_PHASE = more accurate but slower
# Paper used N_PHASE=5 for N<=2000 and N_PHASE=3 for N=3000
N_PHASE  = 5

# Modulation strengths to sweep
# 2.0 is the critical transition point (main result)
LAMBDAS  = [1.5, 2.0, 2.5, 3.0]

# Baseline hopping strength (energy scale)
T = 1.0

# ── Continued fraction tools (for reference/verification) ─────────────────
def get_continued_fraction(x, n_terms=10):
    """Return first n_terms partial quotients of x."""
    cf, xi = [], float(x)
    for _ in range(n_terms):
        a = int(xi)
        cf.append(a)
        remainder = xi - a
        if remainder < 1e-12:
            break
        xi = 1.0 / remainder
    return cf

# ── Hamiltonian construction ───────────────────────────────────────────────
def build_hamiltonian(N, lam, alpha, phi_offset, t=1.0):
    """
    Build the off-diagonal Aubry-André Hamiltonian.

    H = t * sum_{n=0}^{N-2} (1 + lam*cos(2*pi*alpha*n + phi)) * |n+1><n| + h.c.

    Parameters:
        N          : number of sites
        lam        : modulation amplitude (lambda/t ratio when t=1)
        alpha      : irrational modulation frequency
        phi_offset : phase offset (averaged over in the simulation)
        t          : baseline hopping strength

    Returns:
        H : (N, N) numpy array, real symmetric
    """
    n_idx    = np.arange(N - 1)
    hoppings = t * (1.0 + lam * np.cos(2 * PI * alpha * n_idx + phi_offset))
    H = np.zeros((N, N), dtype=float)
    H[n_idx, n_idx + 1] = hoppings
    H[n_idx + 1, n_idx] = hoppings
    return H

# ── Spectral metrics ───────────────────────────────────────────────────────
def localization_length(evecs, N, central_fraction=0.10):
    """
    Compute mean localization length for central eigenstates.

    Uses IPR (Inverse Participation Ratio) as localization proxy:
        IPR = sum_j |psi_j|^4
        xi  = 1 / IPR

    Central eigenstates (near E=0) are used because they are most
    sensitive to the localization transition in chiral systems.

    Parameters:
        evecs            : (N, N) eigenvector matrix from eigh()
        N                : system size
        central_fraction : fraction of eigenstates to average over

    Returns:
        mean localization length (float)
    """
    mid  = N // 2
    half = max(1, int(central_fraction * N / 2))
    selected = evecs[:, mid - half : mid + half]
    iprs     = np.sum(selected**4, axis=0)
    valid    = iprs[iprs > 0]
    return float(np.mean(1.0 / valid)) if len(valid) > 0 else 0.0

def r_statistic(evals):
    """
    Compute mean level-spacing r-statistic.

    r_n = min(s_n, s_{n+1}) / max(s_n, s_{n+1})
    where s_n = E_{n+1} - E_n are consecutive level spacings.

    Reference values:
        Poisson (localized):   <r> = 0.386
        GOE    (delocalized):  <r> = 0.530
        Values > 0.530 indicate anomalous level repulsion.

    Parameters:
        evals : 1D array of eigenvalues

    Returns:
        mean r value (float), or nan if insufficient levels
    """
    ev       = np.sort(evals)
    spacings = np.diff(ev)
    spacings = spacings[spacings > 1e-12]
    if len(spacings) < 2:
        return np.nan
    r = (np.minimum(spacings[:-1], spacings[1:]) /
         np.maximum(spacings[:-1], spacings[1:]))
    return float(np.mean(r))

def spectral_fractal_dimension(evals, n_boxes=(10, 20, 50, 100, 200, 500)):
    """
    Estimate box-counting fractal dimension of the eigenvalue spectrum.

    Counts occupied boxes at multiple resolutions and fits the
    log(count) vs log(resolution) slope.

    A Cantor-set-like spectrum has D < 1.
    A fully filled spectrum has D = 1.

    Parameters:
        evals   : 1D array of eigenvalues
        n_boxes : tuple of box counts to use for fitting

    Returns:
        fractal dimension D (float)
    """
    ev   = np.sort(evals)
    span = ev[-1] - ev[0]
    if span < 1e-10:
        return 0.0
    counts = []
    for nb in n_boxes:
        box_size = span / nb
        occupied = len(set(int((e - ev[0]) / box_size) for e in ev))
        counts.append(occupied)
    log_res    = np.log(np.array(n_boxes, dtype=float))
    log_counts = np.log(np.array(counts,  dtype=float))
    if len(set(log_counts)) < 2:
        return 1.0
    slope, _ = np.polyfit(log_res, log_counts, 1)
    return float(slope)

# ── Main simulation loop ───────────────────────────────────────────────────
def run_simulation(n_sizes=N_SIZES, n_phase=N_PHASE, lambdas=LAMBDAS):
    """
    Run the full simulation across all system sizes, parameters, and lambda values.

    Returns:
        results : dict with structure results[alpha_key][metric][N_index, lambda_index]
    """
    results = {
        name: {
            "loc_len":  np.zeros((len(n_sizes), len(lambdas))),
            "loc_std":  np.zeros((len(n_sizes), len(lambdas))),
            "r_stat":   np.zeros((len(n_sizes), len(lambdas))),
            "r_std":    np.zeros((len(n_sizes), len(lambdas))),
            "frac_dim": np.zeros((len(n_sizes), len(lambdas))),
            "frac_std": np.zeros((len(n_sizes), len(lambdas))),
        }
        for name in ALPHAS
    }

    total_runs = len(n_sizes) * len(ALPHAS) * len(lambdas) * n_phase
    completed  = 0
    t_start    = time.time()

    print("=" * 70)
    print("PI ANOMALY SIMULATION")
    print(f"System sizes:       {n_sizes}")
    print(f"Lambda values:      {lambdas}")
    print(f"Phase realisations: {n_phase}")
    print(f"Total eigensolves:  {total_runs}")
    print(f"Random seed:        {RANDOM_SEED}")
    print("=" * 70)

    print("\nContinued fractions of modulation parameters:")
    for name, alpha in ALPHAS.items():
        cf = get_continued_fraction(alpha, 8)
        print(f"  {ALPHA_LABELS[name]:20s}: {cf}")
    print()

    for ni, N in enumerate(n_sizes):
        print(f"\nN = {N}")
        print("-" * 40)

        for name, alpha in ALPHAS.items():
            t0 = time.time()

            for li, lam in enumerate(lambdas):
                lls, rrs, fds = [], [], []

                for _ in range(n_phase):
                    phi = rng.uniform(0, 2 * PI)
                    H   = build_hamiltonian(N, lam, alpha, phi, T)

                    # Full diagonalisation
                    evals, evecs = eigh(H)

                    lls.append(localization_length(evecs, N))
                    rrs.append(r_statistic(evals))
                    fds.append(spectral_fractal_dimension(evals))
                    completed += 1

                results[name]["loc_len"][ni, li]  = np.mean(lls)
                results[name]["loc_std"][ni, li]   = np.std(lls)
                results[name]["r_stat"][ni, li]    = np.nanmean(rrs)
                results[name]["r_std"][ni, li]     = np.nanstd(rrs)
                results[name]["frac_dim"][ni, li]  = np.mean(fds)
                results[name]["frac_std"][ni, li]  = np.std(fds)

            elapsed = time.time() - t0
            lam_idx = lambdas.index(2.0) if 2.0 in lambdas else 0
            print(f"  {ALPHA_LABELS[name]:20s}  "
                  f"ξ={results[name]['loc_len'][ni,lam_idx]:.2f}  "
                  f"r={results[name]['r_stat'][ni,lam_idx]:.4f}  "
                  f"D={results[name]['frac_dim'][ni,lam_idx]:.4f}  "
                  f"({elapsed:.1f}s)")
            sys.stdout.flush()

    print(f"\nTotal simulation time: {time.time()-t_start:.1f}s")
    return results

# ── Results printing ───────────────────────────────────────────────────────
def print_summary(results, n_sizes, lambdas):
    """Print key results to console and save to text file."""

    lam_idx = lambdas.index(2.0) if 2.0 in lambdas else 0
    lam_val = lambdas[lam_idx]

    lines = []
    lines.append("=" * 70)
    lines.append(f"KEY RESULTS AT λ/t = {lam_val}")
    lines.append("=" * 70)

    # Localization length table
    lines.append(f"\nLOCALIZATION LENGTH")
    header = f"{'α':>10}" + "".join(f"  {'N='+str(N):>8}" for N in n_sizes)
    lines.append(header)
    lines.append("-" * len(header))
    for name in ALPHAS:
        row = f"{ALPHA_LABELS[name]:>10}"
        for ni in range(len(n_sizes)):
            row += f"  {results[name]['loc_len'][ni, lam_idx]:>8.2f}"
        lines.append(row)

    # π/φ ratio
    lines.append(f"\nπ/φ LOCALIZATION LENGTH RATIO")
    for ni, N in enumerate(n_sizes):
        pi_val  = results["pi"]["loc_len"][ni, lam_idx]
        phi_val = results["phi"]["loc_len"][ni, lam_idx]
        ratio   = pi_val / phi_val if phi_val > 0 else 0
        adv     = (ratio - 1) * 100
        lines.append(f"  N={N:>5}: π={pi_val:.2f}  φ={phi_val:.2f}  "
                     f"ratio={ratio:.3f}  advantage={adv:+.1f}%")

    # Fractal dimension
    lines.append(f"\nSPECTRAL FRACTAL DIMENSION")
    for name in ALPHAS:
        row = f"  {ALPHA_LABELS[name]:>20}"
        for ni in range(len(n_sizes)):
            row += f"  {results[name]['frac_dim'][ni, lam_idx]:.4f}"
        lines.append(row)

    # r-statistic
    lines.append(f"\nr-STATISTIC  (Poisson=0.386, GOE=0.530)")
    for name in ALPHAS:
        row = f"  {ALPHA_LABELS[name]:>20}"
        for ni in range(len(n_sizes)):
            val  = results[name]['r_stat'][ni, lam_idx]
            flag = " *** SUPER-GOE ***" if val > 0.55 else ""
            row += f"  {val:.4f}{flag}"
        lines.append(row)

    # Verdict
    lines.append("\n" + "=" * 70)
    lines.append("VERDICT")
    lines.append("=" * 70)

    if len(n_sizes) >= 2:
        ratios = []
        for ni in range(len(n_sizes)):
            pi_val  = results["pi"]["loc_len"][ni, lam_idx]
            phi_val = results["phi"]["loc_len"][ni, lam_idx]
            ratios.append(pi_val / phi_val if phi_val > 0 else 0)

        if all(ratios[i] <= ratios[i+1] for i in range(len(ratios)-1)):
            lines.append("✓ SCALING CONFIRMED: π/φ ratio grows monotonically")
            lines.append(f"  Range: {ratios[0]:.3f} → {ratios[-1]:.3f} "
                         f"across N={n_sizes[0]}→{n_sizes[-1]}")
        else:
            lines.append("⚠ SCALING NON-MONOTONIC: check results carefully")

    pi_r = results["pi"]["r_stat"][-1, lam_idx]
    if pi_r > 0.55:
        lines.append(f"✓ r-STATISTIC ANOMALY CONFIRMED: "
                     f"π shows r={pi_r:.4f} > GOE ceiling at N={n_sizes[-1]}")

    output = "\n".join(lines)
    print("\n" + output)

    with open("pi_anomaly_summary.txt", "w") as f:
        f.write(output)
    print("\nSummary saved to: pi_anomaly_summary.txt")

# ── CSV export ─────────────────────────────────────────────────────────────
def save_csv(results, n_sizes, lambdas):
    """Save all results to CSV for external analysis."""
    with open("pi_anomaly_data.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "alpha", "alpha_value", "N", "lambda_over_t",
            "loc_len_mean", "loc_len_std",
            "r_stat_mean",  "r_stat_std",
            "frac_dim_mean","frac_dim_std",
        ])
        for name, alpha in ALPHAS.items():
            for ni, N in enumerate(n_sizes):
                for li, lam in enumerate(lambdas):
                    writer.writerow([
                        ALPHA_LABELS[name], f"{alpha:.10f}", N, lam,
                        f"{results[name]['loc_len'][ni,li]:.6f}",
                        f"{results[name]['loc_std'][ni,li]:.6f}",
                        f"{results[name]['r_stat'][ni,li]:.6f}",
                        f"{results[name]['r_std'][ni,li]:.6f}",
                        f"{results[name]['frac_dim'][ni,li]:.6f}",
                        f"{results[name]['frac_std'][ni,li]:.6f}",
                    ])
    print("Data saved to: pi_anomaly_data.csv")

# ── Plotting ───────────────────────────────────────────────────────────────
def plot_results(results, n_sizes, lambdas):
    """Generate the main results figure."""

    lam_idx  = lambdas.index(2.0) if 2.0 in lambdas else 0
    lam_full = np.array(lambdas)
    N_arr    = np.array(n_sizes)

    fig = plt.figure(figsize=(18, 14), facecolor="#0a0a1a")
    gs  = gridspec.GridSpec(3, 2, figure=fig, hspace=0.50, wspace=0.35)
    axes = [fig.add_subplot(gs[r, c]) for r in range(3) for c in range(2)]

    for ax in axes:
        ax.set_facecolor("#0f0f2a")
        ax.tick_params(colors="white", labelsize=9)
        for sp in ax.spines.values():
            sp.set_edgecolor("#333355")

    def tl(ax, title, xlabel, ylabel):
        ax.set_title(title, color="white", fontsize=10, pad=6)
        ax.set_xlabel(xlabel, color="#aaaacc", fontsize=9)
        ax.set_ylabel(ylabel, color="#aaaacc", fontsize=9)

    LW = {n: 2.5 if n == "pi" else 1.6 for n in ALPHAS}

    # Plot 1: Localization length vs lambda at N=largest available
    ni_max = len(n_sizes) - 1
    N_max  = n_sizes[ni_max]
    tl(axes[0], f"Localization Length vs λ/t  (N={N_max})",
       "λ/t", "Localization length ξ")
    for name in ALPHAS:
        axes[0].plot(lam_full, results[name]["loc_len"][ni_max],
                     color=COLORS[name], lw=LW[name],
                     label=ALPHA_LABELS[name], alpha=0.85)
    axes[0].axvline(2.0, color="#ff4466", lw=1.2, ls="--", alpha=0.7, label="λ/t=2")
    axes[0].legend(fontsize=8, facecolor="#0a0a1a",
                   edgecolor="#333355", labelcolor="white")

    # Plot 2: Localization length vs N at lambda=2.0
    tl(axes[1], "π/φ Scaling: Localization Length vs N  (λ/t=2.0)\nKey result",
       "System size N", "Localization length ξ")
    for name in ALPHAS:
        vals = results[name]["loc_len"][:, lam_idx]
        errs = results[name]["loc_std"][:, lam_idx]
        axes[1].errorbar(N_arr, vals, yerr=errs,
                         color=COLORS[name], lw=LW[name],
                         marker="o", ms=6, capsize=3,
                         label=ALPHA_LABELS[name], alpha=0.85)
    axes[1].legend(fontsize=8, facecolor="#0a0a1a",
                   edgecolor="#333355", labelcolor="white")

    # Plot 3: π/φ ratio vs N
    tl(axes[2], "π/φ Localization Length Ratio vs N\n(super-linear growth = key claim)",
       "N", "π/φ ratio")
    ratios = [results["pi"]["loc_len"][ni, lam_idx] /
              results["phi"]["loc_len"][ni, lam_idx]
              for ni in range(len(n_sizes))]
    axes[2].plot(N_arr, ratios, color=COLORS["pi"], lw=2.5,
                 marker="o", ms=8, label="π/φ ratio")
    axes[2].axhline(1.0, color="white", lw=0.8, ls=":", alpha=0.4)
    for ni, (N, r) in enumerate(zip(n_sizes, ratios)):
        axes[2].annotate(f"{r:.2f}", (N, r),
                         textcoords="offset points", xytext=(0, 10),
                         color="white", fontsize=8, ha="center")
    axes[2].legend(fontsize=8, facecolor="#0a0a1a",
                   edgecolor="#333355", labelcolor="white")

    # Plot 4: Fractal dimension vs N
    tl(axes[3], "Spectral Fractal Dimension vs N  (λ/t=2.0)",
       "N", "Box-counting dimension D")
    for name in ALPHAS:
        vals = results[name]["frac_dim"][:, lam_idx]
        axes[3].plot(N_arr, vals, color=COLORS[name], lw=LW[name],
                     marker="s", ms=5, label=ALPHA_LABELS[name], alpha=0.85)
    axes[3].legend(fontsize=8, facecolor="#0a0a1a",
                   edgecolor="#333355", labelcolor="white")

    # Plot 5: r-statistic vs N
    tl(axes[4], "Level-Spacing r-Statistic vs N  (λ/t=2.0)\n(GOE ceiling = 0.530)",
       "N", "⟨r⟩")
    for name in ALPHAS:
        vals = results[name]["r_stat"][:, lam_idx]
        axes[4].plot(N_arr, vals, color=COLORS[name], lw=LW[name],
                     marker="^", ms=5, label=ALPHA_LABELS[name], alpha=0.85)
    axes[4].axhline(0.386, color="white",  lw=1.0, ls=":", alpha=0.5, label="Poisson")
    axes[4].axhline(0.530, color="yellow", lw=1.0, ls=":", alpha=0.5, label="GOE")
    axes[4].legend(fontsize=8, facecolor="#0a0a1a",
                   edgecolor="#333355", labelcolor="white")

    # Plot 6: π deviation from φ across all three metrics
    tl(axes[5], "π Relative Advantage over φ\nacross all three metrics",
       "N", "Relative advantage (%)")
    for metric, label, color in [
        ("loc_len",  "Loc length",   "#00d4ff"),
        ("frac_dim", "Fractal dim",  "#a8ff78"),
        ("r_stat",   "r-statistic",  "#f7971e"),
    ]:
        adv = [(results["pi"][metric][ni, lam_idx] -
                results["phi"][metric][ni, lam_idx]) /
               (results["phi"][metric][ni, lam_idx] + 1e-10) * 100
               for ni in range(len(n_sizes))]
        axes[5].plot(N_arr, adv, color=color, lw=1.8,
                     marker="o", ms=5, label=label, alpha=0.85)
    axes[5].axhline(0, color="white", lw=0.8, ls=":", alpha=0.4)
    axes[5].legend(fontsize=8, facecolor="#0a0a1a",
                   edgecolor="#333355", labelcolor="white")

    fig.text(0.5, 0.985,
             "Pi Anomaly in Off-Diagonal Aubry-André Systems — Full Results",
             ha="center", va="top", color="white",
             fontsize=13, fontweight="bold")
    fig.text(0.5, 0.960,
             f"Bradley Joseph Downey Jr, Independent Researcher (2026)  ·  "
             f"Seed={RANDOM_SEED}  ·  N_phase={N_PHASE}",
             ha="center", va="top", color="#aaaacc", fontsize=9)

    plt.savefig("pi_anomaly_results.png", dpi=150,
                bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close()
    print("Figure saved to: pi_anomaly_results.png")

# ── Entry point ────────────────────────────────────────────────────────────
if __name__ == "__main__":
    results = run_simulation(N_SIZES, N_PHASE, LAMBDAS)
    print_summary(results, N_SIZES, LAMBDAS)
    save_csv(results, N_SIZES, LAMBDAS)
    plot_results(results, N_SIZES, LAMBDAS)

    print("\n" + "=" * 70)
    print("SIMULATION COMPLETE")
    print("Output files:")
    print("  pi_anomaly_results.png  — main figure")
    print("  pi_anomaly_data.csv     — raw data for external analysis")
    print("  pi_anomaly_summary.txt  — key results summary")
    print("=" * 70)
