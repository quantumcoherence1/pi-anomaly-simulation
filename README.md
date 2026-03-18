# Pi Anomaly Simulation Package
## Off-Diagonal Aubry-André Systems: Irrational Parameter Dependence

**Author:** Bradley Joseph Downey Jr, Independent Researcher (2026)  
**Preprint:** [arXiv link — to be added]  
**Patent:** USPTO Provisional Application (filed 2026)

---

## What This Is

Numerical simulation code comparing five irrational modulation parameters
(π, φ, √2, e, generic) in the off-diagonal Aubry-André model.

**Key finding:** At the localization transition λ/t = 2.0, α=π produces
super-linearly growing localization length advantages over α=φ:

| N    | π      | φ     | π/φ ratio |
|------|--------|-------|-----------|
| 200  | ~15    | ~13   | 1.16      |
| 500  | ~39    | ~23   | 1.67      |
| 1000 | ~82    | ~32   | 2.59      |
| 2000 | ~149   | ~53   | 2.83      |
| 3000 | ~210   | ~65   | 3.23      |

No plateau observed. Spectral fractal dimension gap widens simultaneously.
Anomalous r-statistic r=0.827 at N=3000 exceeds GOE ceiling of 0.530.

---

## Files

```
pi_anomaly_simulation.py   — Full simulation (N=200 to N=3000, ~50 min)
pi_anomaly_quickstart.py   — Quick verification run (N=200 only, ~30 sec)
README.md                  — This file
```

---

## Requirements

```
Python 3.8+
numpy
scipy
matplotlib
```

Install:
```bash
pip install numpy scipy matplotlib
```

---

## Quick Start (30 seconds)

To verify the core finding immediately:

```bash
python pi_anomaly_quickstart.py
```

This runs N=200 only with 3 phase realisations. You should see π showing
higher localization length and fractal dimension than all other irrationals.

---

## Full Simulation (~50 minutes)

```bash
python pi_anomaly_simulation.py
```

Produces three output files:
- `pi_anomaly_results.png`   — main figure with all metrics
- `pi_anomaly_data.csv`      — raw data for external analysis
- `pi_anomaly_summary.txt`   — printed summary of key results

---

## Reproducing Specific System Sizes

Edit the N_SIZES list in pi_anomaly_simulation.py:

```python
# Run only N=500 and N=1000
N_SIZES = [500, 1000]
```

---

## The Model

Off-diagonal Aubry-André Hamiltonian (open boundary conditions):

```
H = t * Σ_{n=0}^{N-2} (1 + λ·cos(2π·α·n + φ)) * |n+1><n| + h.c.
```

- t = 1.0 (energy scale)
- λ/t = 2.0 (transition point, main result)
- α ∈ {π, φ, √2, e, generic}
- φ averaged over N_PHASE uniform realisations in [0, 2π)

Off-diagonal (hopping) modulation is used rather than on-site modulation
because it preserves chiral symmetry and is physically realistic for
laser-modulated photonic and quantum systems.

---

## Three Metrics

**Localization length ξ:**  
ξ = 1/IPR, where IPR = Σ_j |ψ_j|⁴  
Averaged over central 10% of eigenstates (near E=0)  
Higher = better coherence for quantum devices

**Level-spacing r-statistic:**  
r_n = min(s_n, s_{n+1}) / max(s_n, s_{n+1})  
Poisson (localized) = 0.386  
GOE (delocalized) = 0.530  
π shows r=0.827 at N=3000 — anomalously above GOE ceiling

**Spectral fractal dimension D:**  
Box-counting dimension of the eigenvalue set  
D < 1 = Cantor-set-like spectrum  
π shows D≈0.52-0.55 consistently, others show D=0.27-0.40

---

## Falsifiable Predictions

This code can directly test:

- **P1:** π/φ ratio > 4.0 at N=5000  
  → Edit N_SIZES = [5000] and run

- **P2:** D(π) > 0.55, D(φ) < 0.25 at N=5000  
  → Same run, check frac_dim output

- **P3:** r(π) > 0.80 confirmed with N_phase ≥ 20 at N=3000  
  → Edit N_PHASE = 20, N_SIZES = [3000] and run

---

## Contact

Bradley Joseph Downey Jr  
Independent Researcher  
[email]
