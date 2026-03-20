# Pi Anomaly Simulation Package
## Off-Diagonal Aubry-André Systems: Irrational Parameter Dependence

**Author:** Bradley Joseph Downey Jr, Independent Researcher (2026)  
**Preprint:** https://doi.org/10.5281/zenodo.19101668
**Patents:** USPTO 64/008,843 and 64/009,939 (filed March 18, 2026)

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

Off-diagonal Aubry-André Hamiltonian (open ## What This Is

Numerical simulation code comparing irrational modulation parameters
in the off-diagonal Aubry-André model. Establishes the Spectral Gap
Resonance (SGR) mechanism governing finite-size localization across
24 distinct irrational parameters and 2,147 phase realisations.

**Key finding:** At the localization transition λ/t = 2.0, the choice
of modulation parameter α substantially affects spectral properties
at all experimentally relevant finite system sizes. The conventional
choice α = φ is suboptimal. Engineered irrationals of the form
[1;K,1,K,...] outperform both π and φ by factors of 1.4× to 3.3×
in localization length, confirmed at 99% confidence.

**Design principle:** For a quasi-periodic device of N elements:
- N < 50: use K = floor(N/2), α = [1;K,1,K,...]
- N ≥ 50: use K = N−2, α = [1;K,1,K,...]

Minimum effective device size: N = 6.

---

## Results Summary

| N | Best α | vs φ | 99% CI |
|---|--------|------|--------|
| 11 | [1;4,1,4,...] ≈ 0.8284 | 1.49× | >1.46× |
| 20 | [1;10,1,10,...] ≈ 0.9161 | 1.65× | >1.61× |
| 50 | [1;3,1,3,...] ≈ 0.7913 | 2.43× | >2.38× |
| 200–1500 | [1;100,1,100,...] ≈ 0.9902 | 1.84–3.18× | confirmed |

---

## Files

- `pi_anomaly_simulation.py` — full simulation, N = 200 to 3000
- `pi_anomaly_quickstart.py` — 30-second verification run
- `hsl_preprint_v10.pdf` — full preprint (v2.0)
- `README.md` — this file
- `PATENT_NOTICE.md` — patent and license information

---

## Citation

Downey, B.J. (2026). Irrational Modulation Parameter Dependence
in Aubry-André Systems: Spectral Gap Resonance, Super-GOE Level
Statistics, and a Scale-Dependent Device Design Principle from
Continued Fraction Structure. Zenodo v2.0.
https://doi.org/10.5281/zenodo.19101668

---

## Patent Notice

The methods in this repository are protected under:

- USPTO Provisional Application No. **64/008,843** (filed March 18, 2026)
  Covers: π modulation parameter advantage, super-GOE level statistics

- USPTO Provisional Application No. **64/009,939** (filed March 18, 2026)
  Covers: Spectral Gap Resonance mechanism, engineered irrationals
  [1;K,1,K,...], two-regime design principle, device design method

Academic and non-commercial use permitted without restriction.
Commercial use requires a license. Contact: freeusevids@gmail.com
