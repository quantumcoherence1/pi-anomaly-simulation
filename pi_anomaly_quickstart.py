"""
================================================================================
PI ANOMALY — QUICK VERIFICATION SCRIPT
================================================================================

Runs N=200 only with 3 phase realisations.
Takes approximately 30 seconds.
Should show π leading on localization length and fractal dimension.

For full simulation across N=200 to N=3000, run:
    python pi_anomaly_simulation.py
================================================================================
"""

import numpy as np
from scipy.linalg import eigh
import time

PI    = np.pi
PHI   = (1 + np.sqrt(5)) / 2
SQRT2 = np.sqrt(2)
E     = np.e
GEN   = 1.2020569 / 10 * 3

ALPHAS = {"π": PI, "φ": PHI, "√2": SQRT2, "e": E, "generic": GEN}
COLORS_TEXT = {"π": "***", "φ": "   ", "√2": "   ", "e": "   ", "generic": "   "}

rng = np.random.default_rng(42)

def build_H(N, lam, alpha, phi):
    n = np.arange(N - 1)
    h = 1.0 + lam * np.cos(2 * PI * alpha * n + phi)
    H = np.zeros((N, N))
    H[n, n+1] = h; H[n+1, n] = h
    return H

def loc_len(evecs, N):
    mid = N // 2; half = max(1, int(0.10 * N / 2))
    sel = evecs[:, mid-half:mid+half]
    ipr = np.sum(sel**4, axis=0)
    v = ipr[ipr > 0]
    return float(np.mean(1.0/v)) if len(v) else 0.0

def r_stat(evals):
    ev = np.sort(evals); sp = np.diff(ev); sp = sp[sp > 1e-12]
    if len(sp) < 2: return np.nan
    return float(np.mean(np.minimum(sp[:-1],sp[1:])/np.maximum(sp[:-1],sp[1:])))

def frac_dim(evals, boxes=(10,20,50,100,200)):
    ev = np.sort(evals); span = ev[-1]-ev[0]
    if span < 1e-10: return 0.0
    counts = [len(set(int((e-ev[0])/span*nb) for e in ev)) for nb in boxes]
    ls = np.log(np.array(boxes,float)); lc = np.log(np.array(counts,float))
    return float(np.polyfit(ls,lc,1)[0]) if len(set(lc))>1 else 1.0

print("=" * 60)
print("PI ANOMALY — QUICK VERIFICATION (N=200, λ/t=2.0)")
print("=" * 60)
print(f"\n{'α':>10}  {'ξ (loc len)':>12}  {'r-stat':>8}  {'D (frac)':>10}")
print("-" * 50)

results = {}
for name, alpha in ALPHAS.items():
    lls, rrs, fds = [], [], []
    for _ in range(3):
        phi = rng.uniform(0, 2*PI)
        H = build_H(200, 2.0, alpha, phi)
        ev, evec = eigh(H)
        lls.append(loc_len(evec, 200))
        rrs.append(r_stat(ev))
        fds.append(frac_dim(ev))
    results[name] = {"xi": np.mean(lls), "r": np.nanmean(rrs), "D": np.mean(fds)}
    flag = " <-- π leads" if name == "π" else ""
    print(f"{name:>10}  {results[name]['xi']:>12.2f}  "
          f"{results[name]['r']:>8.4f}  {results[name]['D']:>10.4f}{flag}")

print("\n" + "=" * 60)
pi_phi_ratio = results["π"]["xi"] / results["φ"]["xi"]
print(f"π/φ localization length ratio at N=200: {pi_phi_ratio:.3f}")
print(f"Expected from paper: ~1.16")
print(f"Match: {'YES ✓' if 1.0 < pi_phi_ratio < 1.5 else 'CHECK - may need more phase averaging'}")
print("\nTo run full simulation across N=200 to N=3000:")
print("    python pi_anomaly_simulation.py")
print("=" * 60)
