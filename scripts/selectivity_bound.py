#!/usr/bin/env python3
"""WAVE 24: honest pH-selectivity bound -- Henderson-Hasselbalch protonation x contact free energy, with
the electronic-continuum-correction (ECC) charge-scaling applied to each chemistry's charged partners.

Selectivity (tumor-ON fold change) = K(6.5)/K(7.4) = exp[(f(6.5)-f(7.4)) * sum(dG_contact) / RT], where
f(pH)=1/(1+10^(pH-pKa)). A salt bridge scales BOTH charges under ECC (~0.8*0.8=0.64); a cation-pi scales only
the His (~0.8, and fixed-charge FFs UNDER-count cation-pi, so this is conservative). Shows why the robust
chemistry cannot beat the salt bridge's RAW selectivity, only (at dual) tie its ECC-honest selectivity.
"""
import numpy as np
RT = 0.616  # kcal/mol, 310 K
f = lambda pH, pKa: 1.0 / (1.0 + 10 ** (pH - pKa))
sel = lambda dG, pKa: float(np.exp((f(6.5, pKa) - f(7.4, pKa)) * dG / RT))
DG_SALT, ECC_SALT = 3.5, 0.64
DG_CPI, ECC_CPI = 1.4, 0.80
cols = ["1xSalt raw", "1xSalt ECC", "1xCatpi raw", "1xCatpi ECC", "2xCatpi raw", "2xCatpi ECC", "3xCatpi ECC"]
print(f"{'pKa':>5} | " + " ".join(f"{c:>11}" for c in cols))
print("-" * 98)
for pKa in (5.5, 6.0, 6.3, 6.5, 6.95, 7.4):
    v = [sel(DG_SALT, pKa), sel(DG_SALT * ECC_SALT, pKa), sel(DG_CPI, pKa), sel(DG_CPI * ECC_CPI, pKa),
         sel(2 * DG_CPI, pKa), sel(2 * DG_CPI * ECC_CPI, pKa), sel(3 * DG_CPI * ECC_CPI, pKa)]
    print(f"{pKa:5.2f} | " + " ".join(f"{x:10.1f}x" for x in v))
print("\nBar to beat (PASS20 inv01/inv02 single salt bridge, ECC-honest) ~5.7x at pKa 6.95.")
print("1x cation-pi: loses. 2x cation-pi: ties on honest selectivity, far more ROBUST. 3x cation-pi: ~13x (far better) but needs 3 on-face conserved aromatics.")
