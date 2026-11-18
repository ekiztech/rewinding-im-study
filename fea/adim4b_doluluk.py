# -*- coding: utf-8 -*-
"""
Adim 4b - Etkin sac doluluk faktoru taramasi (doyma bolgesi kalibrasyonu)
g_eff = 0.44 mm, k_H = 1, LAM_FILL = 0.93 / 0.90 / 0.87
Bosta: 4.27 A (doymus), 2.85 A (dirsek), 1.42 A (dogrusal)
Calistirma: py adim4b_doluluk.py   (~15 dakika)   Cikti: sonuclar_adim4b.csv
"""
import math, csv, time
import femm
import im_geometry as geo
import adim2_dogrulama as A

f = open(A.P("sonuclar_adim4b.csv"), "w", newline=""); W = csv.writer(f)
W.writerow(["g_eff_mm", "lam_fill", "I_rms", "V_olc", "V_hes", "fark_%", "L_hes_H"])
t0 = time.time()
nl = A.meas["NL"]
A.F_SUPPLY = 50.0
print("g_eff = 0.44 mm, etkin sac doluluk faktoru taramasi (bosta, 50 Hz)")
print("  dol.   I[A]   V_olc   V_hes   Fark[%]", flush=True)
for fill in (0.93, 0.90, 0.87):
    geo.R_ROTOR = 62.8 - 0.44
    A.LAM_FILL = fill
    A._shown = False
    A.build_model()
    for i in (0, 4, 10):
        Im, Vm = nl["I"][i], nl["V"][i]
        lam, Pb = A.run_case(Im, 1e-4)
        V, Pin, Z = A.phase_values(Im, lam, Pb, nl["Rs"], 1e-4)
        e = 100 * (V - Vm) / Vm
        print("  %.2f  %.3f  %6.1f  %6.1f  %+6.1f" % (fill, Im, Vm, V, e), flush=True)
        W.writerow([0.44, fill, Im, Vm, round(V, 2), round(e, 2), round(Z.imag / (2 * math.pi * 50), 5)]); f.flush()
    femm.closefemm()
f.close()
print("\nBitti (%.0f s). Sonuclar: sonuclar_adim4b.csv" % (time.time() - t0))
