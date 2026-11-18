# -*- coding: utf-8 -*-
"""
Adim 5 - Kalibre edilmis modelin son dogrulamasi
================================================
Kalibre parametreler:
  g_eff      = 0.44 mm   (bosta, dogrusal bolge)
  lam_fill   = 0.87      (bosta, doyma bolgesi; etkin deger)
  sigma_Al   = 24.5 MS/m (rotor direnci, 2 Hz; TU Graz Rr(20C)=1.012 Ohm)
Testler:
  - Bosta 50 Hz: 11 olcum noktasinin tamami
  - Kilitli rotor 10-50 Hz: Rr(f) ve 3D uc endüktansi L_3D (olculen - 2D)
Calistirma: py adim5_son_dogrulama.py   (~25-30 dakika)   Cikti: sonuclar_adim5.csv
"""
import json, math, csv, time
import femm
import im_geometry as geo
import adim2_dogrulama as A

G_EFF, FILL, SIG = 0.44, 0.87, 24.5
geo.R_ROTOR = 62.8 - G_EFF
A.LAM_FILL = FILL
A.SIGMA_AL = SIG
A.SIGMA_BAR_EFF = SIG / A.K_RING

f = open(A.P("sonuclar_adim5.csv"), "w", newline=""); W = csv.writer(f)
W.writerow(["test", "f_Hz", "I_rms", "V_olc", "V_hes", "fark_%", "P_olc", "P_hes", "L_hes_H", "Rr_olc", "Rr_hes", "X_olc", "X_hes", "L3D_mH"])
t0 = time.time()
A._shown = False
A.build_model()

nl = A.meas["NL"]; A.F_SUPPLY = 50.0
print("BOSTA (50 Hz) - kalibre model")
print("  I[A]   V_olc   V_hes   Fark[%]", flush=True)
for Im, Vm in zip(nl["I"], nl["V"]):
    lam, Pb = A.run_case(Im, 1e-4)
    V, Pin, Z = A.phase_values(Im, lam, Pb, nl["Rs"], 1e-4)
    e = 100 * (V - Vm) / Vm
    print("  %.3f  %6.1f  %6.1f  %+6.1f" % (Im, Vm, V, e), flush=True)
    W.writerow(["bosta", 50, Im, Vm, round(V, 2), round(e, 2), "", "", round(Z.imag / (2 * math.pi * 50), 5), "", "", "", "", ""]); f.flush()

print("\nKILITLI ROTOR (10-50 Hz) - kalibre model")
print("  f[Hz]  I[A]   Rr_olc  Rr_hes   X_olc  X_hes   L_3D[mH]", flush=True)
L3D = []
for d in json.load(open(A.P("lr_freq.json"))):
    fr, I = d["f"], d["I"]
    A.F_SUPPLY = fr
    lam, Pb = A.run_case(I, 1.0)
    V, Pin, Z = A.phase_values(I, lam, Pb, d["Rs"], 1.0)
    Zm = d["V"] / I; Rm = d["P"] / (3 * I * I); Xm = math.sqrt(max(Zm * Zm - Rm * Rm, 0))
    Rh = Pin / (3 * I * I); Xh = Z.imag
    dL = (Xm - Xh) / (2 * math.pi * fr) * 1e3; L3D.append(dL)
    print("  %5.1f  %.3f  %.3f   %.3f   %.3f  %.3f   %.2f" % (fr, I, Rm - d["Rs"], Rh - d["Rs"], Xm, Xh, dL), flush=True)
    W.writerow(["kilitli", fr, I, d["V"], round(V, 2), round(100 * (V - d["V"]) / d["V"], 2), d["P"], round(Pin, 1), "",
                round(Rm - d["Rs"], 4), round(Rh - d["Rs"], 4), round(Xm, 4), round(Xh, 4), round(dL, 3)]); f.flush()
femm.closefemm()
f.close()
print("\nOrtalama L_3D = %.2f mH (min %.2f, max %.2f)" % (sum(L3D) / len(L3D), min(L3D), max(L3D)))
print("Bitti (%.0f s). Sonuclar: sonuclar_adim5.csv" % (time.time() - t0))
