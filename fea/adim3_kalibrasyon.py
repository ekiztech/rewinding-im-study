# -*- coding: utf-8 -*-
"""
Adim 3 - Model kalibrasyonu ve bagimsiz dogrulama
=================================================
Bolum 1: Etkin hava araligi g_eff = 0.43 ve 0.46 mm icin bosta egrisi (5 nokta).
         (0.40 mm sonuclari Adim 2'den var.)
Bolum 2: Rotor kafes iletkenligi Rr(50 Hz) olcumune gore ayarlanir; sonra kilitli rotor
         testi 10-20-30-40-50 Hz'de cozulur. Rr(f) egrisi (deri etkisi) ve eksik kacak
         endüktans (sargi basi) frekansa gore karsilastirilir -> bagimsiz dogrulama.
Calistirma:  py adim3_kalibrasyon.py      (yaklasik 20-30 dakika)
Cikti: sonuclar_adim3.csv
"""
import os, json, math, csv, time
import femm
import im_geometry as geo
import adim2_dogrulama as A

OUT = A.P("sonuclar_adim3.csv")
f = open(OUT, "w", newline=""); W = csv.writer(f)
W.writerow(["bolum", "g_eff_mm", "sigma_Al_MSm", "f_Hz", "I_rms", "V_olc", "V_hes", "P_olc", "P_hes",
            "L_hes_H", "Rr_olc", "Rr_hes", "X_olc", "X_hes", "dL_mH"])
def row(r):
    W.writerow(r); f.flush()

def set_sigma(sig_al):
    A.SIGMA_AL = sig_al
    A.SIGMA_BAR_EFF = sig_al / A.K_RING      # K_RING oranı iletkenlikten bagimsiz

t0 = time.time()
nl = A.meas["NL"]
idx = [0, 2, 4, 6, 10]                       # 4.27, 3.56, 2.85, 2.14, 1.42 A

print("BOLUM 1: etkin hava araligi taramasi (bosta, 50 Hz)")
print("  g[mm]   I[A]   V_olc   V_hes   Fark[%]", flush=True)
for g in (0.43, 0.46):
    geo.R_ROTOR = 62.8 - g
    set_sigma(26.0)
    A._shown = False
    A.build_model()
    for i in idx:
        Im, Vm = nl["I"][i], nl["V"][i]
        lam, Pb = A.run_case(Im, 1e-4)
        V, Pin, Z = A.phase_values(Im, lam, Pb, nl["Rs"], 1e-4)
        L = Z.imag / (2 * math.pi * 50)
        print("  %.2f   %.3f  %6.1f  %6.1f  %+6.1f" % (g, Im, Vm, V, 100 * (V - Vm) / Vm), flush=True)
        row(["bosta", g, 26.0, 50, Im, Vm, round(V, 2), "", "", round(L, 5), "", "", "", "", ""])
    femm.closefemm()

print("\nBOLUM 2: kilitli rotor, 10-50 Hz (g_eff = 0.43 mm, sigma_Al = 23.2 MS/m)")
print("  f[Hz]  I[A]   Rr_olc  Rr_hes   X_olc  X_hes   eksik L [mH]", flush=True)
geo.R_ROTOR = 62.8 - 0.43
set_sigma(23.2)                               # 26 * 1.580/1.769 : Rr(50 Hz) olcumune gore
A._shown = False
A.build_model()
lrf = json.load(open(A.P("lr_freq.json")))
for d in lrf:
    fr = d["f"]
    A.F_SUPPLY = fr
    lam, Pb = A.run_case(d["I"], 1.0)
    V, Pin, Z = A.phase_values(d["I"], lam, Pb, d["Rs"], 1.0)
    I = d["I"]
    Zm = d["V"] / I; Rm = d["P"] / (3 * I * I); Xm = math.sqrt(max(Zm * Zm - Rm * Rm, 0))
    Rh = Pin / (3 * I * I); Xh = Z.imag
    dL = (Xm - Xh) / (2 * math.pi * fr) * 1e3
    print("  %5.1f  %.3f  %.3f   %.3f   %.3f  %.3f   %.2f" % (fr, I, Rm - d["Rs"], Rh - d["Rs"], Xm, Xh, dL), flush=True)
    row(["kilitli", 0.43, 23.2, fr, I, d["V"], round(V, 2), d["P"], round(Pin, 1), "",
         round(Rm - d["Rs"], 4), round(Rh - d["Rs"], 4), round(Xm, 4), round(Xh, 4), round(dL, 3)])
femm.closefemm()
f.close()
print("\nBitti (%.0f s). Sonuclar: sonuclar_adim3.csv" % (time.time() - t0))
