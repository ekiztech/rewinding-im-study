# -*- coding: utf-8 -*-
"""
Adim 4 - Son kalibrasyon
========================
Bolum 1: g_eff = 0.44 mm sabit; demir B-H egrisinde H, k_H katsayisiyla olceklenir
         (kesim kenari bozulmasi -> etkin gecirgenlik dususu). k_H = 1.5, 2.0, 3.0
         Bosta: doymus (4.27 A), dirsek (2.85 A), dogrusal (1.42 A) noktalar.
Bolum 2: Rotor direnci dusuk frekansta (2 Hz ~ nominal kayma frekansi) hesaplanir;
         TU Graz Rr(20 C) = 1.012 Ohm degerine gore kafes iletkenligi bulunur.
Calistirma:  py adim4_kalibrasyon2.py      (yaklasik 15-20 dakika)
Cikti: sonuclar_adim4.csv
"""
import json, math, csv, time
import femm
import im_geometry as geo
import adim2_dogrulama as A

BH0 = [tuple(p) for p in A.bh]
f = open(A.P("sonuclar_adim4.csv"), "w", newline=""); W = csv.writer(f)
W.writerow(["bolum", "g_eff_mm", "k_H", "sigma_Al_MSm", "f_Hz", "I_rms", "V_olc", "V_hes", "fark_%", "L_hes_H", "Rr_hes"])
def row(r):
    W.writerow(r); f.flush()

def setup(g, kH, sig_al):
    geo.R_ROTOR = 62.8 - g
    A.bh = [(b, h * kH) for b, h in BH0]
    A.SIGMA_AL = sig_al
    A.SIGMA_BAR_EFF = sig_al / A.K_RING
    A._shown = False
    A.build_model()

t0 = time.time()
nl = A.meas["NL"]
print("BOLUM 1: g_eff = 0.44 mm, demir H olcekleme taramasi (bosta, 50 Hz)")
print("  k_H   I[A]   V_olc   V_hes   Fark[%]", flush=True)
A.F_SUPPLY = 50.0
for kH in (1.5, 2.0, 3.0):
    setup(0.44, kH, 26.0)
    for i in (0, 4, 10):
        Im, Vm = nl["I"][i], nl["V"][i]
        lam, Pb = A.run_case(Im, 1e-4)
        V, Pin, Z = A.phase_values(Im, lam, Pb, nl["Rs"], 1e-4)
        e = 100 * (V - Vm) / Vm
        print("  %.1f   %.3f  %6.1f  %6.1f  %+6.1f" % (kH, Im, Vm, V, e), flush=True)
        row(["bosta", 0.44, kH, 26.0, 50, Im, Vm, round(V, 2), round(e, 2), round(Z.imag / (2 * math.pi * 50), 5), ""])
    femm.closefemm()

print("\nBOLUM 2: dusuk frekans rotor direnci (2 Hz), sigma_Al = 26 MS/m")
setup(0.44, 1.0, 26.0)
A.F_SUPPLY = 2.0
I = 4.622
lam, Pb = A.run_case(I, 1.0)
Pb = abs(Pb.real if isinstance(Pb, complex) else Pb)
Rr = Pb / (3 * I * I)
sig_new = 26.0 * Rr / 1.012
print("  Rr_model(2 Hz) = %.3f Ohm   -> hedef 1.012 Ohm icin sigma_Al = %.1f MS/m" % (Rr, sig_new), flush=True)
row(["rotor_2Hz", 0.44, 1.0, 26.0, 2, I, "", "", "", "", round(Rr, 4)])
femm.closefemm()
f.close()
print("\nBitti (%.0f s). Sonuclar: sonuclar_adim4.csv" % (time.time() - t0))
