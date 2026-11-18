# -*- coding: utf-8 -*-
"""
Adim 6 - Bobin adimi varyantlari (kalibre model)
================================================
Orijinal adim 7/9. Burada 8/9 ve 9/9 (tam adim) icin:
  - Bosta 50 Hz: 6 akim noktasi -> senkron endüktans egrisi
  - Kilitli rotor 2, 10, 50 Hz (I = 4.622 A) -> 2D kacak endüktans ve rotor direnci
Sarim sayisi FEA'da 18/katman tutulur; farkli sarim sayilari hibrit modelde olceklenir.
Calistirma: py adim6_adim_varyantlari.py   (~30-35 dakika)   Cikti: sonuclar_adim6.csv
"""
import math, csv, time
import femm
import im_geometry as geo
import adim2_dogrulama as A

geo.R_ROTOR = 62.8 - 0.44
A.LAM_FILL = 0.87
A.SIGMA_AL = 24.5
A.SIGMA_BAR_EFF = 24.5 / A.K_RING

f = open(A.P("sonuclar_adim6.csv"), "w", newline=""); W = csv.writer(f)
W.writerow(["pitch_slots", "test", "f_Hz", "I_rms", "L_hes_H", "X_hes", "Rr_hes"])
t0 = time.time()
nl = A.meas["NL"]
for pitch in (8, 9):
    A.PITCH_SLOTS = pitch
    A._shown = False
    A.build_model()
    print("\nADIM %d/9 - bosta (50 Hz)" % pitch, flush=True)
    A.F_SUPPLY = 50.0
    for i in (0, 2, 4, 6, 8, 10):
        Im = nl["I"][i]
        lam, Pb = A.run_case(Im, 1e-4)
        V, Pin, Z = A.phase_values(Im, lam, Pb, 0.0, 1e-4)
        L = Z.imag / (2 * math.pi * 50)
        print("  I=%.3f A  L=%.4f H" % (Im, L), flush=True)
        W.writerow([pitch, "bosta", 50, Im, round(L, 6), "", ""]); f.flush()
    print("ADIM %d/9 - kilitli rotor" % pitch, flush=True)
    for fr in (2.0, 10.0, 50.0):
        A.F_SUPPLY = fr
        I = 4.622
        lam, Pb = A.run_case(I, 1.0)
        V, Pin, Z = A.phase_values(I, lam, Pb, 0.0, 1.0)
        Rr = Pin / (3 * I * I)
        print("  f=%4.1f Hz  X=%.4f Ohm  Rr=%.4f Ohm" % (fr, Z.imag, Rr), flush=True)
        W.writerow([pitch, "kilitli", fr, I, "", round(Z.imag, 5), round(Rr, 5)]); f.flush()
    femm.closefemm()
f.close()
print("\nBitti (%.0f s). Sonuclar: sonuclar_adim6.csv" % (time.time() - t0))
