# -*- coding: utf-8 -*-
"""
Adim 7 - Kafes iletkenligi duzeltmesi (sigma_Al = 29.3 MS/m) ile kilitli rotor tekrar
=====================================================================================
Neden: Adim 4'te rotor direnci P/3I^2 ile (gorunur deger) kalibre edilmisti; miknatislanma
kolu ayristirilinca dogru iletkenlik ~29.3 MS/m cikti. Bosta sonuclari degismez.
Adim 7/9: 2, 10, 20, 30, 40, 50 Hz (olcum akimlariyla) ; Adim 8/9 ve 9/9: 2, 10, 50 Hz
Calistirma: py adim7_rotor_yeniden.py   (~20-25 dakika)   Cikti: sonuclar_adim7.csv
"""
import json, math, csv, time
import femm
import im_geometry as geo
import adim2_dogrulama as A

SIG = 29.3
geo.R_ROTOR = 62.8 - 0.44
A.LAM_FILL = 0.87
A.SIGMA_AL = SIG
A.SIGMA_BAR_EFF = SIG / A.K_RING

f = open(A.P("sonuclar_adim7.csv"), "w", newline=""); W = csv.writer(f)
W.writerow(["pitch_slots", "f_Hz", "I_rms", "R_app_hes", "X_hes"])
t0 = time.time()
lrf = json.load(open(A.P("lr_freq.json")))
for pitch in (7, 8, 9):
    A.PITCH_SLOTS = pitch
    A._shown = False
    A.build_model()
    pts = [(2.0, 4.622)] + ([(d["f"], d["I"]) for d in lrf] if pitch == 7 else [(10.0, 4.622), (50.0, 4.622)])
    print("\nADIM %d/9 - kilitli rotor (sigma_Al = %.1f MS/m)" % (pitch, SIG), flush=True)
    for fr, I in pts:
        A.F_SUPPLY = fr
        lam, Pb = A.run_case(I, 1.0)
        V, Pin, Z = A.phase_values(I, lam, Pb, 0.0, 1.0)
        R = Pin / (3 * I * I)
        print("  f=%5.2f Hz  I=%.3f  R=%.4f  X=%.4f Ohm" % (fr, I, R, Z.imag), flush=True)
        W.writerow([pitch, fr, I, round(R, 5), round(Z.imag, 5)]); f.flush()
    femm.closefemm()
f.close()
print("\nBitti (%.0f s). Sonuclar: sonuclar_adim7.csv" % (time.time() - t0))
