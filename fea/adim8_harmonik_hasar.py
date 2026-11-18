# -*- coding: utf-8 -*-
"""
Adim 8 - Hava araligi harmonikleri (3 adim) + hasarli cekirdek gecirgenligi
==========================================================================
Bolum A: Adim 7/9, 8/9, 9/9 icin nominal miknatislanma akiminda (3.2 A) manyetostatik cozum;
         hava araligi ortasinda (r = 62.58 mm) 1440 noktada radyal aki yogunlugu B_r kaydedilir.
         (Harmonik analizi Claude yapar.)
Bolum B: Cekirdek hasari - etkin sac dolulugu 0.84 ve 0.80 (kalibre: 0.87), adim 7/9,
         6 bosta noktasi -> miknatislanma egrisi.
Calistirma: py adim8_harmonik_hasar.py   (~25-30 dakika)
Cikti: sonuclar_adim8_bgap.csv, sonuclar_adim8_hasar.csv
"""
import math, csv, time
import femm
import im_geometry as geo
import adim2_dogrulama as A

geo.R_ROTOR = 62.8 - 0.44
A.SIGMA_AL = 29.3; A.SIGMA_BAR_EFF = 29.3 / A.K_RING
t0 = time.time()
R_MID = 62.8 - 0.22
NPTS = 1440

# ---------------- Bolum A ----------------
fb = open(A.P("sonuclar_adim8_bgap.csv"), "w", newline=""); WB = csv.writer(fb)
WB.writerow(["pitch_slots", "angle_deg", "Br_T", "Bt_T"])
A.LAM_FILL = 0.87
A.F_SUPPLY = 50.0
for pitch in (7, 8, 9):
    A.PITCH_SLOTS = pitch
    A._shown = False
    A.build_model()
    I = 3.201; Ipk = math.sqrt(2) * I
    femm.mi_probdef(0, "millimeters", "planar", 1e-8, A.L_STACK, 30, 0)
    for ph, z in (("A", Ipk), ("B", -Ipk / 2), ("C", -Ipk / 2)):
        getattr(femm, "callfemm_noeval", femm.callfemm)('mi_modifycircprop("%s",1,%.10g)' % (ph, z))
    femm.mi_analyze(1); femm.mi_loadsolution()
    print("ADIM %d/9: hava araligi ornekleniyor ..." % pitch, flush=True)
    for k in range(NPTS):
        th = 2 * math.pi * k / NPTS
        x, y = R_MID * math.cos(th), R_MID * math.sin(th)
        bx, by = femm.mo_getb(x, y)
        br = bx * math.cos(th) + by * math.sin(th); bt = -bx * math.sin(th) + by * math.cos(th)
        WB.writerow([pitch, round(math.degrees(th), 4), round(br, 6), round(bt, 6)])
    fb.flush()
    femm.mo_close(); femm.closefemm()
fb.close()

# ---------------- Bolum B ----------------
fh = open(A.P("sonuclar_adim8_hasar.csv"), "w", newline=""); WH = csv.writer(fh)
WH.writerow(["lam_fill", "I_rms", "L_hes_H"])
A.PITCH_SLOTS = 7
nl = A.meas["NL"]
for fill in (0.84, 0.80):
    A.LAM_FILL = fill
    A._shown = False
    A.build_model()
    print("\nCekirdek hasari: etkin doluluk %.2f" % fill, flush=True)
    for i in (0, 2, 4, 6, 8, 10):
        Im = nl["I"][i]
        lam, Pb = A.run_case(Im, 1e-4)
        V, Pin, Z = A.phase_values(Im, lam, Pb, 0.0, 1e-4)
        L = Z.imag / (2 * math.pi * 50)
        print("  I=%.3f A  L=%.4f H" % (Im, L), flush=True)
        WH.writerow([fill, Im, round(L, 6)]); fh.flush()
    femm.closefemm()
fh.close()
print("\nBitti (%.0f s). Sonuclar: sonuclar_adim8_bgap.csv, sonuclar_adim8_hasar.csv" % (time.time() - t0))
