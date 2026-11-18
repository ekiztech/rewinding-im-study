# -*- coding: utf-8 -*-
"""
Adim 2 - CREATOR asenkron motorunun FEMM modeli ve dogrulamasi
================================================================
Ne yapar:
  1) Geometriyi (im_geometry.py, DXF'ten cikarilan temiz parametrik model) FEMM'de cizer, malzeme / sargi / sinir kosullarini atar.
  2) 50 Hz bosta (s ~ 0) ve kilitli rotor (s = 1) testlerini, olculen akimlarda simule eder.
  3) Hesaplanan gerilim ve gucu TU Graz olcumleriyle karsilastirir.
Cikti:
  - IM_model.fem               (FEMM model dosyasi, istersen GUI'de acabilirsin)
  - sonuclar_adim2.csv         (tum karsilastirma tablosu)
Calistirma:  py adim2_dogrulama.py
Veri kaynagi: CREATOR Case IM, TU Graz, DOI 10.3217/kh3v7-dhn98 (CC BY-NC 4.0)
"""
import os, json, math, cmath, csv, time
import femm
import im_geometry as geo

HERE = os.path.dirname(os.path.abspath(__file__))
P = lambda f: os.path.join(HERE, f).replace("\\", "/")

# ----------------------------------------------------------------------------
# 1. MOTOR VERILERI (CREATOR veri setinden)
# ----------------------------------------------------------------------------
F_SUPPLY   = 50.0          # Hz
POLE_PAIRS = 2
N_SLOT     = 36
N_BAR      = 28
L_STACK    = 100.0         # mm (2D derinlik)
TURNS_LAYER = 18           # her oluk katmanindaki sarim sayisi
PITCH_SLOTS = 7            # adim 7/9
SIGMA_AL   = 26.0          # MS/m, Al 6082 (veri seti 24-32 MS/m araligi veriyor)
LAM_FILL   = 0.97          # sac doluluk faktoru (varsayim)

# Uc halkasi duzeltmesi: halka direncini cubuk iletkenligine yansitiyoruz
A_BAR  = 68.3e-6           # m2, DXF'ten hesaplanan cubuk kesiti
D_RING = (0.081 + 0.1242) / 2
A_RING = 0.014 * (0.1242 - 0.081) / 2
R_bar  = (L_STACK / 1000) / (SIGMA_AL * 1e6 * A_BAR)
R_seg  = math.pi * D_RING / (N_BAR * SIGMA_AL * 1e6 * A_RING)
R_ring_eq = R_seg / (2 * math.sin(math.pi * POLE_PAIRS / N_BAR) ** 2)
K_RING = 1 + R_ring_eq / R_bar
SIGMA_BAR_EFF = SIGMA_AL / K_RING

meas = json.load(open(P("meas50.json")))
bh = json.load(open(P("bh_m800_50a.json")))

# ----------------------------------------------------------------------------
# 2. SARGI YERLESIMI (cift katman, q=3, adim 7/9)
# ----------------------------------------------------------------------------
BELTS = [("A", 1), ("C", -1), ("B", 1), ("A", -1), ("C", 1), ("B", -1)]
def upper_layer(k):            # ust katman (hava araligina yakin)
    return BELTS[(k // 3) % 6]
def lower_layer(k):            # alt katman = 7 oluk geriden gelen bobinin donus kenari
    ph, sg = upper_layer((k - PITCH_SLOTS) % N_SLOT)
    return ph, -sg

# ----------------------------------------------------------------------------
# 3. MODELI KUR
# ----------------------------------------------------------------------------
def build_model():
    nodes, segs, arcs, labels = geo.build()
    femm.openfemm(1)                      # 1: pencere gizli calisir
    femm.newdocument(0)
    femm.mi_probdef(0.001, "millimeters", "planar", 1e-8, L_STACK, 30, 0)
    getattr(femm, "callfemm_noeval", femm.callfemm)("mi_smartmesh(0)")

    print("  %d dugum, %d dogru, %d yay ciziliyor ..." % (len(nodes), len(segs), len(arcs)), flush=True)
    for x, y in nodes:
        femm.mi_addnode(x, y)
    for (x1, y1), (x2, y2) in segs:
        femm.mi_addsegment(x1, y1, x2, y2)
    for (x1, y1), (x2, y2), sweep in arcs:
        chord = math.hypot(x2 - x1, y2 - y1)
        radius = chord / (2 * math.sin(math.radians(sweep) / 2))
        # her yay parcasi ~0.2 mm (hava araligi yaylari) ile ~0.5 mm arasinda; kucuk yuvarlatmalarda en fazla 15 derece
        target = 0.2 if 60 < radius < 64 else 0.5
        maxseg = min(15.0, max(0.25, math.degrees(target / radius)))
        femm.mi_addarc(x1, y1, x2, y2, sweep, maxseg)

    # Malzemeler
    femm.mi_addmaterial("Air", 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
    femm.mi_addmaterial("Coil", 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)  # direnc analitik
    femm.mi_addmaterial("M800-50A", 2300, 2300, 0, 0, 0, 0.5, 0, LAM_FILL, 0, 0, 0, 0, 0)
    for b, h in bh[1:]:
        femm.mi_addbhpoint("M800-50A", b, h)
    femm.mi_addmaterial("Al_bar", 1, 1, 0, 0, SIGMA_BAR_EFF, 0, 0, 1, 0, 0, 0, 0, 0)
    femm.mi_addmaterial("Shaft", 100, 100, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
    for ph in "ABC":
        femm.mi_addcircprop(ph, 0, 1)    # seri devre
    femm.mi_addboundprop("A0", 0, 0, 0, 0, 0, 0, 0, 0, 0)

    mat = {"stator_iron": ("M800-50A", 2.5), "rotor_iron": ("M800-50A", 2.0),
           "airgap": ("Air", 0.2), "liner": ("Air", 0.4), "shaft": ("Shaft", 5.0),
           "bar": ("Al_bar", 0.7)}
    for role, x, y, k in labels:
        femm.mi_addblocklabel(x, y)
        femm.mi_selectlabel(x, y)
        if role in ("coil_up", "coil_lo"):
            ph, sg = upper_layer(k) if role == "coil_up" else lower_layer(k)
            femm.mi_setblockprop("Coil", 0, 1.5, ph, 0, 1, sg * TURNS_LAYER)
        else:
            m, ms = mat[role]
            grp = 2 if role == "bar" else 0
            femm.mi_setblockprop(m, 0, ms, "<None>", 0, grp, 0)
        femm.mi_clearselected()

    # Dis sinir: A = 0
    femm.mi_seteditmode("arcsegments")
    for ang in (90, 270):
        femm.mi_selectarcsegment(100 * math.cos(math.radians(ang)), 100 * math.sin(math.radians(ang)))
    femm.mi_setarcsegmentprop(2, "A0", 0, 1)
    femm.mi_clearselected()
    femm.mi_saveas(P("IM_model.fem"))

# ----------------------------------------------------------------------------
# 4. TEK CALISMA: verilen RMS akim ve kayma icin
# ----------------------------------------------------------------------------
_shown = False
def run_case(I_rms, slip, bitmap=None):
    """Bosta (slip ~ 0): rotor akimi yok -> manyetostatik dogrusal olmayan cozum (hizli, saglam yakinsama).
       Faz A akimi tepe degerinde iken (Ia=Ipk, Ib=Ic=-Ipk/2) faz A aki halkalanmasi = senkron endüktans.
       Kilitli rotor (slip = 1): 50 Hz zaman-harmonik cozum (cubuk akimlari dahil)."""
    Ipk = math.sqrt(2) * I_rms
    if slip < 0.01:
        femm.mi_probdef(0, "millimeters", "planar", 1e-8, L_STACK, 30, 0)
        cur = (("A", Ipk), ("B", -Ipk / 2), ("C", -Ipk / 2))
    else:
        femm.mi_probdef(slip * F_SUPPLY, "millimeters", "planar", 1e-8, L_STACK, 30, 0)
        a = cmath.exp(-2j * math.pi / 3)
        cur = (("A", Ipk), ("B", Ipk * a), ("C", Ipk * a * a))
    for ph, z in cur:
        z = complex(z)
        getattr(femm, "callfemm_noeval", femm.callfemm)('mi_modifycircprop("%s",1,%.10g+I*(%.10g))' % (ph, z.real, z.imag))
    t = time.time()
    femm.mi_analyze(1)
    femm.mi_loadsolution()
    global _shown
    if not _shown:
        print("    ag: %d eleman" % femm.mo_numelements(), flush=True); _shown = True
    print("    (cozum %.0f s)" % (time.time() - t), flush=True)
    I_A, V_A, lam_A = femm.mo_getcircuitproperties("A")
    if slip < 0.01:
        P_bar = 0.0
    else:
        femm.mo_groupselectblock(2)
        P_bar = femm.mo_blockintegral(4)            # tum cubuklar (uc halkasi dahil), W
        femm.mo_clearblock()
    femm.mo_close()
    return lam_A, P_bar

def phase_values(I_rms, lam_A, P_bar, Rs, slip):
    w = 2 * math.pi * F_SUPPLY
    Ipk = math.sqrt(2) * I_rms
    Z = Rs + 1j * w * lam_A / Ipk            # faz empedansi (2D, sargi basi kacagi haric)
    V = abs(Z) * I_rms
    P_in = 3 * Rs * I_rms ** 2 + abs(P_bar.real if isinstance(P_bar, complex) else P_bar)
    return V, P_in, Z

# ----------------------------------------------------------------------------
# 5. ANA AKIS
# ----------------------------------------------------------------------------
if __name__ == "__main__":
    t0 = time.time()
    print("Uc halkasi katsayisi K_ring = %.3f  -> cubuk iletkenligi %.2f MS/m" % (K_RING, SIGMA_BAR_EFF))
    print("Model kuruluyor ...")
    build_model()
    rows = []
    fcsv = open(P("sonuclar_adim2.csv"), "w", newline="")
    wcsv = csv.writer(fcsv)
    wcsv.writerow(["test", "I_rms", "V_olc", "V_hes", "V_fark_%", "P_olc", "P_hes", "R_faz", "X_faz"])
    def add(row):
        rows.append(row); wcsv.writerow(row); fcsv.flush()
    nl, lr = meas["NL"], meas["LR"]
    i_bmp = min(range(len(nl["I"])), key=lambda i: abs(nl["V"][i] - 231))
    print("\nBOSTA (50 Hz, s~0)")
    print("  I_olc[A]  V_olc[V]  V_hes[V]  Fark[%]   Lm_hes[H]")
    for i, (Im, Vm) in enumerate(zip(nl["I"], nl["V"])):
        lam, Pb = run_case(Im, 1e-4)
        V, Pin, Z = phase_values(Im, lam, Pb, nl["Rs"], 1e-4)
        L = Z.imag / (2 * math.pi * F_SUPPLY)
        err = 100 * (V - Vm) / Vm
        print("  %7.3f  %8.1f  %8.1f  %+7.1f   %.4f" % (Im, Vm, V, err, L), flush=True)
        add(["bosta", Im, Vm, round(V, 2), round(err, 2), nl["P"][i], "", round(Z.real, 4), round(Z.imag, 4)])

    print("\nKILITLI ROTOR (50 Hz, s=1)  [2D model: sargi basi kacak enduktansi YOK]")
    print("  I_olc[A]  V_olc[V]  V_hes[V]  Fark[%]   P_olc[W]  P_hes[W]  Fark[%]")
    for i, (Im, Vm, Pm) in enumerate(zip(lr["I"], lr["V"], lr["P"])):
        lam, Pb = run_case(Im, 1.0)
        V, Pin, Z = phase_values(Im, lam, Pb, lr["Rs"], 1.0)
        eV, eP = 100 * (V - Vm) / Vm, 100 * (Pin - Pm) / Pm
        print("  %7.3f  %8.1f  %8.1f  %+7.1f   %8.1f  %8.1f  %+7.1f" % (Im, Vm, V, eV, Pm, Pin, eP), flush=True)
        add(["kilitli", Im, Vm, round(V, 2), round(eV, 2), Pm, round(Pin, 1), round(Z.real, 4), round(Z.imag, 4)])

    fcsv.close()
    femm.closefemm()
    print("\nBitti (%.0f s). Sonuclar: sonuclar_adim2.csv" % (time.time() - t0))
