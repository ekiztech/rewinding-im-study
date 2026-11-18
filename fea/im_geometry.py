# -*- coding: utf-8 -*-
"""
CREATOR asenkron motoru - temiz parametrik 2D geometri
Olculer TU Graz DXF dosyasindan (IM.dxf) tek oluk icin cikarildi, sonra dondurularak cogaltildi.
Tum uzunluklar mm.
"""
import math

N_SLOT, N_BAR = 36, 28
R_OUT, R_BORE, R_ROTOR, R_SHAFT = 100.0, 62.8, 62.4, 22.0
STATOR_SLOT0 = 5.0                    # derece, ilk stator olugunun merkez acisi
ROTOR_SLOT0 = 360.0 / N_BAR / 2       # 6.4286 derece

def rot(p, deg):
    a = math.radians(deg); c, s = math.cos(a), math.sin(a)
    return (p[0] * c - p[1] * s, p[0] * s + p[1] * c)

def mirror(p, axis_deg):
    q = rot(p, -axis_deg); return rot((q[0], -q[1]), axis_deg)

def ang(c, p):
    return math.degrees(math.atan2(p[1] - c[1], p[0] - c[0]))

# --- Stator olugu (merkez 5 derece), alt yari, orijinal koordinatlar ---
S1 = (62.6576, 4.2271)      # delik (bore) uzerinde, agiz kenari
S2 = (62.9644, 4.2539)
C271 = (65.2508, 5.7087)    # R=2.71 ust yuvarlatma merkezi
S3 = (64.2016, 3.2100)      # kama/bobin ust siniri bu yay uzerinde
S4 = (65.2506, 2.9987)
S5 = (72.2506, 2.9987)      # ust bobin alt siniri
S6 = (72.7632, 2.9987)      # alt bobin ust siniri
S7 = (78.7726, 2.9987)
C12 = (78.7337, 4.1981)     # R=1.2 dip kose
S8 = (79.9251, 4.3414)

# --- Rotor olugu (merkez 6.4286 derece), alt yari ---
Q1 = (41.2290, 3.6554)
C099 = (41.2291, 4.6454)    # R=0.99 dip yay (iki yari tek yay)
Q2 = (51.1560, 3.6539)
C211 = (51.1563, 5.7639)
Q3 = (53.2493, 5.4966)
Q4 = (55.8472, 5.7893)
C285 = (58.5793, 6.6003)
Q5 = (61.4235, 6.4176)
Q6_NOM = (62.0616, 6.4895)  # rotor yuzeyinde (r=62.4); cubuk yuzeye kadar uzatildi (0.05 mm sadelestirme)

def ccw_sweep(c, p_from, p_to):
    a = (ang(c, p_to) - ang(c, p_from)) % 360.0
    return a

def build():
    """nodes: liste, segs: [(p,q)], arcs: [(p_from, p_to, sweep_deg_ccw)], labels: [(role, x, y, extra)]"""
    segs, arcs, labels = [], [], []

    # ---------- STATOR ----------
    ax = STATOR_SLOT0
    m = lambda p: mirror(p, ax)
    for k in range(N_SLOT):
        d = 10.0 * k
        R = lambda p: rot(p, d)
        # alt yari
        segs.append((R(S1), R(S2)))
        arcs.append((R(S2), R(S3), ccw_sweep(C271, S2, S3)))
        arcs.append((R(S3), R(S4), ccw_sweep(C271, S3, S4)))
        segs += [(R(S4), R(S5)), (R(S5), R(S6)), (R(S6), R(S7))]
        arcs.append((R(S7), R(S8), ccw_sweep(C12, S7, S8)))
        # ust yari (aynalanmis; yay yonu ters)
        mc271, mc12 = m(C271), m(C12)
        segs.append((R(m(S1)), R(m(S2))))
        arcs.append((R(m(S3)), R(m(S2)), ccw_sweep(mc271, m(S3), m(S2))))
        arcs.append((R(m(S4)), R(m(S3)), ccw_sweep(mc271, m(S4), m(S3))))
        segs += [(R(m(S4)), R(m(S5))), (R(m(S5)), R(m(S6))), (R(m(S6)), R(m(S7)))]
        arcs.append((R(m(S8)), R(m(S7)), ccw_sweep(mc12, m(S8), m(S7))))
        # dip ve ara cizgiler
        segs.append((R(S8), R(m(S8))))
        segs.append((R(S3), R(m(S3))))      # ust bobin ust siniri (kama alti)
        segs.append((R(S5), R(m(S5))))      # ust bobin alt siniri
        segs.append((R(S6), R(m(S6))))      # alt bobin ust siniri
        # delik yayi: bu olugun ust agiz noktasindan bir sonraki olugun alt agiz noktasina
        p_from = R(m(S1)); p_to = rot(S1, 10.0 * ((k + 1) % N_SLOT))
        arcs.append((p_from, p_to, ccw_sweep((0, 0), p_from, p_to)))
        a = math.radians(ax + d)
        for role, r in (("coil_up", 68.0), ("liner", 72.46), ("coil_lo", 76.5)):
            labels.append((role, r * math.cos(a), r * math.sin(a), k))

    # ---------- ROTOR ----------
    # Q6 rotor yuzeyinde; R_ROTOR degisirse (etkin hava araligi kalibrasyonu) radyal olarak olceklenir
    sc = R_ROTOR / 62.4
    Q6 = (Q6_NOM[0] * sc, Q6_NOM[1] * sc)
    ax = ROTOR_SLOT0
    m = lambda p: mirror(p, ax)
    pitch = 360.0 / N_BAR
    for j in range(N_BAR):
        d = pitch * j
        R = lambda p: rot(p, d)
        arcs.append((R(m(Q1)), R(Q1), ccw_sweep(C099, m(Q1), Q1)))      # dip yay
        segs.append((R(Q1), R(Q2)))
        arcs.append((R(Q2), R(Q3), ccw_sweep(C211, Q2, Q3)))
        segs.append((R(Q3), R(Q4)))
        arcs.append((R(Q4), R(Q5), ccw_sweep(C285, Q4, Q5)))
        segs.append((R(Q5), R(Q6)))
        mc211, mc285 = m(C211), m(C285)
        segs.append((R(m(Q1)), R(m(Q2))))
        arcs.append((R(m(Q3)), R(m(Q2)), ccw_sweep(mc211, m(Q3), m(Q2))))
        segs.append((R(m(Q3)), R(m(Q4))))
        arcs.append((R(m(Q5)), R(m(Q4)), ccw_sweep(mc285, m(Q5), m(Q4))))
        segs.append((R(m(Q5)), R(m(Q6))))
        # cubuk ust yayi (rotor yuzeyinde) ve dis arasi yuzey yayi
        arcs.append((R(Q6), R(m(Q6)), ccw_sweep((0, 0), Q6, m(Q6))))
        p_from = R(m(Q6)); p_to = rot(Q6, pitch * ((j + 1) % N_BAR))
        arcs.append((p_from, p_to, ccw_sweep((0, 0), p_from, p_to)))
        a = math.radians(ax + d)
        labels.append(("bar", 50.0 * math.cos(a), 50.0 * math.sin(a), j))

    # ---------- CEMBERLER ----------
    for r in (R_OUT, R_SHAFT):
        arcs.append(((r, 0), (-r, 0), 180.0))
        arcs.append(((-r, 0), (r, 0), 180.0))

    labels += [("stator_iron", 90.0, 0.0, None), ("airgap", 62.6, 0.0, None),
               ("rotor_iron", 30.0, 0.0, None), ("shaft", 0.0, 0.0, None)]

    nodes = []
    for p, q in segs:
        nodes += [p, q]
    for p, q, _ in arcs:
        nodes += [p, q]
    # benzersiz dugumler (1e-6 mm)
    uniq = {}
    for p in nodes:
        uniq[(round(p[0], 6), round(p[1], 6))] = p
    return list(uniq.values()), segs, arcs, labels
