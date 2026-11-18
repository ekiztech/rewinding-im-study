# -*- coding: utf-8 -*-
"""
Adim 9 - Hakem revizyonu icin ek FEA kosulari
=============================================
Bolum A (~1-1.5 saat): Kilitli rotor, 50 Hz, uc bobin adimi (7/9, 8/9, 9/9)
   - Akimlar: 4.622 A (olcum seviyesi), 15 A ve nominal gerilimdeki kalkis akimi (doymus)
   - Ayrica 2 Hz, 4.622 A (nominal kaymaya yakin)
   - Her cozumde manyetik enerji: stator bobinleri, katman ayiricisi, stator oluk agizlari
     (nokta ornekleme), rotor cubuklari, hava araligi -> stator/rotor kacak paylasimi FEA'dan
   Cikti: sonuclar_adim9A.csv
Bolum B (~2.5-3 saat): Nominal yukte rotoru adim adim dondurerek manyetostatik cozumler
   - Stator akimlari ve temel rotor cubuk akimlari nominal yukteki degerlerinde (hibrit modelden)
   - 1/6 elektriksel periyot (30 derece mekanik) 24 adimda taranir
   - Her adimda her cubugun ortalama vektor potansiyeli kaydedilir -> harmonik cubuk EMK'lari
     ve harmonik kafes kayiplari (analizi Claude yapar)
   Cikti: sonuclar_adim9B2.csv, adim9B_kalibrasyon.json

Calistirma:
   py adim9_ek_fea.py A      (sadece Bolum A)
   py adim9_ek_fea.py B      (sadece Bolum B)
   py adim9_ek_fea.py        (ikisi de, ~4 saat)
Kesilirse ayni komutla tekrar baslat: biten satirlar atlanir, kaldigi yerden devam eder.
"""
import os, sys, json, math, cmath, csv, time
import femm
import im_geometry as geo
import adim2_dogrulama as A

geo.R_ROTOR = 62.8 - 0.44
A.LAM_FILL = 0.87
A.SIGMA_AL = 29.3
A.SIGMA_BAR_EFF = 29.3 / A.K_RING
ROTOR0 = geo.ROTOR_SLOT0
MU0 = 4e-7 * math.pi
CALL = getattr(femm, "callfemm_noeval", femm.callfemm)

# Hibrit modelden nominal yuk (3.7 kW, 18 sarim, 2x0.75 mm, 75 C); fazorlar faz gerilimine gore (derece)
RATED = {7: dict(I1=7.7029, a1=-29.78, I2=6.9015, a2=-9.67, Ist=30.91),
         8: dict(I1=7.7118, a1=-28.78, I2=7.0102, a2=-10.77, Ist=28.62),
         9: dict(I1=7.8651, a1=-29.31, I2=7.1836, a2=-12.76, Ist=26.03)}


def done_keys(fname, ncols):
    keys = set()
    if os.path.exists(A.P(fname)):
        with open(A.P(fname)) as f:
            r = csv.reader(f); next(r, None)
            for row in r:
                if len(row) >= ncols: keys.add(tuple(row[:ncols]))
    return keys


def open_csv(fname, header):
    new = not os.path.exists(A.P(fname))
    f = open(A.P(fname), "a", newline=""); w = csv.writer(f)
    if new: w.writerow(header); f.flush()
    return f, w


def regroup(labels, bar_circuits=False):
    """Bobinler grup 3, katman ayiricisi grup 4; istenirse her cubuga ayri devre (b0..b27)."""
    if bar_circuits:
        for j in range(geo.N_BAR):
            femm.mi_addcircprop("b%d" % j, 0, 1)
    for role, x, y, k in labels:
        if role in ("coil_up", "coil_lo", "liner") or (bar_circuits and role == "bar"):
            femm.mi_selectlabel(x, y)
            if role == "liner":
                femm.mi_setblockprop("Air", 0, 0.4, "<None>", 0, 4, 0)
            elif role == "bar":
                femm.mi_setblockprop("Al_bar", 0, 0.7, "b%d" % k, 0, 2, 1)
            else:
                ph, sg = A.upper_layer(k) if role == "coil_up" else A.lower_layer(k)
                femm.mi_setblockprop("Coil", 0, 1.5, ph, 0, 3, sg * A.TURNS_LAYER)
            femm.mi_clearselected()


def local_to_global(xl, yl, deg):
    a = math.radians(deg); c, s = math.cos(a), math.sin(a)
    return xl * c - yl * s, xl * s + yl * c


def sampled_energy(x0, x1, nx, hy, ny, slots, harmonic=True):
    """Oluk yerel koordinatlarinda dikdortgen icinde, sadece hava/bakir (mu<1.5) noktalarda enerji [J]."""
    dx, dy = (x1 - x0) / nx, 2 * hy / ny
    W = 0.0
    for k in slots:
        deg = geo.STATOR_SLOT0 + 10.0 * k
        for i in range(nx):
            xl = x0 + (i + 0.5) * dx
            for j in range(ny):
                yl = -hy + (j + 0.5) * dy
                x, y = local_to_global(xl, yl, deg)
                v = femm.mo_getpointvalues(x, y)
                if not v: continue
                mu = abs(complex(v[9]))
                if mu > 1.5: continue
                b2 = abs(complex(v[1])) ** 2 + abs(complex(v[2])) ** 2
                W += b2 / ((4 if harmonic else 2) * MU0) * dx * dy * 1e-6 * A.L_STACK * 1e-3
    return W


def block_energy(x, y):
    femm.mo_selectblock(x, y); e = femm.mo_blockintegral(2); femm.mo_clearblock(); return e


def group_integral(g, kind):
    femm.mo_groupselectblock(g); v = femm.mo_blockintegral(kind); femm.mo_clearblock(); return v


# ============================================================ BOLUM A
def part_A():
    hdr = ["pitch_slots", "f_Hz", "I_rms", "R_app", "X", "W_coils", "W_sep", "W_open", "W_bars",
           "W_airgap", "W_up0_block", "W_up0_sampled", "t_s"]
    done = done_keys("sonuclar_adim9A.csv", 3)
    f, w = open_csv("sonuclar_adim9A.csv", hdr)
    for pitch in (7, 8, 9):
        cases = [(50.0, 4.622), (50.0, 15.0), (50.0, RATED[pitch]["Ist"]), (2.0, 4.622)]
        todo = [c for c in cases if (str(pitch), str(c[0]), str(c[1])) not in done]
        if not todo: continue
        A.PITCH_SLOTS = pitch; A._shown = False
        nodes, segs, arcs, labels = geo.build()
        A.build_model(); regroup(labels)
        lab = {(r, k): (x, y) for r, x, y, k in labels}
        print("\nBOLUM A - adim %d/9" % pitch, flush=True)
        for fr, I in todo:
            t = time.time()
            A.F_SUPPLY = fr
            Ipk = math.sqrt(2) * I
            femm.mi_probdef(fr, "millimeters", "planar", 1e-8, A.L_STACK, 30, 0)
            a = cmath.exp(-2j * math.pi / 3)
            for ph, z in (("A", Ipk), ("B", Ipk * a), ("C", Ipk * a * a)):
                z = complex(z); CALL('mi_modifycircprop("%s",1,%.10g+I*(%.10g))' % (ph, z.real, z.imag))
            femm.mi_analyze(1); femm.mi_loadsolution()
            _, _, lam = femm.mo_getcircuitproperties("A")
            Pb = abs(complex(group_integral(2, 4)))
            Z = 1j * 2 * math.pi * fr * complex(lam) / Ipk
            W_coils = abs(complex(group_integral(3, 2)))
            W_sep = abs(complex(group_integral(4, 2)))
            W_bars = abs(complex(group_integral(2, 2)))
            W_gap = abs(complex(block_energy(62.6, 0.0)))
            W_open = sampled_energy(62.81, 64.15, 8, 2.6, 14, range(geo.N_SLOT))
            xu, yu = lab[("coil_up", 0)]
            W_up_b = abs(complex(block_energy(xu, yu)))
            W_up_s = sampled_energy(64.35, 72.15, 16, 3.4, 14, [0])
            femm.mo_close()
            dt = time.time() - t
            R = Pb / (3 * I * I)
            print("  f=%4.1f Hz I=%6.2f A  R=%.4f X=%.4f  W[coil %.3g, open %.3g, bar %.3g]  (%.0f s)"
                  % (fr, I, R, Z.imag, W_coils, W_open, W_bars, dt), flush=True)
            w.writerow([pitch, fr, I, round(R, 6), round(Z.imag, 6), W_coils, W_sep, W_open, W_bars,
                        W_gap, W_up_b, W_up_s, round(dt)]); f.flush()
        femm.closefemm()
    f.close()


# ============================================================ BOLUM B
R_MID = 62.8 - 0.22
NSTEP = 24


def b1_complex(npts=720):
    c = 0j
    for n in range(npts):
        th = 2 * math.pi * n / npts
        bx, by = femm.mo_getb(R_MID * math.cos(th), R_MID * math.sin(th))
        br = bx * math.cos(th) + by * math.sin(th)
        c += br * cmath.exp(-2j * th)
    return 2 * c / npts


def set_stator(Ipk, phi_deg):
    for k, ph in enumerate("ABC"):
        CALL('mi_modifycircprop("%s",1,%.10g)' % (ph, Ipk * math.cos(math.radians(phi_deg - 120 * k))))


def set_bars(Ib, delta_deg, rot_deg):
    for j in range(geo.N_BAR):
        th = math.radians(ROTOR0 + rot_deg + 360.0 / geo.N_BAR * j)
        CALL('mi_modifycircprop("b%d",1,%.10g)' % (j, Ib * math.cos(2 * th - math.radians(delta_deg))))


_ROT = [0.0]
_orig_build = geo.build


def _build_rotated():
    """Standart geometriyi kurar, sonra rotor tarafini (r < 62.5 mm) butun olarak dondurur.
    (ROTOR_SLOT0 degistirmek oluk seklini bozuyordu: Q noktalari donmuyor, sadece ayna ekseni kayiyordu.)"""
    nodes, segs, arcs, labels = _orig_build()
    d = _ROT[0]
    if abs(d) < 1e-12:
        return nodes, segs, arcs, labels
    f = lambda p: geo.rot(p, d) if math.hypot(p[0], p[1]) < 62.5 else p
    segs = [(f(p), f(q)) for p, q in segs]
    arcs = [(f(p), f(q), sw) for p, q, sw in arcs]
    labels = [(r,) + f((x, y)) + (k,) for r, x, y, k in labels]
    uniq = {}
    for p, q in segs:
        for t in (p, q): uniq[(round(t[0], 6), round(t[1], 6))] = t
    for p, q, _ in arcs:
        for t in (p, q): uniq[(round(t[0], 6), round(t[1], 6))] = t
    return list(uniq.values()), segs, arcs, labels


geo.build = _build_rotated


def build_B(pitch, rot_deg):
    _ROT[0] = rot_deg
    A.PITCH_SLOTS = pitch
    nodes, segs, arcs, labels = geo.build()
    A.build_model(); regroup(labels, bar_circuits=True)
    femm.mi_probdef(0, "millimeters", "planar", 1e-8, A.L_STACK, 30, 0)
    return [(x, y) for r, x, y, k in sorted([l for l in labels if l[0] == "bar"], key=lambda l: l[3])]


def solve():
    femm.mi_analyze(1); femm.mi_loadsolution()


def calibrate(pitch):
    fn = A.P("adim9B_kalibrasyon.json")
    cal = json.load(open(fn)) if os.path.exists(fn) else {}
    if str(pitch) in cal: return cal[str(pitch)]
    print("  kalibrasyon (4 dogrusal cozum) ...", flush=True)
    build_B(pitch, 0.0)
    Ic, Ibc = 3.0, 100.0
    res = {}
    for name, st, br in (("s0", (Ic, 0), None), ("s90", (Ic, 90), None),
                         ("b0", None, (Ibc, 0)), ("b90", None, (Ibc, 90))):
        set_stator(*(st if st else (0.0, 0)))
        set_bars(*(br if br else (0.0, 0)), rot_deg=0.0)
        solve(); c = b1_complex(); femm.mo_close()
        res[name] = [c.real, c.imag]
    femm.closefemm()
    c = {k: complex(*v) for k, v in res.items()}
    s_s = 1 if abs(cmath.phase(c["s90"] / c["s0"]) - math.pi / 2) < math.pi / 4 else -1
    s_b = 1 if abs(cmath.phase(c["b90"] / c["b0"]) - math.pi / 2) < math.pi / 4 else -1
    r = RATED[pitch]
    Ib = Ibc * abs(c["s0"]) / abs(c["b0"]) * math.sqrt(2) * r["I2"] / Ic
    arg_mI2 = r["a2"] + 180.0
    d0 = s_b * math.degrees(cmath.phase(c["s0"]) - cmath.phase(c["b0"]) + s_s * math.radians(arg_mI2))
    out = dict(raw=res, s_s=s_s, s_b=s_b, Ib_pk=Ib, delta0_deg=d0)
    cal[str(pitch)] = out; json.dump(cal, open(fn, "w"), indent=1)
    print("  s_s=%d s_b=%d  cubuk tepe akimi %.1f A  delta0=%.1f derece" % (s_s, s_b, Ib, d0), flush=True)
    return out


def part_B():
    hdr = ["pitch_slots", "step", "tau_deg", "rot_deg", "B1_abs", "lamA", "lamB", "lamC"] + \
          ["Abar%d" % j for j in range(geo.N_BAR)]
    done = done_keys("sonuclar_adim9B2.csv", 2)
    f, w = open_csv("sonuclar_adim9B2.csv", hdr)
    for pitch in (7, 8, 9):
        if all((str(pitch), str(j)) in done for j in range(NSTEP)): continue
        print("\nBOLUM B - adim %d/9" % pitch, flush=True)
        cal = calibrate(pitch); r = RATED[pitch]
        for j in range(NSTEP):
            if (str(pitch), str(j)) in done: continue
            t = time.time()
            tau = 60.0 * j / NSTEP                        # elektriksel derece
            rot = -cal["s_s"] * tau / 2.0                 # mekanik derece, alan yonunde
            bars = build_B(pitch, rot)
            set_stator(math.sqrt(2) * r["I1"], tau + r["a1"])
            delta = cal["delta0_deg"] + cal["s_b"] * cal["s_s"] * tau
            set_bars(cal["Ib_pk"], delta, rot)
            solve()
            lam = [complex(femm.mo_getcircuitproperties(ph)[2]).real for ph in "ABC"]
            b1 = abs(b1_complex(360))
            Ab = []
            for x, y in bars:
                femm.mo_selectblock(x, y)
                Ab.append(femm.mo_blockintegral(1) / femm.mo_blockintegral(5))
                femm.mo_clearblock()
            femm.mo_close(); femm.closefemm()
            w.writerow([pitch, j, round(tau, 4), round(rot, 5), round(b1, 6)] + lam + Ab); f.flush()
            print("  adim %2d/%d  B1=%.3f T  (%.0f s)" % (j + 1, NSTEP, b1, time.time() - t), flush=True)
    f.close()
    _ROT[0] = 0.0


if __name__ == "__main__":
    t0 = time.time()
    which = sys.argv[1].upper() if len(sys.argv) > 1 else "AB"
    if "A" in which: part_A()
    if "B" in which: part_B()
    print("\nBitti (%.0f dk)." % ((time.time() - t0) / 60))
