"""Cage losses from stator MMF space harmonics (classical harmonic equivalent circuit, per Alger)."""
import math, json, numpy as np
import hybrid2 as H
def kw(p, nu): return math.sin(nu * p / 9 * math.pi / 2) * math.sin(nu * math.pi / 6) / (3 * math.sin(nu * math.pi / 18))
def ksk(nu, skew_el_deg):
    if skew_el_deg == 0: return 1.0
    x = math.radians(nu * skew_el_deg) / 2; return math.sin(x) / x
W = 2 * math.pi * 50
Lm = 0.245; Rr50 = 1.908; Lr50 = 9.38e-3       # measured, de-embedded, 20 C
def harm_loss(pitch, I1, s, skew_el=0.0, T=75.0):
    kw1 = kw(pitch, 1); tot = 0; per = {}
    for nu in (-5, 7, -11, 13, -17, 19):      # sign = rotation direction
        a = abs(nu)
        if a in (17, 19): continue            # slot harmonics: same for all pitches, excluded from the pitch comparison
        Xm_nu = W * Lm * (kw(pitch, a) / (a * kw1)) ** 2 * ksk(a, skew_el) ** 2
        s_nu = 1 - nu * (1 - s)
        fr = abs(s_nu) * 50
        Rr = Rr50 * math.sqrt(fr / 50) * (1 + 0.0039 * (T - 20))   # deep-bar extrapolation ~ sqrt(f)
        Xr = W * Lr50 * (50 / fr) ** 0.5 * 1.0                       # leakage falls with skin effect
        Xr_nu = Xr / a ** 0                                          # referred to the harmonic field (approx.)
        Zr = Rr / s_nu + 1j * Xr_nu * (W * abs(s_nu) * 0 + 1)        # stator-frequency impedance view
        Zm = 1j * Xm_nu
        Zp = Zr * Zm / (Zr + Zm)
        P = 3 * I1 ** 2 * Zp.real            # air-gap power of harmonic ...
        Pcu = abs(P * s_nu)                  # ... rotor copper loss = s_nu * P_ag,nu
        per[nu] = Pcu; tot += Pcu
    return tot, per
Rs = H.Rs_of(18, 2, 0.75)
res = {}
for sk in (0.0, 20.0):
    print('skew = %s' % ('none' if sk == 0 else 'one stator slot pitch (20 deg el.)'))
    for p in (7, 8, 9):
        Nl = 18 if p == 7 else 17
        r = H.at_power(3700, H.Rs_of(Nl, 2, 0.75, pitch=p), pitch=p, kN=Nl / 18)
        tot, per = harm_loss(p, r['I1'], r['s'], sk)
        res['%d_%d' % (p, sk)] = tot
        print('  pitch %d/9 (N=%d): I=%.2f A, belt-harmonic cage loss = %.1f W  ' % (p, Nl, r['I1'], tot) + ' '.join('nu%+d:%.1f' % (k, v) for k, v in per.items()))
json.dump(res, open('harmloss.json', 'w'))
