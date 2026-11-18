"""Hybrid model v2: de-embedded (true) circuit parameters from the calibrated 2D FEA, three coil pitches."""
import math, cmath, json
import numpy as np
D = json.load(open('deembedded2.json'))
F0, P_POLE, V_PH = 50.0, 2, 400.0 / math.sqrt(3)
W_S = 2 * math.pi * F0
L3D_S, L3D_R = 2.91e-3, 3.15e-3                      # 3D end effects split from measured vs. 2D
KP = {7: math.sin(math.radians(70)), 8: math.sin(math.radians(80)), 9: 1.0}
MLT = {7: 0.553, 8: 0.553 + 0.0297, 9: 0.553 + 2 * 0.0297}   # end-turn arc grows with coil span (assumption)
ALPHA_CU, ALPHA_AL, RHO20 = 0.00393, 0.0039, 1.724e-8
N0 = 18; NPH0 = 216
A_SLOT = 2 * N0 * 2 * math.pi * (0.75e-3) ** 2 / 4 / 0.30

class Pitch:
    def __init__(self, p):
        d = D[str(p)]
        self.p = p
        self.Lss = d['Ls'] + L3D_S
        nl = sorted(d['nl']); self.I = np.array([a for a, b in nl]); self.Ls_nl = np.array([b for a, b in nl])
        self.psi = (self.Ls_nl - self.Lss) * self.I            # magnetizing flux linkage vs current (18 turns)
        lr = sorted(d['lr']); self.f = np.array([x[0] for x in lr]); self.Rr = np.array([x[1] for x in lr])
        self.Lsr = np.array([x[2] for x in lr]) + L3D_R
    def psi_m(self, I0):
        if I0 <= self.I[0]: return self.psi[0] / self.I[0] * I0
        if I0 >= self.I[-1]:
            sl = (self.psi[-1] - self.psi[-2]) / (self.I[-1] - self.I[-2]); return self.psi[-1] + sl * (I0 - self.I[-1])
        return float(np.interp(I0, self.I, self.psi))
PITCH = {p: Pitch(p) for p in (7, 8, 9)}

PFE_C, PFW_F, PLL_F = 1.0, 1.0, 1.0
def pfe(psi_phys, k_fe=1.0): return PFE_C * k_fe * max(252.0 * psi_phys ** 2 - 6.3, 0.0)
def pfw(n): return PFW_F * 161.0 * (n / 1500.0) ** 2

def solve(s, Rs, pitch=7, kN=1.0, Tr=75.0, k_fe=1.0, fr=None, V=V_PH):
    P = PITCH[pitch]; fr = s * F0 if fr is None else fr
    Xs = W_S * P.Lss * kN ** 2
    Xr = W_S * float(np.interp(fr, P.f, P.Lsr)) * kN ** 2
    Rr = float(np.interp(fr, P.f, P.Rr)) * (1 + ALPHA_AL * (Tr - 20)) * kN ** 2
    Zr = Rr / s + 1j * Xr
    Lm = 0.25 * kN ** 2
    for _ in range(300):
        Zm = 1j * W_S * Lm; Zp = 1 / (1 / Zm + 1 / Zr)
        I1 = V / (Rs + 1j * Xs + Zp); E = I1 * Zp; Im = abs(E) / (W_S * Lm)
        Lm_new = kN * P.psi_m(kN * Im) / Im
        if abs(Lm_new - Lm) < 1e-8: Lm = Lm_new; break
        Lm = 0.5 * (Lm + Lm_new)
    I2 = E / Zr; psi = abs(E) / W_S
    Pfe = pfe(psi / (kN * KP[pitch] / KP[7]), k_fe)   # physical flux relative to original winding
    Pcu1 = 3 * abs(I1) ** 2 * Rs; Pag = 3 * abs(I2) ** 2 * Rr / s; Pcu2 = s * Pag
    n = (1 - s) * 60 * F0 / P_POLE
    P1 = 3 * (V * I1.conjugate()).real + Pfe
    Pm = Pag - Pcu2; Pf = pfw(n)
    PLL = PLL_F * P1 * (0.025 - 0.005 * math.log10(max(Pm - Pf, 1) / 1000))
    P2 = Pm - Pf - PLL
    T2 = P2 / (2 * math.pi * n / 60) if n > 0 else Pag / (W_S / P_POLE)
    return dict(s=s, n=n, I1=abs(I1), pf=math.cos(cmath.phase(I1)), P1=P1, P2=P2, eta=P2 / P1, T2=T2,
                Tem=Pag / (W_S / P_POLE), P_cu1=Pcu1, P_cu2=Pcu2, P_fe=Pfe, P_fw=Pf, P_LL=PLL, psi=psi)

def at_power(P2t, Rs, **kw):
    lo, hi = 1e-4, 0.3
    for _ in range(70):
        m = 0.5 * (lo + hi)
        if solve(m, Rs, **kw)['P2'] < P2t: lo = m
        else: hi = m
    return solve(0.5 * (lo + hi), Rs, **kw)

def Rs_of(Nl, nstr, d_mm, pitch=7, mlt_factor=1.0, T=75.0):
    A = nstr * math.pi * (d_mm * 1e-3) ** 2 / 4
    return RHO20 * (1 + ALPHA_CU * (T - 20)) * (Nl * NPH0 / N0) * MLT[pitch] * mlt_factor / A
def fill_of(Nl, nstr, d_mm): return 2 * Nl * nstr * math.pi * (d_mm * 1e-3) ** 2 / 4 / A_SLOT
def start_bd(Rs, **kw):
    st = solve(1.0, Rs, fr=50.0, **kw)
    bd = max(solve(s, Rs, **kw)['Tem'] for s in np.linspace(0.04, 0.45, 42))
    return st['I1'], st['Tem'], bd

if __name__ == '__main__':
    for p in (7, 8, 9):
        P = PITCH[p]; print('pitch %d/9: Lss %.2f mH, Lsr(2Hz) %.2f mH, Rr(2Hz) %.3f, Lm_lin %.4f H, MLT %.3f' % (p, P.Lss*1e3, P.Lsr[0]*1e3, P.Rr[0], P.psi[0]/P.I[0], MLT[p]))
    Rs = Rs_of(18, 2, 0.75)
    for T in (20, 75):
        r = at_power(3700, Rs_of(18,2,0.75,T=T), Tr=T)
        print('T=%d: n=%.0f I=%.2f pf=%.3f eta=%.2f%% | Pcu1 %.0f Pcu2 %.0f Pfe %.0f Pfw %.0f PLL %.0f' % (T, r['n'], r['I1'], r['pf'], 100*r['eta'], r['P_cu1'], r['P_cu2'], r['P_fe'], r['P_fw'], r['P_LL']))
    Ist, Tst, Tbd = start_bd(Rs); print('start I %.1f A, T_start %.1f Nm, T_bd %.1f Nm' % (Ist, Tst, Tbd))
