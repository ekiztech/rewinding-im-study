"""FEA run 9B: harmonic cage losses at rated load from rotating-rotor magnetostatic solutions.
Per bar, the mean vector potential over 1/6 electrical period (rotor at synchronous speed, fundamental bar
currents imposed) gives the open-circuit harmonic bar EMFs. A 2D DFT (time x bars) separates them into
(rotor frequency, cage mode) pairs; each mode drives bar currents through the cage impedance
R_b,ac + R_e/(2 sin^2(pi m/N_b)) + j w (L_cage(m) + L_b), with L_cage(m) the air-gap (differential) inductance of
the cage mode summed over all field orders m + k N_b."""
import csv, json, math, sys, numpy as np
sys.path.insert(0, '../fea'); import im_geometry as g
NB, L, MU0, F = 28, 0.100, 4e-7*math.pi, 50.0
R_ROT, G_EFF = 0.06236, 0.44e-3
SIG = 29.3e6/(1 + 0.0039*55)                               # Al bar conductivity at 75 C
# bar cross-section from the slot outline (polygon with arcs sampled)
def arcpts(c, p, q, n=24, ccw=True):
    a0 = math.atan2(p[1]-c[1], p[0]-c[0]); a1 = math.atan2(q[1]-c[1], q[0]-c[0]); r = math.hypot(p[0]-c[0], p[1]-c[1])
    d = (a1 - a0) % (2*math.pi) if ccw else -((a0 - a1) % (2*math.pi))
    return [(c[0]+r*math.cos(a0+d*i/n), c[1]+r*math.sin(a0+d*i/n)) for i in range(n+1)]
m_ = lambda p: g.mirror(p, g.ROTOR_SLOT0); Q6 = g.Q6_NOM
half = arcpts(g.C099, m_(g.Q1), g.Q1)[len(arcpts(g.C099, m_(g.Q1), g.Q1))//2:] + [g.Q2] + arcpts(g.C211, g.Q2, g.Q3) + [g.Q4] + arcpts(g.C285, g.Q4, g.Q5) + [Q6]
poly = half + [m_(p) for p in reversed(half)]
A_BAR = 0.5*abs(sum(poly[i][0]*poly[(i+1)%len(poly)][1] - poly[(i+1)%len(poly)][0]*poly[i][1] for i in range(len(poly))))*1e-6
H_BAR = (math.hypot(*Q6) - (math.hypot(*g.C099) - 0.99))*1e-3
R_B = L/(SIG*A_BAR); K_RING = 1.263
R_E = (K_RING - 1)*2*math.sin(math.pi*2/NB)**2*R_B           # one end-ring segment, from K_ring of Eq. (1)
def r_ac(f):                                                  # deep-bar AC resistance, simple skin-depth model
    d = math.sqrt(2/(2*math.pi*f*MU0*SIG)); return R_B*max(1.0, H_BAR/d)
def l_bar(f):                                                 # bar slot leakage, only the skin layer carries current
    d = math.sqrt(2/(2*math.pi*f*MU0*SIG)); w = A_BAR/H_BAR; h = min(H_BAR, d)
    return MU0*L*h/(3*w)
def l_cage(m, carter):
    gpp = G_EFF*carter; s = sum(1.0/(m + k*NB)**2 for k in range(-40, 41) if m + k*NB != 0)
    return MU0*NB*R_ROT*L/(2*math.pi*gpp)*s
# radial width profile of the bar -> scale factor between bar-mean potential and the potential in the skin layer
from matplotlib.path import Path as _P
_path = _P(poly); _ax = math.radians(g.ROTOR_SLOT0)
_rs = np.linspace(40.0, 62.4, 450); _w = []
for r in _rs:
    th = np.linspace(_ax - 0.12, _ax + 0.12, 400); pts = np.c_[r*np.cos(th), r*np.sin(th)]
    _w.append(_path.contains_points(pts).sum()*(0.24/400)*r)
_w = np.array(_w)
def depth_factor(n, f):
    d = 1e3*math.sqrt(2/(2*math.pi*f*MU0*SIG)); prof = (_rs/62.4)**abs(n)
    mean_bar = np.sum(prof*_w)/np.sum(_w); top = (_rs >= 62.4 - d)
    return float(np.sum(prof[top]*_w[top])/np.sum(_w[top])/mean_bar)
def analyse(fname='../results/sonuclar_adim9B2.csv', carters=(1.1, 1.25, 1.4), skin=False):
    rows = list(csv.DictReader(open(fname))); out = {}
    for p in sorted({int(r['pitch_slots']) for r in rows}):
        R = sorted([r for r in rows if int(r['pitch_slots']) == p], key=lambda r: int(r['step']))
        N = len(R)
        if N < 4: continue
        A = np.array([[float(r['Abar%d' % k]) for k in range(NB)] for r in R])      # (time, bar)
        X = np.fft.fft2(A)/(N*NB)
        res = {'steps': N, 'modes': []}
        for c in carters:
            P = 0.0
            for q in range(1, N//2 + 1):
                f = 6*F*q; w = 2*math.pi*f
                for mi in range(NB):
                    a = 2*abs(X[q, mi]) if q < N/2 else abs(X[q, mi])               # A amplitude of that travelling wave
                    if a < 1e-12: continue
                    m = mi if mi <= NB//2 else mi - NB
                    if m == 0: continue
                    E = w*L*a*(depth_factor(m, f) if skin else 1.0)
                    Z = r_ac(f) + R_E/(2*math.sin(math.pi*m/NB)**2) + 1j*w*(l_cage(abs(m), c) + l_bar(f))
                    I = E/abs(Z); Pm = NB*0.5*I*I*Z.real
                    P += Pm
                    if c == carters[1] and Pm > 0.01: res['modes'].append(dict(f=f, m=m, E=E, I=I, P=Pm))
            res['P_carter%.2f' % c] = P
        A1 = 2*abs(X[0, 2]) if True else 0
        res['A1_fund'] = A1
        res['modes'] = sorted(res['modes'], key=lambda d: -d['P'])[:8]
        out[p] = res
        print('pitch %d/9 (%d steps): harmonic cage loss %s W' % (p, N, ' / '.join('%.1f' % res['P_carter%.2f' % c] for c in carters)),
              '| top modes:', ', '.join('%.0f Hz m=%+d %.1f W' % (d['f'], d['m'], d['P']) for d in res['modes'][:4]))
    return out
if __name__ == '__main__':
    print('bar area %.1f mm2, depth %.1f mm, R_b %.1f uOhm, R_e %.2f uOhm' % (A_BAR*1e6, H_BAR*1e3, R_B*1e6, R_E*1e6))
    fn = sys.argv[1] if len(sys.argv) > 1 else '../results/sonuclar_adim9B2.csv'
    r0 = analyse(fn); print('-- with skin-layer EMF (exciting field order = cage mode, upper estimate) --'); r1 = analyse(fn, skin=True)
    json.dump({'bar_mean': r0, 'skin_layer': r1}, open('fea9b.json', 'w'), indent=1, default=float)
