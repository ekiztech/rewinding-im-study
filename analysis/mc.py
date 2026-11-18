import json, math, numpy as np, copy, time, sys
import hybrid2 as H
rng = np.random.default_rng(2026)
P7 = H.PITCH[7]; ORIG = dict(psi=P7.psi.copy(), Rr=P7.Rr.copy(), Lsr=P7.Lsr.copy(), Lss=P7.Lss)
D7 = H.D['7']; L2D_s = D7['Ls']; L2D_r = np.array([x[2] for x in sorted(D7['lr'])])
L3D_TOT = H.L3D_S + H.L3D_R
def set_sample(p):
    P7.psi = ORIG['psi'] * p['cm']
    P7.Rr = ORIG['Rr'] * p['cr']
    l3s = L3D_TOT * p['split']; l3r = L3D_TOT - l3s
    P7.Lss = L2D_s * p['cl'] + l3s
    P7.Lsr = L2D_r * p['cl'] + l3r
    # magnetizing curve was derived with the nominal Lss; keep psi_m(I) unchanged except scale (stator leakage change is absorbed)
    H.PFE_C, H.PFW_F, H.PLL_F = p['cfe'], p['cfw'], p['cll']
def sample():
    return dict(cm=rng.uniform(0.97, 1.03), cr=rng.normal(1, 0.05), cl=rng.uniform(0.9, 1.1), split=rng.uniform(0.3, 0.7),
                cfe=rng.normal(1, 0.10), cfw=rng.normal(1, 0.15), cll=rng.uniform(0.5, 1.5), T=rng.uniform(60, 95))
def req_area(dk, mlt, T, eta0, Rs0):
    lo, hi = 0.9, 1.8
    for _ in range(22):
        m = 0.5 * (lo + hi)
        if H.at_power(3700, Rs0 * mlt / m, Tr=T, k_fe=1 + dk)['eta'] < eta0: lo = m
        else: hi = m
    return 0.5 * (lo + hi)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
out = []; t0 = time.time()
for i in range(N):
    p = sample(); set_sample(p); T = p['T']
    Rs0 = H.Rs_of(18, 2, 0.75, T=T)
    b = H.at_power(3700, Rs0, Tr=T); eta0 = b['eta']
    r = dict(p, eta0=eta0, ratio=b['P_fe'] / b['P_cu1'])
    for dk in (0.2, 0.5, 1.0):
        for mlt in (1.0, 1.05):
            r['a_%g_%g' % (dk, mlt)] = req_area(dk, mlt, T, eta0, Rs0)
    for lab, (n, d) in (('d285', (2, 0.85)), ('d375', (3, 0.75))):
        for dmg, k in (('D0', 1.0), ('D1', 1.2)):
            r['deta_%s_%s' % (lab, dmg)] = 100 * (H.at_power(3700, H.Rs_of(18, n, d, T=T), Tr=T, k_fe=k)['eta'] - eta0)
    out.append(r)
json.dump(out, open('mc.json', 'w'))
print('N=%d in %.0f s' % (N, time.time() - t0))
keys = [k for k in out[0] if k.startswith('a_') or k.startswith('deta') or k in ('eta0', 'ratio')]
for k in keys:
    v = np.array([o[k] for o in out]); print('%-14s median %.4f  95%% [%.4f, %.4f]' % (k, np.median(v), np.percentile(v, 2.5), np.percentile(v, 97.5)))
