"""Monte Carlo v2: parameter uncertainty (as v1) + model-form uncertainty (stator-copper-loss level / rated-current mismatch,
temperature 60-115 C, magnetizing curve recomputed with perturbed leakage) + shop-side uncertainty of the core-loss increase
(loop-test reading vs. whole-core loss, factor 0.7-1.3)."""
import json, math, numpy as np, time, sys
import hybrid2 as H
rng = np.random.default_rng(2026)
P7 = H.PITCH[7]; ORIG = dict(Rr=P7.Rr.copy(), Lsr=P7.Lsr.copy())
D7 = H.D['7']; L2D_s = D7['Ls']; L2D_r = np.array([x[2] for x in sorted(D7['lr'])]); L3D_TOT = H.L3D_S + H.L3D_R
def set_sample(p):
    l3s = L3D_TOT*p['split']; l3r = L3D_TOT - l3s
    P7.Lss = L2D_s*p['cl'] + l3s; P7.Lsr = L2D_r*p['cl'] + l3r
    P7.psi = (P7.Ls_nl - P7.Lss)*P7.I*p['cm']                 # magnetizing curve consistent with perturbed leakage
    P7.Rr = ORIG['Rr']*p['cr']; H.PFE_C, H.PFW_F, H.PLL_F = p['cfe'], p['cfw'], p['cll']
def sample():
    return dict(cm=rng.uniform(0.97,1.03), cr=rng.normal(1,0.05), cl=rng.uniform(0.9,1.1), split=rng.uniform(0.3,0.7),
                cfe=rng.normal(1,0.10), cfw=rng.normal(1,0.15), cll=rng.uniform(0.5,1.5), T=rng.uniform(60,115),
                cI=rng.uniform((6.9/7.7)**2, 1.0), kloop=rng.uniform(0.7,1.3))
def req(dk, mlt, T, eta0, Rs0):
    lo, hi = 0.9, 2.0
    for _ in range(24):
        m = 0.5*(lo+hi)
        if H.at_power(3700, Rs0*mlt/m, Tr=T, k_fe=1+dk)['eta'] < eta0: lo = m
        else: hi = m
    return 0.5*(lo+hi)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 500
out = []; t0 = time.time()
for i in range(N):
    p = sample(); set_sample(p); T = p['T']; Rs0 = H.Rs_of(18,2,0.75,T=T)*p['cI']
    b = H.at_power(3700, Rs0, Tr=T); e0 = b['eta']
    r = dict(p, eta0=e0, ratio=b['P_fe']/b['P_cu1'])
    r['a_D1'] = req(0.2, 1.0, T, e0, Rs0); r['a_D2'] = req(0.5, 1.0, T, e0, Rs0); r['a_D3'] = req(0.2, 1.05, T, e0, Rs0)
    r['a_D1_shop'] = req(0.2*p['kloop'], 1.0, T, e0, Rs0)
    out.append(r)
json.dump(out, open('mc2.json','w'))
print('N=%d in %.0f s' % (N, time.time()-t0))
from scipy.stats import spearmanr
for k in ('eta0','ratio','a_D1','a_D2','a_D3','a_D1_shop'):
    v = np.array([o[k] for o in out]); print('%-10s median %.4f 95%% [%.4f, %.4f]' % (k, np.median(v), np.percentile(v,2.5), np.percentile(v,97.5)))
for tgt in ('a_D1','a_D1_shop'):
    y = [o[tgt] for o in out]
    print(tgt, ' '.join('%s %.2f' % (k, spearmanr([o[k] for o in out], y)[0]) for k in ('cm','cr','cl','split','cfe','cfw','cll','T','cI','kloop')))
