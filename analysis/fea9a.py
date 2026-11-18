"""FEA run 9A: energy-based stator/rotor slot-leakage split, saturation at full-voltage start,
and sensitivity of the main results to the leakage split used in the separation step."""
import csv, json, math, numpy as np
import hybrid2 as H
from deembed import deembed
rows = list(csv.DictReader(open('../results/sonuclar_adim9A.csv')))
out = {'energy': [], 'sat': {}}
for r in rows:
    v = {k: float(x) for k, x in r.items()}
    kap = v['W_up0_block'] / v['W_up0_sampled']                 # block integral vs. point sampling (same region)
    Ws = v['W_coils'] + v['W_sep'] + kap * v['W_open']; Wr = v['W_bars']
    out['energy'].append(dict(pitch=int(v['pitch_slots']), f=v['f_Hz'], I=v['I_rms'], kappa=kap, W_stator_slot=Ws, W_rotor_slot=Wr,
                              W_airgap=v['W_airgap'], fs_slot=Ws/(Ws+Wr), X=v['X'], R=v['R_app']))
for e in out['energy']:
    print('pitch %d f=%4.1f I=%6.2f  kappa %.3f  stator slot share %.3f  airgap/(slots) %.2f  X %.3f R %.3f' % (
        e['pitch'], e['f'], e['I'], e['kappa'], e['fs_slot'], e['W_airgap']/(e['W_stator_slot']+e['W_rotor_slot']), e['X'], e['R']))
# saturation at start: X and R at the full-voltage starting current vs. 4.622 A, 50 Hz
for p in (7, 8, 9):
    E = {round(e['I'], 2): e for e in out['energy'] if e['pitch'] == p and e['f'] == 50.0}
    lo, hi = E[4.62], E[max(E)]
    out['sat'][p] = dict(I_start=max(E), X_ratio=hi['X']/lo['X'], R_ratio=hi['R']/lo['R'], X15=E[15.0]['X']/lo['X'])
    print('pitch %d: at %.1f A, X/X(4.6 A) = %.3f, R/R(4.6 A) = %.3f (15 A: X ratio %.3f)' % (p, max(E), hi['X']/lo['X'], hi['R']/lo['R'], E[15.0]['X']/lo['X']))

# ---- sensitivity of results to the leakage split fs = L_sig_s/(L_sig_s + L_sig_r') at 50 Hz
R5 = list(csv.DictReader(open('../results/sonuclar_adim5.csv'))); R7 = list(csv.DictReader(open('../results/sonuclar_adim7.csv')))
meas = sorted((float(r['f_Hz']), float(r['Rr_olc']), float(r['X_olc'])) for r in R5 if r['test'] == 'kilitli')
D0 = json.loads(json.dumps(H.D))
def split_fs(lr, Llow, fs):
    f50, r50, x50 = [t for t in lr if t[0] > 49][0]
    g = lambda Ls: deembed(f50, r50, x50, Llow-Ls, Ls)[1]*fs/(1-fs) - Ls
    a, b = 1e-4, x50/(2*math.pi*f50)
    for _ in range(100):
        m = 0.5*(a+b)
        if g(a)*g(m) <= 0: b = m
        else: a = m
    Ls = 0.5*(a+b); Lm = Llow-Ls
    return Ls, Lm, [(f,)+deembed(f, r, x, Lm, Ls) for f, r, x in lr]
def run_fs(fs):
    Lsm, Lmm, mres = split_fs(meas, 0.2545, fs)
    new = {}; sc = None
    for p in '789':
        lr = sorted((float(r['f_Hz']), float(r['R_app_hes']), float(r['X_hes'])) for r in R7 if r['pitch_slots'] == p)
        Llow = max(L for I, L in D0[p]['nl'] if I < 2.2)
        Ls, Lm, res = split_fs(lr, Llow, fs)
        if p == '7':
            sc = 1.012/np.interp(2.0, [x[0] for x in res], [x[1] for x in res])         # re-calibrate cage conductivity
            L3s = Lsm - Ls
            fm = [x[0] for x in mres]; L3r = float(np.mean([x[2] - np.interp(f, [y[0] for y in res], [y[2] for y in res]) for f, _, x2 in [(m[0], 0, 0) for m in mres] for x in [mres[fm.index(f)]]]))
        new[p] = dict(Ls=Ls, Lm=Lm, lr=[(f, R*sc, L) for f, R, L in res], nl=D0[p]['nl'])
    H.D.update(new); H.L3D_S, H.L3D_R = L3s, L3r; H.PITCH = {p: H.Pitch(p) for p in (7, 8, 9)}
    Rs0 = H.Rs_of(18, 2, 0.75); b = H.at_power(3700, Rs0)
    def areq(kfe, mlt):
        lo, hi = 0.9, 1.6
        for _ in range(35):
            a = 0.5*(lo+hi)
            if H.at_power(3700, Rs0*mlt/a, k_fe=kfe)['eta'] < b['eta']: lo = a
            else: hi = a
        return 0.5*(lo+hi)
    res = dict(fs=fs, Lss_meas=Lsm, L3D_S=L3s, L3D_R=L3r, sigma_scale=sc, n=b['n'], eta=b['eta'], I1=b['I1'], r=b['P_fe']/b['P_cu1'],
               D1=areq(1.2, 1.0), D2=areq(1.5, 1.0), D3=areq(1.2, 1.05))
    Ist0 = H.start_bd(Rs0)[0]; res['Ist'] = Ist0
    return res
out['split'] = []
for fs in (0.30, 0.35, 0.40, 0.45, 0.50):
    r = run_fs(fs); out['split'].append(r)
    print('fs %.2f: Lss(meas) %.2f mH, L3D s/r %.2f/%.2f mH, sigma x%.3f | n %.0f rpm, eta %.2f%%, I %.2f A, r %.3f | A_req D1 %.3f D2 %.3f D3 %.3f | Ist %.1f A' % (
        fs, 1e3*r['Lss_meas'], 1e3*r['L3D_S'], 1e3*r['L3D_R'], r['sigma_scale'], r['n'], 100*r['eta'], r['I1'], r['r'], r['D1'], r['D2'], r['D3'], r['Ist']))
json.dump(out, open('fea9a.json', 'w'), indent=1, default=float)

# ---- saturated starting current of the candidate designs (2D leakage at 50 Hz scaled by the FEA ratio at start)
import copy
H.D.clear(); H.D.update(copy.deepcopy(D0)); H.L3D_S, H.L3D_R = 2.91e-3, 3.15e-3
def make(sat):
    P = {}
    for p in (7, 8, 9):
        q = H.Pitch(p)
        if sat:
            k = out['sat'][p]['X_ratio']; d = H.D[str(p)]
            q.Lss = d['Ls']*k + H.L3D_S
            q.Lsr = np.array([x[2]*(k if x[0] > 49 else 1.0) for x in sorted(d['lr'])]) + H.L3D_R
        P[p] = q
    return P
cands = [(7,18,3,0.67),(7,18,2,0.85),(8,17,2,0.85),(8,17,3,0.71),(9,16,2,0.85),(9,17,3,0.71)]
out['start_sat'] = []
for sat in (False, True):
    H.PITCH = make(sat); Rs0 = H.Rs_of(18,2,0.75); I0, T0, B0 = H.start_bd(Rs0)
    for p, N, ns, d in cands:
        Rs = H.Rs_of(N, ns, d, pitch=p); I, T, B = H.start_bd(Rs, pitch=p, kN=N/18)
        out['start_sat'].append(dict(sat=sat, pitch=p, N=N, ns=ns, d=d, Ist=I/I0, Tst=T/T0, Tbd=B/B0, Ist_abs=I))
        print('%-5s %d/9 N%d %dx%.2f: Ist %.3f (%.1f A) Tst %.3f Tbd %.3f' % ('sat' if sat else 'unsat', p, N, ns, d, I/I0, I, T/T0, B/B0))
H.D.clear(); H.D.update(copy.deepcopy(D0)); H.PITCH = make(False)
json.dump(out, open('fea9a.json', 'w'), indent=1, default=float)
