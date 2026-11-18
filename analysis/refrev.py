"""Reference-motor analyses requested in review: insulated-wire fill, start-limit and MLT sensitivity,
thermal feedback, part load, nameplate-current consistency, design counts."""
import itertools, json, math, numpy as np
import hybrid2 as H
from rule2 import rule_v1, rule_v2
WIRES = [0.50,0.53,0.56,0.60,0.63,0.67,0.71,0.75,0.80,0.85,0.90,0.95,1.00,1.06,1.12]
DINS = lambda d: float(np.interp(d, [0.50,0.75,1.00,1.12], [0.566,0.826,1.093,1.217]))   # grade-2 max overall dia. (approx.)
A_COIL, PERIM, T_LINER = 99.71, 55.6, 0.25            # mm2, mm, mm (from slot geometry; liner thickness assumed)
A_NET = A_COIL - PERIM*T_LINER
def kins(Nl, ns, d): return 2*Nl*ns*DINS(d)**2/A_NET     # insulated lay factor, d^2 basis
OUT = {}
Rs0 = H.Rs_of(18,2,0.75); base = H.at_power(3700, Rs0); eta0 = base['eta']; Ist0, Tst0, Tbd0 = H.start_bd(Rs0)
OUT['net_area'] = dict(A_coil=A_COIL, A_net=A_NET, k_orig=kins(18,2,0.75), cu_fill_net_orig=72*math.pi*0.75**2/4/A_NET)
print('net area %.1f mm2; original winding k_ins=%.3f, bare-Cu fill on net area %.3f' % (A_NET, kins(18,2,0.75), OUT['net_area']['cu_fill_net_orig']))

# ---------------- 1. scan with insulated fill, start-limit sensitivity, MLT-per-pitch sensitivity
SCEN = (('D0',1.0,1.0),('D1',1.2,1.0),('D2',1.5,1.0),('D3',1.2,1.05))
def scan(mlt_inc=0.0297):
    H.MLT[8] = 0.553 + mlt_inc; H.MLT[9] = 0.553 + 2*mlt_inc
    rows = []
    for p in (7,8,9):
        for dmg,kfe,mlt in SCEN:
            for Nl,ns,d in itertools.product(range(14,21),(1,2,3),WIRES):
                k = kins(Nl,ns,d)
                if k > 0.80 or k < 0.40: continue
                Rs = H.Rs_of(Nl,ns,d,pitch=p,mlt_factor=mlt); kN = Nl/18
                r = H.at_power(3700,Rs,pitch=p,kN=kN,k_fe=kfe); Ist,Tst,Tbd = H.start_bd(Rs,pitch=p,kN=kN,k_fe=kfe)
                rows.append(dict(pitch=p,dmg=dmg,Nl=Nl,ns=ns,d=d,k=k,deta=100*(r['eta']-eta0),Ist=Ist/Ist0,Tst=Tst/Tst0,Tbd=Tbd/Tbd0))
    H.MLT[8] = 0.553 + 0.0297; H.MLT[9] = 0.553 + 2*0.0297
    return rows
def best(rows, p, dmg, klim, ilim=1.10):
    c = [r for r in rows if r['pitch']==p and r['dmg']==dmg and r['k']<=klim+1e-9 and r['Ist']<=ilim and r['Tst']>=0.9 and r['Tbd']>=0.9]
    return max(c, key=lambda r: r['deta']) if c else None
rows = scan()
nwind = len({(r['pitch'],r['Nl'],r['ns'],r['d']) for r in rows})
OUT['counts'] = dict(windings=nwind, evaluations=len(rows), admissible_k070=sum(1 for r in rows if r['k']<=0.70 and r['Ist']<=1.1 and r['Tst']>=0.9 and r['Tbd']>=0.9),
                     admissible_k075=sum(1 for r in rows if r['k']<=0.75 and r['Ist']<=1.1 and r['Tst']>=0.9 and r['Tbd']>=0.9))
print('designs: %d windings, %d evaluations' % (nwind, len(rows)), OUT['counts'])
T4 = {}
for klim in (0.65, 0.70, 0.75):
    for dmg in ('D0','D1','D2','D3'):
        for p in (7,8,9):
            for ilim in (1.05, 1.10, 1.15):
                b = best(rows, p, dmg, klim, ilim)
                T4['%s|%d|%.2f|%.2f' % (dmg,p,klim,ilim)] = b
OUT['best'] = T4
orig = {dmg: [r for r in rows if r['pitch']==7 and r['dmg']==dmg and r['Nl']==18 and r['ns']==2 and r['d']==0.75][0]['deta'] for dmg in ('D0','D1','D2','D3')}
OUT['orig'] = orig
print('original winding deta:', {k: round(v,2) for k,v in orig.items()})
for klim in (0.65,0.70,0.75):
    for ilim in (1.05,1.10,1.15):
        s = ' '.join('%d/9:%s' % (p, ('N%d %dx%.2f k%.2f %+.2f' % (b['Nl'],b['ns'],b['d'],b['k'],b['deta'])) if (b:=T4['D1|%d|%.2f|%.2f'%(p,klim,ilim)]) else '-') for p in (7,8,9))
        print('D1 k<=%.2f Ist<=%.2f: %s' % (klim, ilim, s))
# MLT increment sensitivity for the pitch comparison (D1)
OUT['mlt_sens'] = {}
for inc_pct in (5.4, 7.2, 9.0, 10.8):
    rr = scan(0.553*inc_pct/100)
    OUT['mlt_sens'][str(inc_pct)] = {('%d|%.2f' % (p,kl)): (best(rr,p,'D1',kl) or {}).get('deta') for p in (7,8,9) for kl in (0.70,0.75)}
    print('MLT +%.1f%%/pitch step:' % inc_pct, {k: (round(v,2) if v else None) for k,v in OUT['mlt_sens'][str(inc_pct)].items()})

# ---------------- 2. required area: fixed 75 C vs thermal feedback
def thermal(Rs_fun, kfe, dT0, w=0.6, P2=3700, ref=None):
    T = 25 + dT0
    for _ in range(12):
        r = H.at_power(P2, Rs_fun(T), k_fe=kfe, Tr=T)
        Ptot = r['P1'] - r['P2']; q = w*r['P_cu1'] + (1-w)*(Ptot - r['P_cu1'])
        if ref is None: return r, T, q
        T = 25 + dT0*q/ref
    return r, T, q
OUT['thermal'] = {}
for dT0 in (50, 70, 90):
    rb, Tb, qref = thermal(lambda T: H.Rs_of(18,2,0.75,T=T), 1.0, dT0)
    e0 = rb['eta']; res = {}
    for dmg,kfe,mlt in SCEN[1:]:
        rd, Td, _ = thermal(lambda T: H.Rs_of(18,2,0.75,T=T)*mlt, kfe, dT0, ref=qref)
        lo, hi = 0.9, 1.6
        for _ in range(35):
            a = 0.5*(lo+hi); r, Ta, _ = thermal(lambda T: H.Rs_of(18,2,0.75,T=T)*mlt/a, kfe, dT0, ref=qref)
            if r['eta'] < e0: lo = a
            else: hi = a
        res[dmg] = dict(a=0.5*(lo+hi), T_damaged=Td, T_comp=Ta, deta_uncomp=100*(rd['eta']-e0))
    OUT['thermal'][str(dT0)] = dict(eta0=e0, T0=Tb, **res)
    print('thermal dT0=%d (base %.0f C, eta0 %.2f%%): ' % (dT0, Tb, 100*e0) + ' '.join('%s a=%.3f (T dmg %.0f, comp %.0f, deta %.2f)' % (k, v['a'], v['T_damaged'], v['T_comp'], v['deta_uncomp']) for k,v in res.items()))

# ---------------- 3. part load (fixed temperature, as in main study)
def a_req(P2, kfe, mlt):
    e0 = H.at_power(P2, Rs0)['eta']; lo, hi = 0.9, 2.0
    for _ in range(40):
        a = 0.5*(lo+hi)
        if H.at_power(P2, Rs0*mlt/a, k_fe=kfe)['eta'] < e0: lo = a
        else: hi = a
    return 0.5*(lo+hi)
OUT['partload'] = {}
for load in (0.5, 0.75, 1.0):
    b = H.at_power(3700*load, Rs0); r = b['P_fe']/b['P_cu1']; r2 = b['P_cu2']/b['P_cu1']; dl = b['I1']*Rs0*b['pf']/H.V_PH
    d = {}
    for dmg,kfe,mlt in SCEN[1:]:
        d[dmg] = dict(mod=a_req(3700*load,kfe,mlt), v1=rule_v1(kfe-1,r,mlt), v2=rule_v2(kfe-1,r,r2,dl,mlt))
    OUT['partload'][str(load)] = dict(r=r, r2=r2, delta=dl, eta0=b['eta'], **d)
    print('load %.2f: eta0 %.2f%% r %.3f r2 %.3f delta %.3f | ' % (load, 100*b['eta'], r, r2, dl) + ' '.join('%s mod %.3f v1 %.3f v2 %.3f' % (k, v['mod'], v['v1'], v['v2']) for k,v in d.items()))
# efficiency vs load curves for Fig.
OUT['curves'] = {}
loads = np.linspace(0.25, 1.25, 21)
for lbl, Rs, kfe in (('orig_D0', Rs0, 1.0), ('orig_D1', Rs0, 1.2), ('2x0.85_D1', H.Rs_of(18,2,0.85), 1.2), ('orig_D2', Rs0, 1.5), ('2x0.85_D2', H.Rs_of(18,2,0.85), 1.5)):
    OUT['curves'][lbl] = [100*H.at_power(3700*L, Rs, k_fe=kfe)['eta'] for L in loads]
OUT['loads'] = loads.tolist()

# ---------------- 4. nameplate-current consistency and loss-ratio sensitivity
I1, pf = base['I1'], base['pf']; Iact, Ireac = I1*pf, I1*math.sin(math.acos(pf))
Ineed_act = math.sqrt(max(6.9**2 - Ireac**2, 0)); eta_needed = 3700/(3*H.V_PH*Ineed_act)  # input power at 6.9 A with the model's reactive current
OUT['nameplate'] = dict(I1=I1, pf=pf, I_active=Iact, I_reactive=Ireac, eta_needed_for_6p9A=eta_needed)
print('model I1 %.2f A (active %.2f, reactive %.2f); 6.9 A would need eta >= %.1f%%' % (I1, Iact, Ireac, 100*eta_needed))
OUT['ratio_sens'] = {}
for lbl, Iref in (('7.70 A (model)', 7.70), ('7.2 A (data-provider parameters)', 7.2), ('6.9 A (nameplate)', 6.9)):
    sc = (Iref/7.70)**2; d = {}
    for dmg,kfe,mlt in SCEN[1:]:
        e0 = H.at_power(3700, Rs0*sc)['eta']; lo, hi = 0.9, 1.6
        for _ in range(35):
            a = 0.5*(lo+hi)
            if H.at_power(3700, Rs0*sc*mlt/a, k_fe=kfe)['eta'] < e0: lo = a
            else: hi = a
        d[dmg] = 0.5*(lo+hi)
    b = H.at_power(3700, Rs0*sc); OUT['ratio_sens'][lbl] = dict(r=b['P_fe']/b['P_cu1'], **d)
    print('P_Cu1 scaled to %-34s r=%.3f  ' % (lbl, b['P_fe']/b['P_cu1']) + ' '.join('%s %.3f' % (k,v) for k,v in d.items()))
json.dump(OUT, open('refrev.json','w'), indent=1, default=float)
