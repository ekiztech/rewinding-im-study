"""Design scan with saturated starting parameters: the 2D leakage at 50 Hz is scaled per pitch by the FEA ratio
X(start current)/X(4.6 A) of run 9A before the starting current and torques are evaluated."""
import itertools, json, numpy as np, copy
import hybrid2 as H
import refrev_lib as L
SAT = {int(k): v['X_ratio'] for k, v in json.load(open('fea9a.json'))['sat'].items()}
UNSAT = {p: H.Pitch(p) for p in (7, 8, 9)}
def sat_pitch():
    P = {}
    for p in (7, 8, 9):
        q = H.Pitch(p); d = H.D[str(p)]; k = SAT[p]
        q.Lss = d['Ls']*k + H.L3D_S
        q.Lsr = np.array([x[2]*(k if x[0] > 49 else 1.0) for x in sorted(d['lr'])]) + H.L3D_R
        P[p] = q
    return P
SATP = sat_pitch()
def start_sat(Rs, **kw):
    H.PITCH = SATP; r = H.start_bd(Rs, **kw); H.PITCH = UNSAT; return r
Rs0 = H.Rs_of(18, 2, 0.75); base = H.at_power(3700, Rs0); eta0 = base['eta']; I0, T0, B0 = start_sat(Rs0)
def scan(mlt_inc=0.0297):
    H.MLT[8] = 0.553 + mlt_inc; H.MLT[9] = 0.553 + 2*mlt_inc; rows = []
    for p in (7, 8, 9):
        for dmg, kfe, mlt in L.SCEN:
            for Nl, ns, d in itertools.product(range(14, 21), (1, 2, 3), L.WIRES):
                k = L.kins(Nl, ns, d)
                if k > 0.80 or k < 0.40: continue
                Rs = H.Rs_of(Nl, ns, d, pitch=p, mlt_factor=mlt); kN = Nl/18
                r = H.at_power(3700, Rs, pitch=p, kN=kN, k_fe=kfe); I, T, B = start_sat(Rs, pitch=p, kN=kN, k_fe=kfe)
                rows.append(dict(pitch=p, dmg=dmg, Nl=Nl, ns=ns, d=d, k=k, deta=100*(r['eta']-eta0), Ist=I/I0, Tst=T/T0, Tbd=B/B0))
    H.MLT[8] = 0.553 + 0.0297; H.MLT[9] = 0.553 + 2*0.0297
    return rows
if __name__ == '__main__':
    rows = scan(); out = {'best': {}, 'counts': {}}
    for kl in (0.70, 0.75):
        out['counts'][str(kl)] = sum(1 for r in rows if r['k'] <= kl and r['Ist'] <= 1.1 and r['Tst'] >= 0.9 and r['Tbd'] >= 0.9)
        for il in (1.05, 1.10, 1.15):
            for dmg in ('D0', 'D1', 'D2', 'D3'):
                for p in (7, 8, 9):
                    out['best']['%s|%d|%.2f|%.2f' % (dmg, p, kl, il)] = L.best(rows, p, dmg, kl, il)
    for kl in (0.70, 0.75):
        for il in (1.05, 1.10, 1.15):
            print('D1 k<=%.2f Ist<=%.2f: ' % (kl, il) + ' '.join('%d/9:%s' % (p, ('N%d %dx%.2f %+.2f (Ist %.3f)' % (b['Nl'], b['ns'], b['d'], b['deta'], b['Ist'])) if (b := out['best']['D1|%d|%.2f|%.2f' % (p, kl, il)]) else '-') for p in (7, 8, 9)))
    print('counts', out['counts'])
    for dmg in ('D0', 'D1', 'D2', 'D3'):
        print(dmg, ' '.join('k%.2f: %s' % (kl, (lambda b: 'N%d %dx%.2f %+.2f Ist %.3f Tbd %.3f' % (b['Nl'], b['ns'], b['d'], b['deta'], b['Ist'], b['Tbd']))(out['best']['%s|7|%.2f|1.10' % (dmg, kl)])) for kl in (0.70, 0.75)))
    json.dump(out, open('scan_sat.json', 'w'), indent=1, default=float)
