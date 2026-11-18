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
