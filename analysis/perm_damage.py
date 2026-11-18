"""Core damage with degraded permeability (effective lamination fill 0.84/0.80) + core-loss increase."""
import csv, json, numpy as np
import hybrid2 as H
P=H.PITCH[7]; curves={0.87:(P.I.copy(),P.psi.copy())}
R=list(csv.DictReader(open('../results/sonuclar_adim8_hasar.csv')))
for fill in (0.84,0.80):
    d=sorted((float(r['I_rms']),float(r['L_hes_H'])) for r in R if abs(float(r['lam_fill'])-fill)<1e-6)
    I=np.array([a for a,b in d]); L=np.array([b for a,b in d]); curves[fill]=(I,(L-P.Lss)*I)
Rs0=H.Rs_of(18,2,0.75); eta0=H.at_power(3700,Rs0)['eta']
def use(f): P.I,P.psi=curves[f]
def req(dk,f):
    use(f); lo,hi=0.9,1.8
    for _ in range(24):
        m=0.5*(lo+hi)
        if H.at_power(3700,Rs0/m,k_fe=1+dk)['eta']<eta0: lo=m
        else: hi=m
    return 0.5*(lo+hi)
out=[]
for f in (0.87,0.84,0.80):
    for dk in (0.0,0.2,0.5):
        use(f); r=H.at_power(3700,Rs0,k_fe=1+dk); a=req(dk,f) if (f!=0.87 or dk>0) else 1.0
        out.append(dict(fill=f,dk=dk,eta=r['eta'],deta=100*(r['eta']-eta0),I1=r['I1'],pf=r['pf'],a_req=a,P_cu1=r['P_cu1'],P_fe=r['P_fe']))
        print(out[-1])
use(0.87); json.dump(out,open('perm_damage.json','w'))
