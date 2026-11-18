import itertools, json, math, numpy as np
from hybrid2 import *
WIRES=[0.50,0.53,0.56,0.60,0.63,0.67,0.71,0.75,0.80,0.85,0.90,0.95,1.00,1.06,1.12]
Rs0=Rs_of(18,2,0.75); base=at_power(3700,Rs0); Ist0,Tst0,Tbd0=start_bd(Rs0)
eta0=base['eta']
rows=[]
for p in (7,8,9):
    for dmg,kfe,mlt in (('D0',1.0,1.0),('D1',1.2,1.0),('D2',1.5,1.0),('D3',1.2,1.05)):
        for Nl,ns,d in itertools.product(range(14,21),(1,2,3),WIRES):
            fl=fill_of(Nl,ns,d)
            if fl>0.47 or fl<0.22: continue
            kN=Nl/18
            Rs=Rs_of(Nl,ns,d,pitch=p,mlt_factor=mlt)
            r=at_power(3700,Rs,pitch=p,kN=kN,k_fe=kfe)
            Ist,Tst,Tbd=start_bd(Rs,pitch=p,kN=kN,k_fe=kfe)
            rows.append(dict(pitch=p,dmg=dmg,Nl=Nl,ns=ns,d=d,fill=fl,eta=r['eta'],deta=100*(r['eta']-eta0),I1=r['I1'],pf=r['pf'],n=r['n'],
                             Ist=Ist,Tst=Tst,Tbd=Tbd,ok=bool(Ist<=1.10*Ist0 and Tbd>=0.90*Tbd0 and Tst>=0.90*Tst0),
                             Pcu1=r['P_cu1'],Pcu2=r['P_cu2'],Pfe=r['P_fe'],PLL=r['P_LL'],Pfw=r['P_fw']))
json.dump(dict(base=dict(eta=eta0,Ist=Ist0,Tst=Tst0,Tbd=Tbd0,**{k:base[k] for k in ('I1','pf','n','P_cu1','P_cu2','P_fe','P_fw','P_LL')}),rows=rows),open('scan2.json','w'))
print('base eta %.2f%%, Ist %.1f, Tst %.1f, Tbd %.1f'%(100*eta0,Ist0,Tst0,Tbd0))
for dmg in ('D0','D1','D2','D3'):
    print('\n',dmg)
    for lim in (0.40,0.45):
        for p in (7,8,9):
            c=[r for r in rows if r['pitch']==p and r['dmg']==dmg and r['ok'] and r['fill']<=lim+1e-9]
            b=max(c,key=lambda r:r['eta'])
            print('  fill<=%.0f%% pitch %d/9: best N=%d %dx%.2f fill %.1f%%  deta %+.2f  I %.2f pf %.3f Ist %.1f Tst %.1f Tbd %.1f'%(100*lim,p,b['Nl'],b['ns'],b['d'],100*b['fill'],b['deta'],b['I1'],b['pf'],b['Ist'],b['Tst'],b['Tbd']))
    o=[r for r in rows if r['pitch']==7 and r['dmg']==dmg and r['Nl']==18 and r['ns']==2 and r['d']==0.75][0]
    print('  original winding: deta %+.2f'%o['deta'])
