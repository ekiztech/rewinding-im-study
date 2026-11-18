import json, math, numpy as np
from hybrid2 import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':9,'axes.linewidth':0.6,'axes.spines.top':False,'axes.spines.right':False,
    'axes.grid':True,'grid.color':'#e3e2dd','grid.linewidth':0.5,'savefig.dpi':300,'legend.frameon':False})
S=json.load(open('scan2.json')); B=S['base']; eta0=B['eta']
COL={7:'#2a78d6',8:'#eb6834',9:'#1baf7a'}; MK={7:'o',8:'s',9:'^'}
fills=np.arange(0.30,0.4601,0.01)
def cont_best(p,lim,kfe=1.2):
    best=None
    for Nl in range(14,21):
        kN=Nl/18; A=lim*A_SLOT/(2*Nl)
        Rs=RHO20*(1+ALPHA_CU*55)*(Nl*NPH0/N0)*MLT[p]/A
        r=at_power(3700,Rs,pitch=p,kN=kN,k_fe=kfe); Ist,Tst,Tbd=start_bd(Rs,pitch=p,kN=kN,k_fe=kfe)
        if Ist<=1.1*B['Ist'] and Tbd>=0.9*B['Tbd'] and Tst>=0.9*B['Tst']:
            d=100*(r['eta']-eta0)
            if best is None or d>best[0]: best=(d,Nl)
    return best
res={}
for p in (7,8,9):
    raw=[cont_best(p,f) for f in fills]; env=[]; cur=None
    for b in raw:
        if b is not None and (cur is None or b[0]>cur[0]): cur=b
        env.append(cur)
    res[p]=env
json.dump({str(p):[[float(f),b[0],b[1]] for f,b in zip(fills,res[p])] for p in res},open('pitch_cont.json','w'))
fig,ax=plt.subplots(figsize=(3.6,2.9))
for p in (7,8,9):
    y=np.array([b[0] for b in res[p]])
    pen={7:(0.2,1.4),8:(1.1,3.8),9:(6.0,22.6)}[p]          # belt-harmonic cage loss, W (estimate x0.5..x2)
    P1=4730.0
    ax.plot(100*fills,y,'-',color=COL[p],lw=2,label='pitch %d/9'%p)
    if p!=7: ax.fill_between(100*fills,y-100*pen[1]/P1,y-100*pen[0]/P1,color=COL[p],alpha=0.25,lw=0)
    pts=[r for r in S['rows'] if r['pitch']==p and r['dmg']=='D1' and r['ok']]
    xs,ys=[],[]
    for f in fills:
        c=[r for r in pts if r['fill']<=f+1e-9]
        if c: xs.append(100*f); ys.append(max(r['deta'] for r in c))
    ax.step(xs,ys,where='post',color=COL[p],lw=0.9,ls=':')
ax.axhline(0,color='#7a7970',lw=0.7)
ax.axvline(40,color='#22221f',lw=0.6,ls='--'); ax.text(39.7,-0.1,'routine hand-\nwinding limit',fontsize=6.3,va='top',ha='right')
ax.set_xlabel('Maximum copper fill of the rewind (%)'); ax.set_ylabel(r'$\Delta\eta$ vs. undamaged original (pp)')
ax.set_title('Core loss +20 % (EASA limit): best design per pitch',fontsize=8.5,loc='left')
ax.legend(fontsize=6.8,loc='upper left'); ax.text(45.8,0.15,'dotted: IEC 60317 wire steps; band: with\nbelt-harmonic cage loss',fontsize=6.3,ha='right',color='#22221f')
ax.set_xlim(30,46)
fig.savefig('../figures/fig_pitch_variants.png',bbox_inches='tight'); fig.savefig('../figures/fig_pitch_variants.pdf',bbox_inches='tight')
for p in (7,8,9): print(p,[('%.0f%%'%(100*f),round(b[0],2),b[1]) for f,b in list(zip(fills,res[p]))[::2]])
