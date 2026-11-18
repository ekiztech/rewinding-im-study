import json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':9,'axes.linewidth':0.6,'axes.spines.top':False,'axes.spines.right':False,
    'axes.grid':True,'grid.color':'#e3e2dd','grid.linewidth':0.5,'savefig.dpi':300,'legend.frameon':False})
BLUE,ORANGE,INK='#2a78d6','#eb6834','#22221f'
M=json.load(open('mc.json'))
def rank(x): return np.argsort(np.argsort(x))
def spear(a,b): return np.corrcoef(rank(np.array(a)),rank(np.array(b)))[0,1]
par=[('cm','Magnetizing curve (±3 %)'),('cr','Rotor resistance (σ 5 %)'),('cl','2D leakage (±10 %)'),('split','L3D stator/rotor split'),
     ('cfe','Core-loss data (σ 10 %)'),('cfw','Friction & windage (σ 15 %)'),('cll','Additional load loss (×0.5–1.5)'),('T','Temperature 60–95 °C')]
fig,ax=plt.subplots(1,2,figsize=(7.2,2.9),gridspec_kw={'wspace':0.9,'width_ratios':[1,1.15]})
a=ax[0]; dks=[0,20,50,100]
for mlt,col,lab in ((1.0,BLUE,'MLT unchanged'),(1.05,ORANGE,'MLT +5 %')):
    med=[mlt]; lo=[mlt]; hi=[mlt]
    for dk in (0.2,0.5,1.0):
        v=np.array([m['a_%g_%g'%(dk,mlt)] for m in M]); med.append(np.median(v)); lo.append(np.percentile(v,2.5)); hi.append(np.percentile(v,97.5))
    a.fill_between(dks,lo,hi,color=col,alpha=0.22,lw=0); a.plot(dks,med,'-o',color=col,lw=2,ms=4,label=lab+' (median, 95 %)')
a.axvline(20,color=INK,ls='--',lw=0.7); a.text(22,0.995,'EASA limit',fontsize=6.5)
a.set_xlabel('Core-loss increase (%)'); a.set_ylim(0.98,1.42); a.set_ylabel(r'Required $A_{Cu}/A_{Cu,0}$'); a.legend(fontsize=6.3,loc='upper left',bbox_to_anchor=(0,1.0))
a.set_title('(a) Compensation need, N = %d'%len(M),fontsize=9,loc='left')
b=ax[1]; tgt=[m['a_0.5_1'] for m in M]; eta=[m['eta0'] for m in M]
rs=[(lab,spear([m[k] for m in M],tgt),spear([m[k] for m in M],eta)) for k,lab in par]
rs.sort(key=lambda t:abs(t[1]))
y=np.arange(len(rs)); h=0.38
b.barh(y+h/2,[r[1] for r in rs],h,color=BLUE,label=r'required $A_{Cu}$ ratio (+50 % core loss)')
b.barh(y-h/2,[r[2] for r in rs],h,color='#b8b7ae',label=r'absolute efficiency $\eta_0$')
b.set_yticks(y); b.set_yticklabels([r[0] for r in rs],fontsize=7); b.axvline(0,color=INK,lw=0.6); b.set_xlim(-1,1)
b.set_xlabel('Spearman rank correlation'); b.legend(fontsize=6.3,loc='upper center',bbox_to_anchor=(0.3,-0.2),ncol=1)
b.set_title('(b) What drives the uncertainty',fontsize=9,loc='left')
fig.savefig('../figures/fig_uncertainty.png',bbox_inches='tight'); fig.savefig('../figures/fig_uncertainty.pdf',bbox_inches='tight')
for r in sorted(rs,key=lambda t:-abs(t[1])): print('%-32s a: %+.2f  eta0: %+.2f'%r)
