import json, numpy as np, math
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':9,'axes.linewidth':0.6,'axes.spines.top':False,'axes.spines.right':False,
    'axes.grid':True,'grid.color':'#e3e2dd','grid.linewidth':0.5,'savefig.dpi':300,'legend.frameon':False})
INK,MUTED='#22221f','#7a7970'
R=[r for r in json.load(open('general.json')) if r['a_mod']<2.49 and not math.isnan(r['a_rule'])]
fig,ax=plt.subplots(1,2,figsize=(7.2,3.0),gridspec_kw={'wspace':0.32})
a=ax[0]; cm=plt.get_cmap('viridis'); nrm=Normalize(0,1.25)
lim=(1,2.3); a.fill_between(lim,[x*0.98 for x in lim],[x*1.02 for x in lim],color='#e3e2dd',lw=0,label='±2 %')
a.plot(lim,lim,color=MUTED,lw=0.8)
for fam,mk in (('core','o'),('copper','s')):
    rr=[r for r in R if r['family']==fam]
    sc=a.scatter([r['a_mod'] for r in rr],[r['a_rule'] for r in rr],c=[r['ratio'] for r in rr],cmap=cm,norm=nrm,marker=mk,s=22,edgecolor=INK,linewidth=0.4,
                 label='varied core loss' if fam=='core' else 'varied stator resistance')
a.axvline(1.5,color=INK,ls='--',lw=0.7); a.text(1.52,1.05,'beyond practical\nfill headroom',fontsize=6.5)
a.set_xlim(*lim); a.set_ylim(*lim); a.set_xlabel(r'Required $A_{Cu}/A_{Cu,0}$, full model'); a.set_ylabel(r'Required $A_{Cu}/A_{Cu,0}$, shop rule')
a.legend(fontsize=6.5,loc='upper left'); cb=fig.colorbar(sc,ax=a,pad=0.02); cb.set_label(r'$P_{fe,0}/P_{Cu1,0}$',fontsize=8)
a.set_title('(a) Rule vs. model, 13 loss splits',fontsize=9,loc='left')
b=ax[1]; ratio=np.linspace(0.08,1.3,200)
for amax,ls in ((1.05,':'),(1.10,'-.'),(1.20,'--'),(1.33,'-')):
    dk=(1-1/amax)/ratio
    b.plot(ratio,100*dk,ls,color=INK,lw=1.3,label=r'$A_{max}/A_0$ = %.2f'%amax)
b.axhline(20,color='#eb6834',lw=0.9); b.text(1.28,22,'EASA 20 % limit',fontsize=6.5,ha='right',color=INK)
b.axvline(0.19,color='#2a78d6',lw=0.9); b.text(0.21,150,'this motor\n(0.19)',fontsize=6.5,color=INK)
b.set_ylim(0,200); b.set_xlim(0.08,1.3); b.set_xlabel(r'Core-to-stator-copper loss ratio $P_{fe,0}/P_{Cu1,0}$'); b.set_ylabel('Max. compensable core-loss increase (%)')
b.legend(fontsize=6.5,loc='upper right',title='copper-area headroom',title_fontsize=6.5); b.set_title('(b) Damage a rewind can absorb (MLT unchanged)',fontsize=9,loc='left')
fig.savefig('../figures/fig_rule_generality.png',bbox_inches='tight'); fig.savefig('../figures/fig_rule_generality.pdf',bbox_inches='tight')
ok=[r for r in R if r['a_mod']<=1.5]; print('n=%d within a<=1.5, max |err| = %.2f%%, mean |err| = %.2f%%'%(len(ok),max(abs(r['err']) for r in ok),np.mean([abs(r['err']) for r in ok])))
