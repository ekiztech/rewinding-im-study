import json, math, numpy as np, importlib
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':9,'axes.linewidth':0.6,'axes.spines.top':False,'axes.spines.right':False,
    'axes.grid':True,'grid.color':'#e3e2dd','grid.linewidth':0.5,'savefig.dpi':300,'legend.frameon':False})
BLUE,ORANGE,INK,MUTED='#2a78d6','#eb6834','#22221f','#7a7970'
D=json.load(open('deembedded2.json')); M=D['meas']
# analytic: apparent vs true with measured circuit (Lm, Lss, Lsr(f), Rr(f) interpolated/extrapolated)
fm=np.array([x[0] for x in M['lr']]); Rm=np.array([x[1] for x in M['lr']]); Lrm=np.array([x[2] for x in M['lr']])
c=np.polyfit(fm,Rm,2); cl=np.polyfit(fm,Lrm,1)
f=np.linspace(1,50,200); eR=[]; eL=[]
for ff in f:
    w=2*math.pi*ff; Rr=np.polyval(c,ff); Lr=np.polyval(cl,ff); Zr=Rr+1j*w*Lr; Zm=1j*w*M['Lm']
    Zp=Zr*Zm/(Zr+Zm); eR.append(100*(Zp.real-Rr)/Rr); eL.append(100*((Zp.imag/w)-Lr)/Lr)
# points from data
pts=[(x[0],100*(r_app-x[1])/x[1]) for x,r_app in zip(M['lr'],[1.047,1.176,1.363,1.569,1.768])]
fig,ax=plt.subplots(1,2,figsize=(7.2,2.7),gridspec_kw={'wspace':0.35})
a=ax[0]; a.plot(f,eR,color=BLUE,lw=2,label=r"$R_r'$ read as $P/3I^2-R_s$"); a.plot(f,eL,color=ORANGE,lw=2,ls='--',label=r"$L_{\sigma r}$ read as $X/\omega - L_{\sigma s}$")
a.plot([p[0] for p in pts],[p[1] for p in pts],'o',mfc='white',mec=INK,ms=5,label='TU Graz measurements')
a.axhline(0,color=MUTED,lw=0.6); a.set_xlabel('Locked-rotor test frequency (Hz)'); a.set_ylabel('Error of apparent value (%)')
a.set_ylim(-40,60); a.legend(fontsize=6.5,loc='upper right'); a.set_title('(a) Apparent vs. true rotor parameters',fontsize=9,loc='left')
# impact on rated-point prediction: calibrate Rr from apparent 2 Hz value -> true Rr overestimated by factor k
import hybrid2 as H
Rs=H.Rs_of(18,2,0.75); base=H.at_power(3700,Rs)
P=H.PITCH[7]; R0=P.Rr.copy()
ks=np.linspace(1.0,1.25,11); ds=[]; de=[]; dn=[]
for k in ks:
    P.Rr=R0*k; r=H.at_power(3700,Rs); ds.append(100*(r['s']-base['s'])/base['s']); de.append(100*(r['eta']-base['eta'])); dn.append(r['n']-base['n'])
P.Rr=R0
b=ax[1]; b.plot(100*(ks-1),ds,color=BLUE,lw=2,label='rated slip error (%)'); b.plot(100*(ks-1),[10*x for x in de],color=ORANGE,lw=2,ls='--',label=r'efficiency error ($\times$10, pp)')
b.axvline(20,color=INK,ls=':',lw=0.8); b.text(19.5,5,'2 Hz case\n(+20 %)',fontsize=6.5,ha='right')
b.set_xlabel(r"Overestimate of $R_r'$ from skipping de-embedding (%)"); b.set_ylabel('Prediction error'); b.legend(fontsize=6.5,loc='upper left')
b.set_title('(b) Effect on rated-point prediction',fontsize=9,loc='left')
fig.savefig('../figures/fig_deembedding.png',bbox_inches='tight'); fig.savefig('../figures/fig_deembedding.pdf',bbox_inches='tight')
i=np.argmin(abs(ks-1.2)); print('at +20%%: slip %+.1f%%, speed %+.1f rpm, eta %+.2f pp'%(ds[i],dn[i],de[i]))
print('analytic R error at 2,10,50 Hz:',[round(float(np.interp(x,f,eR)),1) for x in (2,10,50)], ' L error:',[round(float(np.interp(x,f,eL)),1) for x in (2,10,50)])
