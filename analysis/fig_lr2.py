import json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':9,'axes.linewidth':0.6,'axes.spines.top':False,'axes.spines.right':False,
    'axes.grid':True,'grid.color':'#e3e2dd','grid.linewidth':0.5,'savefig.dpi':300,'legend.frameon':False})
BLUE,ORANGE,AQUA,INK,MUTED='#2a78d6','#eb6834','#1baf7a','#22221f','#7a7970'
D=json.load(open('deembedded2.json')); M=D['meas']; S=D['7']
fm=[x[0] for x in M['lr']]; Rm=[x[1] for x in M['lr']]; Lrm=[x[2]*1e3 for x in M['lr']]
fs=[x[0] for x in S['lr']]; Rs=[x[1] for x in S['lr']]; Lrs=[x[2]*1e3 for x in S['lr']]
L3s,L3r=2.91,3.15
fig,ax=plt.subplots(1,3,figsize=(7.2,2.5),gridspec_kw={'wspace':0.45})
a=ax[0]; a.plot(fs,Rs,'-',color=BLUE,lw=2,marker='o',ms=4,label='2D FEA (calibrated)')
a.plot(fm,Rm,'o',ms=5,mfc='white',mec=INK,mew=1,label='Measurement (de-embedded)')
a.plot([2],[1.012],'D',ms=5,mfc='white',mec=AQUA,mew=1.3,label=r'TU Graz $R_{r,20^\circ C}$ (calibration)')
a.set_xlim(0,53); a.set_xlabel('Rotor frequency (Hz)'); a.set_ylabel(r"$R_r'$ at 20 °C ($\Omega$)"); a.legend(fontsize=6.3,loc='upper left')
a.set_title('(a) Rotor resistance',fontsize=9,loc='left')
a=ax[1]
tot_m=[M['Ls']*1e3+l for l in Lrm]; tot_2d=[S['Ls']*1e3+l for l in Lrs]; tot_c=[t+L3s+L3r for t in tot_2d]
a.plot(fs,tot_2d,'--',color=ORANGE,lw=1.5,marker='s',ms=3.5,label='2D FEA only')
a.plot(fs,tot_c,'-',color=BLUE,lw=2,marker='o',ms=4,label=r'2D FEA + $L_{3D}$')
a.plot(fm,tot_m,'o',ms=5,mfc='white',mec=INK,mew=1,label='Measurement (de-embedded)')
a.set_xlabel('Rotor frequency (Hz)'); a.set_ylabel(r'$L_{\sigma s}+L_{\sigma r}$ (mH)'); a.legend(fontsize=6.3,loc='lower left'); a.set_ylim(10,22)
a.set_title('(b) Total leakage inductance',fontsize=9,loc='left')
a=ax[2]
a.plot(fm,[m-s for m,s in zip(Lrm,[dict((round(f),l) for f,l in zip(fs,Lrs))[round(f)] for f in fm])],'-',color=AQUA,lw=2,marker='o',ms=4,label=r'rotor side $L_{3D,r}$')
a.axhline(L3s,color=INK,lw=1,ls=':'); a.text(50,L3s-0.18,r'stator side $L_{3D,s}$ = %.2f mH'%L3s,fontsize=6.5,ha='right',va='top')
a.set_ylim(2,4); a.set_xlabel('Rotor frequency (Hz)'); a.set_ylabel(r'$L_{3D}$ (mH)'); a.legend(fontsize=6.5,loc='upper right')
a.set_title('(c) 3D end-effect inductance',fontsize=9,loc='left')
fig.savefig('../figures/fig_lockedrotor_validation.png',bbox_inches='tight'); fig.savefig('../figures/fig_lockedrotor_validation.pdf',bbox_inches='tight')
print('tot err %', [round(100*(c-m)/m,1) for c,m in zip([dict(zip([round(f) for f in fs],tot_c))[round(f)] for f in fm],tot_m)])
