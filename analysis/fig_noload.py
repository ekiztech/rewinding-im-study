import csv, math, json
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':9,'axes.linewidth':0.6,
    'axes.spines.top':False,'axes.spines.right':False,'xtick.major.width':0.6,'ytick.major.width':0.6,
    'axes.grid':True,'grid.color':'#e3e2dd','grid.linewidth':0.5,'legend.frameon':False,'savefig.dpi':300})
BLUE,ORANGE,AQUA,INK,MUTED='#2a78d6','#eb6834','#1baf7a','#22221f','#7a7970'
R=list(csv.DictReader(open('../results/sonuclar_adim5.csv')))
nl=[r for r in R if r['test']=='bosta']; lr=[r for r in R if r['test']=='kilitli']
I=[float(r['I_rms']) for r in nl]; Vm=[float(r['V_olc']) for r in nl]; Vc=[float(r['V_hes']) for r in nl]
Vn=[279.4,273.1,265.3,254.9,240.1,218.1,188.5,172.6,156.4,140.2,123.9]   # nominal (uncalibrated) model, step 2
o=sorted(range(len(I)),key=lambda i:I[i]); I=[I[i] for i in o]; Vm=[Vm[i] for i in o]; Vc=[Vc[i] for i in o]; Vn=[Vn[i] for i in o]

# ---------- Fig 1: no-load ----------
fig,(a,b)=plt.subplots(2,1,figsize=(3.5,4.2),sharex=True,gridspec_kw={'height_ratios':[2.2,1],'hspace':0.12})
a.plot(I,Vn,'--',color=ORANGE,lw=1.5,label='2D FEA, nominal parameters')
a.plot(I,Vc,'-',color=BLUE,lw=2,label='2D FEA, calibrated')
a.plot(I,Vm,'o',ms=5,mfc='white',mec=INK,mew=1,label='Measurement (TU Graz)')
a.set_ylabel('Phase voltage (V rms)'); a.legend(loc='lower right',fontsize=7.5)
b.axhline(0,color=MUTED,lw=0.6)
b.plot(I,[100*(v-m)/m for v,m in zip(Vn,Vm)],'--',marker='s',ms=3.5,color=ORANGE,lw=1.2)
b.plot(I,[100*(v-m)/m for v,m in zip(Vc,Vm)],'-',marker='o',ms=3.5,color=BLUE,lw=1.5)
b.set_ylabel('Error (%)'); b.set_xlabel('Stator current (A rms)'); b.set_ylim(-5,15)
a.set_title('No-load test, 50 Hz',fontsize=9,loc='left')
fig.savefig('../figures/fig_noload_validation.png',bbox_inches='tight'); fig.savefig('../figures/fig_noload_validation.pdf',bbox_inches='tight')

