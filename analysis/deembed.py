import csv, math, json, numpy as np
def deembed(f, R_app, X_app, Lm, Lss):
    w=2*math.pi*f; Zp = R_app + 1j*(X_app - w*Lss); Zm = 1j*w*Lm
    Zr = 1/(1/Zp - 1/Zm); return Zr.real, Zr.imag/w
def split(lr, Llow):
    """find Lsig_s such that Lsig_s = Lsig_r at 50 Hz (equal split at 50 Hz), by bisection"""
    f50,r50,x50=[t for t in lr if t[0]>49][0]
    g=lambda Ls: deembed(f50,r50,x50,Llow-Ls,Ls)[1]-Ls
    a,b=1e-4,x50/(2*math.pi*f50)
    for _ in range(100):
        m=0.5*(a+b)
        if g(a)*g(m)<=0: b=m
        else: a=m
    Ls=0.5*(a+b); Lm=Llow-Ls
    return Ls, Lm, [(f,)+deembed(f,r,x,Lm,Ls) for f,r,x in lr]
R5=list(csv.DictReader(open('../results/sonuclar_adim5.csv'))); R6=list(csv.DictReader(open('../results/sonuclar_adim6.csv')))
sets={'meas':( [(float(r['f_Hz']),float(r['Rr_olc']),float(r['X_olc'])) for r in R5 if r['test']=='kilitli'], 0.2545)}
nl7=[(float(r['I_rms']),float(r['L_hes_H'])) for r in R5 if r['test']=='bosta']
sets['7']=([(float(r['f_Hz']),float(r['Rr_hes']),float(r['X_hes'])) for r in R5 if r['test']=='kilitli'], max(L for I,L in nl7 if I<2.2))
NL={'7':sorted(nl7)}
for p in ('8','9'):
    nl=[(float(r['I_rms']),float(r['L_hes_H'])) for r in R6 if r['pitch_slots']==p and r['test']=='bosta']; NL[p]=sorted(nl)
    sets[p]=([(float(r['f_Hz']),float(r['Rr_hes']),float(r['X_hes'])) for r in R6 if r['pitch_slots']==p and r['test']=='kilitli'], max(L for I,L in nl if I<2.2))
out={}
for k,(lr,Llow) in sets.items():
    Ls,Lm,res=split(sorted(lr),Llow)
    out[k]=dict(Ls=Ls,Lm=Lm,lr=res,nl=NL.get(k))
    print('%-4s Lsig_s=%.2f mH Lm=%.4f | '%(k,Ls*1e3,Lm)+'  '.join('f%.0f: R %.3f Lr %.2f'%(f,R,l*1e3) for f,R,l in res))
json.dump(out,open('deembedded.json','w'))
