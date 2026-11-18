"""De-embed the magnetizing branch from locked-rotor data (measured and FEA) -> deembedded2.json"""
import csv, json
from deembed import split
R5=list(csv.DictReader(open('../results/sonuclar_adim5.csv'))); R6=list(csv.DictReader(open('../results/sonuclar_adim6.csv'))); R7=list(csv.DictReader(open('../results/sonuclar_adim7.csv')))
meas=[(float(r['f_Hz']),float(r['Rr_olc']),float(r['X_olc'])) for r in R5 if r['test']=='kilitli']
nl={'7':sorted((float(r['I_rms']),float(r['L_hes_H'])) for r in R5 if r['test']=='bosta')}
for p in '89': nl[p]=sorted((float(r['I_rms']),float(r['L_hes_H'])) for r in R6 if r['pitch_slots']==p and r['test']=='bosta')
Ls_m,Lm_m,mres=split(sorted(meas),0.2545)
out={'meas':dict(Ls=Ls_m,Lm=Lm_m,lr=mres)}
for p in '789':
    lr=sorted((float(r['f_Hz']),float(r['R_app_hes']),float(r['X_hes'])) for r in R7 if r['pitch_slots']==p)
    Ls,Lm,res=split(lr,max(L for I,L in nl[p] if I<2.2)); out[p]=dict(Ls=Ls,Lm=Lm,lr=res,nl=nl[p])
json.dump(out,open('deembedded2.json','w')); print('written deembedded2.json')
