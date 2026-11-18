import numpy as np, json, math
from hybrid2 import *
T_OP=75.0
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':9,'axes.linewidth':0.6,'axes.spines.top':False,'axes.spines.right':False,
    'savefig.dpi':300,'legend.frameon':False})
INK,MUTED,BLUE,ORANGE='#22221f','#7a7970','#2a78d6','#eb6834'
Rs0 = Rs_of(18,2,0.75)
eta0 = at_power(3700, Rs0)['eta']
kfe = np.linspace(1.0, 2.0, 21); ar = np.linspace(0.9, 1.5, 25)
res = {}
for mlt in (1.0, 1.05):
    D = np.zeros((len(ar), len(kfe)))
    for i,a in enumerate(ar):
        for j,k in enumerate(kfe):
            D[i,j] = 100*(at_power(3700, Rs0*mlt/a, k_fe=k)['eta'] - eta0)
    res[mlt] = D
np.save('compmap.npy', np.stack([res[1.0], res[1.05]]))
# required area ratio for full recovery (dEta = 0)
req = {m: [float(np.interp(0, res[m][:,j], ar)) for j in range(len(kfe))] for m in res}
print('eta0 %.2f%%' % (100*eta0))
for j in (0,4,10,20): print('core +%3.0f%%: req. Cu-area ratio MLT1.00 %.3f | MLT1.05 %.3f' % (100*(kfe[j]-1), req[1.0][j], req[1.05][j]))
json.dump({'kfe':kfe.tolist(),'req':{str(k):v for k,v in req.items()},'eta0':eta0}, open('compmap.json','w'))

designs = [(2,0.75,'2×0.75 (orig.)'),(2,0.80,'2×0.80'),(2,0.85,'2×0.85'),(3,0.71,'3×0.71'),(3,0.75,'3×0.75')]
fig, axs = plt.subplots(1,2,figsize=(7.2,3.1),sharey=True,gridspec_kw={'wspace':0.08})
lev = np.arange(-3,3.01,0.5)
for ax,(mlt,title) in zip(axs,[(1.0,'(a) Original end-winding length'),(1.05,'(b) MLT +5 % (longer end turns)')]):
    X = 100*(kfe-1)
    cf = ax.contourf(X, ar, res[mlt], levels=lev, cmap='RdBu', vmin=-3, vmax=3, extend='both')
    cs = ax.contour(X, ar, res[mlt], levels=[l for l in lev if abs(l)>1e-9], colors=INK, linewidths=0.3, alpha=0.5)
    ax.clabel(cs, fmt='%+.1f', fontsize=6.5, inline_spacing=2)
    c0 = ax.contour(X, ar, res[mlt], levels=[0], colors=INK, linewidths=1.8)
    ax.axvline(20, color=INK, ls='--', lw=0.8); ax.text(21, 0.93, 'EASA 20 %\nlimit', fontsize=6.5, va='bottom', color=INK)
    for yl,lab in ((0.40/0.30,'fill 40 %'),(0.45/0.30,'fill 45 %')):
        ax.axhline(yl, color=MUTED, lw=0.7, ls=':'); ax.text(97, yl-0.008, lab, fontsize=6.5, ha='right', va='top', color=INK)
    for n,d,lab in designs:
        y = fill_of(18,n,d)/0.30
        ax.plot([2], [y], marker='o', ms=4, mfc='white', mec=INK, mew=0.8, clip_on=False)
    ax.set_title(title, fontsize=9, loc='left', pad=6)
    ax.clabel(c0, fmt={0:r'$\Delta\eta=0$'}, fontsize=7, manual=[(70, float(np.interp(70, X, [float(np.interp(0, res[mlt][:,j], ar)) for j in range(len(kfe))])))]); ax.set_xlabel('Core-loss increase after repair (%)')
    ax.set_xlim(0,100); ax.set_ylim(0.9,1.5)
axs[0].set_ylabel('Copper area ratio  $A_{Cu}/A_{Cu,0}$  (N = 18)')
for n,d,lab in designs:
    y = fill_of(18,n,d)/0.30
    axs[0].text(4, y, lab, fontsize=6.3, va='center', color=INK, bbox=dict(fc='white', ec='none', alpha=0.7, pad=0.6))
cb = fig.colorbar(cf, ax=axs, shrink=0.9, pad=0.02); cb.set_label(r'$\Delta\eta$ vs. undamaged original (pp)')
fig.savefig('../figures/fig_compensation_map.png', bbox_inches='tight'); fig.savefig('../figures/fig_compensation_map.pdf', bbox_inches='tight')
