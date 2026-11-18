import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import scan_sat as S, refrev_lib as L
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':9,'axes.linewidth':0.6,'axes.spines.top':False,'axes.spines.right':False,
    'axes.grid':True,'grid.color':'#e3e2dd','grid.linewidth':0.5,'savefig.dpi':300,'legend.frameon':False})
MUTED = '#7a7970'; COL = {7:'#2a78d6', 8:'#eb6834', 9:'#1baf7a'}
ks = np.arange(0.55, 0.801, 0.01)
rows = {inc: S.scan(0.553*inc) for inc in (0.054, 0.108)}
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=True, gridspec_kw={'wspace':0.08})
for ax, il, title in ((axs[0], 1.10, '(a) Starting current +10 % (saturated)'), (axs[1], 1.15, '(b) Starting current +15 % (saturated)')):
    for p in (7, 8, 9):
        env = {}
        for inc in (0.054, 0.108):
            y = [(b['deta'] if (b := L.best(rows[inc], p, 'D1', k, il)) else np.nan) for k in ks]
            env[inc] = np.fmax.accumulate(np.nan_to_num(y, nan=-9))
        if p != 7: ax.fill_between(ks, env[0.108], env[0.054], step='post', color=COL[p], alpha=0.18, lw=0)
        ax.step(ks, env[0.054], where='post', color=COL[p], lw=1.4, label='%d/9' % p)
    ax.axvline(L.kins(18, 2, 0.75), color=MUTED, lw=0.7, ls=':')
    ax.set_xlim(0.55, 0.80); ax.set_ylim(0, 3); ax.set_xlabel('Allowed insulated lay factor $k_{ins}$'); ax.set_title(title, fontsize=9, loc='left')
axs[0].set_xticks([0.55, 0.60, 0.65, 0.70, 0.75]); axs[0].set_ylabel(r'Best $\Delta\eta$ vs. undamaged original (pp)')
axs[0].text(L.kins(18,2,0.75)-0.004, 1.9, 'original winding', fontsize=6.5, color=MUTED, rotation=90, ha='right', va='center')
axs[0].legend(fontsize=7, loc='upper left', bbox_to_anchor=(0.1, 1.0), title='coil pitch', title_fontsize=7)
axs[1].text(0.795, 0.2, 'bands: MLT +5.4…10.8 %\nper slot of pitch', fontsize=6, ha='right', color=MUTED)
fig.savefig('../figures/fig_pitch_variants.png', bbox_inches='tight'); fig.savefig('../figures/fig_pitch_variants.pdf', bbox_inches='tight')
