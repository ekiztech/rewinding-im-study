import json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':9,'axes.linewidth':0.6,'axes.spines.top':False,'axes.spines.right':False,
    'axes.grid':True,'grid.color':'#e3e2dd','grid.linewidth':0.5,'savefig.dpi':300,'legend.frameon':False})
INK, MUTED = '#22221f', '#7a7970'
D = json.load(open('easa.json'))
fig, ax = plt.subplots(figsize=(4.6, 3.2))
r = np.linspace(0.15, 1.15, 200)
for h, ls in ((1.05, ':'), (1.10, '-.'), (1.20, '--')):
    ax.plot(r, 100*(1-1/h)/r, ls, color=MUTED, lw=1.0); ax.text(1.14, 100*(1-1/h)/1.14 + 1.2, '+%d %% Cu' % round(100*(h-1)), fontsize=6.5, ha='right', color=MUTED)
ax.axhline(20, color=INK, lw=0.8, ls='--'); ax.text(0.17, 21.5, 'EASA limit', fontsize=6.5)
ax.axhline(0, color=MUTED, lw=0.6)
for st, mk, c, lab in (('2019', 'D', '#2a78d6', 'IE3 motors, 2019 study'), ('2003', 'o', '#eb6834', '2003 study')):
    s = [o for o in D if o['study'] == st and not o['faulty']]
    ax.scatter([o['r'] for o in s], [100*o['dk'] for o in s], marker=mk, s=24, color=c, edgecolor='white', lw=0.5, label=lab, zorder=3)
f = [o for o in D if o['faulty']][0]
ax.scatter([f['r']], [min(100*f['dk'], 39)], marker='^', s=30, color='white', edgecolor=INK, lw=0.8, label='faulty core (+48 %, clipped)', zorder=3)
ax.scatter([0.192], [0], marker='*', s=70, color=INK, label='reference motor (r = 0.19)', zorder=3)
ax.set_xlim(0.15, 1.15); ax.set_ylim(-10, 40)
ax.set_xlabel(r'Core-to-stator-copper loss ratio $r$ before repair'); ax.set_ylabel('Measured core-loss change (%)')
ax.legend(fontsize=6.5, loc='upper left', bbox_to_anchor=(0.0, 0.93))
fig.savefig('../figures/fig_easa.png', bbox_inches='tight'); fig.savefig('../figures/fig_easa.pdf', bbox_inches='tight')
