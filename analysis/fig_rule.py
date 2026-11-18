import json, math, itertools, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import hybrid2 as H, ecgen as E
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':9,'axes.linewidth':0.6,'axes.spines.top':False,'axes.spines.right':False,
    'axes.grid':True,'grid.color':'#e3e2dd','grid.linewidth':0.5,'savefig.dpi':300,'legend.frameon':False})
INK, MUTED = '#22221f', '#7a7970'
C4 = ['#2a78d6', '#eb6834', '#1baf7a', '#8a5cd0']
FIG = '/home/claude/release/rewinding-im-study/figures/'
def save(fig, name):
    fig.savefig(FIG+name+'.png', bbox_inches='tight'); fig.savefig(FIG+name+'.pdf', bbox_inches='tight')

# ---------- Fig. 8: rule accuracy and compensable damage
R3 = json.load(open('rule3.json'))
S = [r for r in R3['synthetic'] if np.isfinite(r['a_mod']) and np.isfinite(r['a1'])]
G = json.load(open('ecgen.json')); RR = json.load(open('refrev.json'))
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={'wspace':0.3})
a = ax[0]
a.axhspan(-2, 2, color='#efeee9', lw=0)
a.scatter([r['dk']*r['r'] for r in S], [r['e1'] for r in S], s=14, facecolor='none', edgecolor=MUTED, lw=0.7, label='first-order rule, Eq. (5)')
a.scatter([r['dk']*r['r'] for r in S if np.isfinite(r['e3'])], [r['e3'] for r in S if np.isfinite(r['e3'])], s=14, color=INK, lw=0, label='with second-order terms, Eq. (6)')
names = [g['name'] for g in G]
for x in R3['ie']:
    i = names.index(x['name']); r = G[i]['r']
    a.scatter(x['dk']*r, x['e1'], marker='D', s=18, facecolor='white', edgecolor=C4[i], lw=0.9, zorder=3)
    a.scatter(x['dk']*r, x['e3'], marker='D', s=18, color=C4[i], edgecolor='white', lw=0.4, zorder=3, label=x['name'].replace('.0 kW', ' kW') if x['dk'] == 0.2 else None)
for x in R3['partload']:
    if x['load'] == '1.0': continue
    rr = RR['partload'][x['load']]['r']; dk = {'D1':0.2,'D2':0.5,'D3':0.2}[x['dmg']]
    a.scatter(dk*rr, x['e1'], marker='^', s=20, facecolor='white', edgecolor=MUTED, lw=0.8, zorder=3)
    a.scatter(dk*rr, x['e3'], marker='^', s=20, color=INK, zorder=3, label='reference motor, 50/75 % load' if (x['load'], x['dmg']) == ('0.5', 'D1') else None)
a.axhline(0, color=MUTED, lw=0.6); a.set_ylim(-22, 16); a.set_xlim(0, 0.75)
a.set_xlabel(r'$\Delta k_{fe}\,r$'); a.set_ylabel('Criterion error vs. model (%)')
a.text(0.74, -3, 'conservative above 0', fontsize=6.5, ha='right', va='top', color=MUTED)
a.legend(fontsize=6, loc='upper center', bbox_to_anchor=(0.5,-0.2), ncol=2); a.set_title('(a) Criterion error', fontsize=9, loc='left')
b = ax[1]; hd = np.linspace(1.0, 1.25, 11)
ref = []
Rs0 = H.Rs_of(18, 2, 0.75); e0 = H.at_power(3700, Rs0)['eta']
for h in hd:
    lo, hi = 0.0, 3.0
    for _ in range(30):
        d = 0.5*(lo+hi)
        if H.at_power(3700, Rs0/h, k_fe=1+d)['eta'] < e0: hi = d
        else: lo = d
    ref.append(100*0.5*(lo+hi))
b.plot(100*(hd-1), ref, color=INK, lw=1.6, label='reference 3.7 kW (r = 0.19)')
mach = [E.Motor(*m) for m in E.MACH]
for i, m in enumerate(mach):
    b.plot(100*(hd-1), [100*E.dk_max(m, h) for h in hd], color=C4[i], lw=1.3, label='%s (r = %.2f)' % (m.name.replace('.0 kW',' kW'), G[i]['r']))
b.axhline(20, color=INK, ls='--', lw=0.8); b.text(24.5, 22, 'EASA limit', fontsize=6.5, ha='right')
b.set_xlim(0, 25); b.set_ylim(0, 120); b.set_xlabel('Copper-area headroom of the rewind (%)'); b.set_ylabel('Max. compensable core-loss increase (%)')
b.legend(fontsize=6, loc='upper left'); b.set_title('(b) Damage a rewind can absorb', fontsize=9, loc='left')
save(fig, 'fig_rule_generality')

import sys; sys.exit(0)
# ---------- Fig. 6: pitch comparison vs insulated lay factor, MLT band
import refrev_lib as L
fig, ax = plt.subplots(figsize=(3.6, 2.8))
ks = np.arange(0.55, 0.801, 0.01); COL = {7:C4[0], 8:C4[1], 9:C4[2]}
env = {}
for inc in (0.054, 0.108):
    rows = L.scan(0.553*inc)
    for p in (7, 8, 9):
        y = []
        for k in ks:
            b_ = L.best(rows, p, 'D1', k); y.append(b_['deta'] if b_ else np.nan)
        env[(p, inc)] = np.fmax.accumulate(np.nan_to_num(y, nan=-9))
for p in (7, 8, 9):
    lo_, hi_ = env[(p, 0.108)], env[(p, 0.054)]
    if p != 7: ax.fill_between(ks, lo_, hi_, step='post', color=COL[p], alpha=0.18, lw=0)
    ax.step(ks, hi_, where='post', color=COL[p], lw=1.4, label='%d/9' % p)
ax.axvline(L.kins(18, 2, 0.75), color=MUTED, lw=0.7, ls=':'); ax.text(L.kins(18,2,0.75)-0.004, 1.9, 'original winding', fontsize=6.5, color=MUTED, rotation=90, ha='right', va='center')
ax.set_xlim(0.55, 0.80); ax.set_ylim(0, 3); ax.set_xlabel('Allowed insulated lay factor $k_{ins}$'); ax.set_ylabel(r'Best $\Delta\eta$ vs. undamaged original (pp)')
ax.legend(fontsize=7, loc='upper left', bbox_to_anchor=(0.1, 1.0), title='coil pitch', title_fontsize=7)
ax.text(0.795, 0.25, 'bands: MLT +5.4…10.8 %\nper slot of pitch', fontsize=6, ha='right', color=MUTED)
save(fig, 'fig_pitch_variants')

# ---------- Fig. 9: efficiency vs load
fig, ax = plt.subplots(figsize=(3.6, 2.7)); Ld = 100*np.array(RR['loads']); C = RR['curves']
ax.plot(Ld, C['orig_D0'], color=INK, lw=1.5, label='original, undamaged')
ax.plot(Ld, C['orig_D1'], color=C4[1], lw=1.2, ls='--', label='original, core +20 %')
ax.plot(Ld, C['2x0.85_D1'], color=C4[1], lw=1.4, label='2×0.85 mm, core +20 %')
ax.plot(Ld, C['orig_D2'], color=C4[0], lw=1.2, ls='--', label='original, core +50 %')
ax.plot(Ld, C['2x0.85_D2'], color=C4[0], lw=1.4, label='2×0.85 mm, core +50 %')
ax.set_xlabel('Load (% of 3.7 kW)'); ax.set_ylabel('Efficiency (%)'); ax.set_xlim(25, 125); ax.legend(fontsize=6.5, loc='lower right')
save(fig, 'fig_partload')

# ---------- Fig. 7: Monte Carlo v2
M = json.load(open('mc2.json')); from scipy.stats import spearmanr
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8), gridspec_kw={'wspace':0.45, 'width_ratios':[1.1, 1]})
a = ax[0]
for i, (k, lab) in enumerate((('a_D1','core +20 %'), ('a_D3','+20 %, MLT +5 %'), ('a_D2','core +50 %'), ('a_D1_shop','+20 %, incl. loop-test mapping'))):
    v = np.array([o[k] for o in M]); lo_, med, hi_ = np.percentile(v, [2.5, 50, 97.5])
    a.plot([lo_, hi_], [i, i], color=C4[i], lw=5, alpha=0.35, solid_capstyle='butt'); a.plot(med, i, 'o', color=C4[i], ms=5)
    a.text(hi_+0.004, i, '%.3f [%.3f, %.3f]' % (med, lo_, hi_), va='center', fontsize=6.5)
a.set_yticks(range(4)); a.set_yticklabels(['core +20 %', '+20 %, MLT +5 %', 'core +50 %', '+20 %, shop\nmeasurement'], fontsize=7)
a.set_xlim(1.0, 1.24); a.set_xlabel(r'Required $A_{Cu}/A_{Cu,0}$ (median, 95 % interval)'); a.invert_yaxis(); a.grid(axis='y', visible=False)
a.set_title('(a) Compensation need, N = 500', fontsize=9, loc='left')
b = ax[1]; names = [('cI','stator Cu-loss level'), ('cfe','core-loss data'), ('T','temperature'), ('cll','add. load loss'), ('cr','rotor resistance'),
                    ('cm','magnetizing curve'), ('cl','2D leakage'), ('cfw','friction'), ('split','3D-inductance split')]
rho = [spearmanr([o[k] for o in M], [o['a_D1'] for o in M])[0] for k, _ in names]
b.barh(range(len(names)), rho, color=[INK if abs(r) > 0.3 else MUTED for r in rho], height=0.6)
b.set_yticks(range(len(names))); b.set_yticklabels([n for _, n in names], fontsize=7); b.invert_yaxis(); b.axvline(0, color=INK, lw=0.6)
b.set_xlim(-0.8, 0.8); b.set_xlabel(r'Spearman $\rho$ with required area (core +20 %)'); b.grid(axis='y', visible=False)
b.set_title('(b) What drives it', fontsize=9, loc='left')
save(fig, 'fig_uncertainty')
print('figures written')
