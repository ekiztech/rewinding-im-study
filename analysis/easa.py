"""Measured before/after-rewind loss segregations (IEEE 112B) from the EASA/AEMT 2019 report, Tables 1-2/1-4
(10 IE3 motors of the 2019 study, 11 motors of the 2003 study). Loss ratio r, core-loss change, and how the
stator copper loss followed the corrected winding resistance."""
import csv, json, numpy as np
R = list(csv.DictReader(open('../data_ext/easa_rewind_losses.csv')))
out = []
for x in R:
    v = {k: (float(x[k]) if k not in ('study', 'motor') else x[k]) for k in x}
    P2 = v['rating_kW']*1e3; P1 = P2/(v['eta_pre']/100)
    d = dict(study=v['study'], motor=v['motor'], kW=v['rating_kW'], r=v['Pfe_pre']/v['Pcu1_pre'], r2=v['Pcu2_pre']/v['Pcu1_pre'],
             dk=v['Pfe_post']/v['Pfe_pre'] - 1, dR=v['R_corr_post']/v['R_corr_pre'] - 1, dPcu1=v['Pcu1_post']/v['Pcu1_pre'] - 1,
             bal=(v['Pcu1_post'] - v['Pcu1_pre'] + v['Pfe_post'] - v['Pfe_pre'])/P1*100,
             other=(sum(v[k+'_post'] - v[k+'_pre'] for k in ('Pcu2', 'Pfw', 'Pll')))/P1*100,
             deta=v['eta_post'] - v['eta_pre'], faulty=(x['motor'] == '14H'))
    out.append(d)
def s(key, sel): a = np.array([o[key] for o in sel]); return np.median(a), a.min(), a.max()
for st in ('2019', '2003'):
    sel = [o for o in out if o['study'] == st and not o['faulty']]
    print('%s (n=%d): r median %.2f [%.2f, %.2f]; core-loss change median %+.1f%% [%+.1f, %+.1f]' % ((st, len(sel)) + s('r', sel) + tuple(100*np.array(s('dk', sel)))))
ok = [o for o in out if not o['faulty']]
print('all (n=%d): r median %.2f [%.2f, %.2f]' % ((len(ok),) + s('r', ok)))
print('core-loss increase > 20 %%: %s' % [(o['motor'], round(100*o['dk'])) for o in out if o['dk'] > 0.2])
# stator loss vs resistance: lower resistance -> stator loss falls by more than resistance (current and temperature fall)
big = [o for o in ok if abs(o['dR']) >= 0.03]
for o in sorted(big, key=lambda o: o['dR']):
    print('  %-4s %s  dR %+5.1f%%  dPcu1 %+5.1f%%  ratio %.2f' % (o['motor'], o['study'], 100*o['dR'], 100*o['dPcu1'], o['dPcu1']/o['dR']))
dn = [o for o in big if o['dR'] < 0]; up = [o for o in big if o['dR'] > 0]
print('R down (n=%d): mean dPcu1/dR = %.2f ; R up (n=%d): mean dPcu1/dR = %.2f' % (len(dn), np.mean([o['dPcu1']/o['dR'] for o in dn]), len(up), np.mean([o['dPcu1']/o['dR'] for o in up])))
# how much of the efficiency change does the stator+core balance explain?
b = np.array([o['bal'] for o in ok]); oth = np.array([o['other'] for o in ok]); de = np.array([o['deta'] for o in ok])
print('corr(balance, -deta) = %.2f ; corr(other losses, -deta) = %.2f ; corr(balance+other, -deta) = %.2f' % (
    np.corrcoef(b, -de)[0, 1], np.corrcoef(oth, -de)[0, 1], np.corrcoef(b + oth, -de)[0, 1]))
print('std of balance term %.2f %%, of other losses %.2f %% of input' % (b.std(), oth.std()))
json.dump(out, open('easa.json', 'w'), indent=1, default=float)
