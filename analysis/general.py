"""Generalisation test of the shop rule: synthetic loss-split variants of the calibrated hybrid model."""
import json, math, numpy as np
import hybrid2 as H
Rs_nom = H.Rs_of(18, 2, 0.75)

def eta_at(a, mlt, kfe_total, rs_scale):
    return H.at_power(3700, Rs_nom * rs_scale * mlt / a, k_fe=kfe_total)

def required_area(dk, mlt, kb, rs_scale, eta0):
    lo, hi = 0.8, 2.5
    for _ in range(40):
        m = 0.5 * (lo + hi)
        if eta_at(m, mlt, kb * (1 + dk), rs_scale)['eta'] < eta0: lo = m
        else: hi = m
    return 0.5 * (lo + hi)

rows = []
for family, grid in (('core', [(kb, 1.0) for kb in (0.5, 1.0, 1.6, 2.6, 3.6, 5.0, 6.5)]),
                     ('copper', [(1.0, rs) for rs in (1.6, 1.0, 0.7, 0.5, 0.35, 0.25)])):
    for kb, rs in grid:
        b = eta_at(1.0, 1.0, kb, rs); eta0 = b['eta']; ratio = b['P_fe'] / b['P_cu1']
        for dk in (0.2, 0.5, 1.0):
            for mlt in (1.0, 1.05):
                a_mod = required_area(dk, mlt, kb, rs, eta0)
                a_rule = mlt / (1 - dk * ratio) if dk * ratio < 0.95 else float('nan')
                rows.append(dict(family=family, kb=kb, rs=rs, ratio=ratio, eta0=eta0, dk=dk, mlt=mlt, a_mod=a_mod, a_rule=a_rule,
                                 err=100 * (a_rule - a_mod) / a_mod))
        print('%-6s ratio=%.2f eta0=%.2f%% | ' % (family, ratio, 100 * eta0) +
              '  '.join('dk%.1f m%.2f: mod %.3f rule %.3f (%+.1f%%)' % (r['dk'], r['mlt'], r['a_mod'], r['a_rule'], r['err'])
                        for r in rows[-6:] if r['mlt'] == 1.0))
json.dump(rows, open('general.json', 'w'))
