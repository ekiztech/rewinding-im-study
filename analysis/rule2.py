"""Shop rule v1 (first-order loss balance) and v2 (adds the flux feedback of the reduced stator voltage drop),
tested on synthetic loss splits of the hybrid model. Required area from the model by bisection (upper bound 4)."""
import json, math, cmath
import hybrid2 as H
Rs_nom = H.Rs_of(18, 2, 0.75)

def run(a, mlt, kfe, rs, P2=3700):
    return H.at_power(P2, Rs_nom * rs * mlt / a, k_fe=kfe)

def a_model(dk, mlt, kb, rs, eta0, P2=3700):
    lo, hi = 0.8, 4.0
    if run(hi, mlt, kb*(1+dk), rs, P2)['eta'] < eta0: return float('nan')
    for _ in range(45):
        m = 0.5*(lo+hi)
        if run(m, mlt, kb*(1+dk), rs, P2)['eta'] < eta0: lo = m
        else: hi = m
    return 0.5*(lo+hi)

def rule_v1(dk, r, mlt): 
    x = dk*r
    return mlt/(1-x) if x < 0.95 else float('nan')

def rule_v2(dk, r, r2, delta, mlt):
    """u = 1 - MLT/a ; loss balance incl. airgap-EMF rise eps = delta*u -> core +2eps, rotor Cu -2eps."""
    den = 1 - 2*delta*((1+dk)*r - r2)
    # with MLT: saving term uses (1 - mlt/a); extra Cu from MLT already in a = mlt/(1-u)
    u = dk*r/den
    return mlt/(1-u) if 0 < u < 0.95 else float('nan')

def base_quant(b, rs):
    Rs = Rs_nom*rs
    # in-phase stator drop fraction: delta = I1*Rs*cos(phi)/V
    return dict(r=b['P_fe']/b['P_cu1'], r2=b['P_cu2']/b['P_cu1'], delta=b['I1']*Rs*b['pf']/H.V_PH, eta0=b['eta'])

if __name__ == '__main__':
    rows = []
    grid = [('core', kb, 1.0) for kb in (0.5, 1.0, 1.6, 2.6, 3.6, 5.0, 6.5)] + \
           [('copper', 1.0, rs) for rs in (1.6, 0.7, 0.5, 0.35, 0.25)]          # rs=1.0 is the reference (in 'core')
    for fam, kb, rs in grid:
        b = run(1.0, 1.0, kb, rs); q = base_quant(b, rs)
        for dk in (0.2, 0.5, 1.0):
            for mlt in (1.0, 1.05):
                am = a_model(dk, mlt, kb, rs, q['eta0'])
                a1 = rule_v1(dk, q['r'], mlt); a2 = rule_v2(dk, q['r'], q['r2'], q['delta'], mlt)
                rows.append(dict(family=fam, kb=kb, rs=rs, dk=dk, mlt=mlt, a_mod=am, a1=a1, a2=a2,
                                 e1=100*(a1-am)/am, e2=100*(a2-am)/am, Pfe_share=b['P_fe']/b['P1'], **q))
    json.dump(rows, open('rule2.json', 'w'))
    import numpy as np
    def stats(sel, key):
        e = np.array([r[key] for r in sel if np.isfinite(r[key])]); 
        return len(e), np.mean(np.abs(e)), np.max(np.abs(e)), np.min(e), np.max(e)
    for lbl, f in (('rule-pred v1<=1.5', lambda r: np.isfinite(r['a1']) and r['a1'] <= 1.5 and np.isfinite(r['a_mod'])),
                   ('rule-pred v2<=1.5', lambda r: np.isfinite(r['a2']) and r['a2'] <= 1.5 and np.isfinite(r['a_mod'])),
                   ('all finite', lambda r: np.isfinite(r['a_mod']))):
        sel = [r for r in rows if f(r)]
        print(lbl, 'v1 n=%d mean|e|=%.2f max|e|=%.2f range %.1f..%.1f' % stats(sel, 'e1'),
              '| v2 n=%d mean|e|=%.2f max|e|=%.2f range %.1f..%.1f' % stats(sel, 'e2'))
    for r in sorted(rows, key=lambda r: r['dk']*r['r']):
        if r['mlt'] == 1.0:
            print('%-6s r=%.2f r2=%.2f d=%.3f dk=%.1f mod=%.3f v1=%.3f(%+.1f) v2=%.3f(%+.1f) eta0=%.1f' % (
                r['family'], r['r'], r['r2'], r['delta'], r['dk'], r['a_mod'], r['a1'], r['e1'], r['a2'], r['e2'], 100*r['eta0']))
