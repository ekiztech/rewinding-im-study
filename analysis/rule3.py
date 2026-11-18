"""Criterion v3: v2 + rise of the magnetizing current with the air-gap EMF (at break-even the input power and
hence the active current are unchanged, so only the reactive part of I1 changes):
   u = dk r / (1 - 2 delta [ (1+dk) r - r2 + mu sin^2(phi) ]),  a = MLT/(1-u)
mu = dln(I_m)/dln(E) at the operating point (1 for a linear magnetizing branch)."""
import json, math, numpy as np
import hybrid2 as H, rule2 as R2, ecgen as E
def rule_v3(dk, r, r2, delta, sin2, mu, mlt=1.0):
    u = dk*r/(1 - 2*delta*((1+dk)*r - r2 + mu*sin2))
    return mlt/(1-u) if 0 < u < 0.95 else float('nan')
def mu_hybrid(load=1.0, rs=1.0):
    P = H.PITCH[7]; b = H.at_power(3700*load, R2.Rs_nom*rs); psi = b['psi']
    Im = np.interp(psi, P.psi, P.I); dI = np.gradient(np.log(P.I), np.log(P.psi))
    return float(np.interp(psi, P.psi, dI)), b
out = {'synthetic': [], 'ie': [], 'partload': []}
rows = json.load(open('rule2.json'))
cache = {}
for r in rows:
    key = (r['family'], r['kb'], r['rs'])
    if key not in cache:
        H.PFE_C = 1.0
        b = R2.run(1.0, 1.0, r['kb'], r['rs']); P = H.PITCH[7]
        mu = float(np.interp(b['psi'], P.psi, np.gradient(np.log(P.I), np.log(P.psi))))
        cache[key] = (mu, 1 - b['pf']**2)
    mu, s2 = cache[key]
    a3 = rule_v3(r['dk'], r['r'], r['r2'], r['delta'], s2, mu, r['mlt'])
    out['synthetic'].append(dict(r, mu=mu, sin2=s2, a3=a3, e3=100*(a3-r['a_mod'])/r['a_mod']))
for cls, Pk, eta, pf in E.MACH:
    m = E.Motor(cls, Pk, eta, pf); q = E.quant(m); b = m.at(m.P2)
    for dk in (0.2, 0.5):
        am = E.a_req(m, dk); a3 = rule_v3(dk, q['r'], q['r2'], q['delta'], 1-b['pf']**2, 1.0)
        out['ie'].append(dict(name=m.name, dk=dk, mod=am, v1=R2.rule_v1(dk, q['r'], 1.0), v3=a3, e1=100*(R2.rule_v1(dk, q['r'], 1.0)-am)/am, e3=100*(a3-am)/am))
RR = json.load(open('refrev.json'))
for L in ('0.5', '0.75', '1.0'):
    mu, b = mu_hybrid(float(L)); pl = RR['partload'][L]
    for dmg, dk, mlt in (('D1', 0.2, 1.0), ('D2', 0.5, 1.0), ('D3', 0.2, 1.05)):
        a3 = rule_v3(dk, pl['r'], pl['r2'], pl['delta'], 1-b['pf']**2, mu, mlt)
        out['partload'].append(dict(load=L, dmg=dmg, mu=mu, mod=pl[dmg]['mod'], v1=pl[dmg]['v1'], v3=a3, e1=100*(pl[dmg]['v1']-pl[dmg]['mod'])/pl[dmg]['mod'], e3=100*(a3-pl[dmg]['mod'])/pl[dmg]['mod']))
def st(a): a = np.array([x for x in a if np.isfinite(x)]); return '%d cases, mean |e| %.2f %%, range %+.1f..%+.1f %%' % (len(a), np.mean(np.abs(a)), a.min(), a.max())
S = out['synthetic']; sel = [s for s in S if np.isfinite(s['a_mod']) and np.isfinite(s['a3']) and s['a3'] <= 1.5]
print('synthetic, rule-pred <=1.5: v1', st([s['e1'] for s in sel]), '| v2', st([s['e2'] for s in sel]), '| v3', st([s['e3'] for s in sel]))
fin = [s for s in S if np.isfinite(s['a_mod'])]
print('synthetic, all finite: v1', st([s['e1'] for s in fin]), '| v3', st([s['e3'] for s in fin]))
print('IE motors: v1', st([x['e1'] for x in out['ie']]), '| v3', st([x['e3'] for x in out['ie']]))
print('reference part/full load: v1', st([x['e1'] for x in out['partload']]), '| v3', st([x['e3'] for x in out['partload']]))
print('mu at rated (reference):', round(mu_hybrid(1.0)[0], 2), ' at 50 %:', round(mu_hybrid(0.5)[0], 2))
json.dump(out, open('rule3.json', 'w'), indent=1, default=float)
