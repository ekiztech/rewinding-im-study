"""Delta k_fe,max (first-order criterion, MLT unchanged) with input uncertainty: stator-copper and core-loss
shares each lognormal with 15 % scatter, and the loop-test reading related to the whole-core loss by a factor 0.7-1.3."""
import json, numpy as np
rng = np.random.default_rng(2026); N = 20000
G = json.load(open('ecgen.json')); out = {}
for g in G + [dict(name='Reference 3.7 kW', r=0.192)]:
    r = g['r']*np.exp(rng.normal(0, 0.15, N) - rng.normal(0, 0.15, N)); k = rng.uniform(0.7, 1.3, N)
    res = {}
    for h in (1.05, 1.10):
        dk = (1 - 1/h)/r/k                                            # allowed increase as read by the core test
        res[str(h)] = dict(nominal=(1-1/h)/g['r'], p05=float(np.percentile(dk, 5)), p50=float(np.median(dk)), p95=float(np.percentile(dk, 95)))
    out[g['name']] = res
    print('%-17s ' % g['name'] + ' | '.join('+%d %% Cu: nominal %.1f %%, 90 %% interval %.1f-%.1f %%' % (round(100*(h-1)), 100*res[str(h)]['nominal'], 100*res[str(h)]['p05'], 100*res[str(h)]['p95']) for h in (1.05, 1.10)))
json.dump(out, open('dkmax_unc.json', 'w'), indent=1)
