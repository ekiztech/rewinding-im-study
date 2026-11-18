"""Generic T-equivalent-circuit motors built from published loss shares (Aarniovuori et al. 2018, IE2/IE3 fits)
and IEC 60034-30-1 minimum efficiencies; used to test the shop rule on machines other than the reference motor."""
import math, cmath, json
from rule2 import rule_v1, rule_v2
V, F = 400/math.sqrt(3), 50.0
FITS = {  # share of total loss (%) = a2*P^2 + a1*P + a0 (P in kW); 'add' = remainder
 'IE2': dict(cu1=(1e-4,-0.0961,32.9719), cu2=(2e-4,-0.0745,20.6439), fe=(-1e-4,0.0385,15.4159), fw=(-1e-4,0.0324,11.9837)),
 'IE3': dict(cu1=(1e-4,-0.0796,32.3678), cu2=(0,-0.0263,17.5878), fe=(0,-0.0067,18.9585), fw=(0,0.0264,10.8933))}
MACH = [('IE2', 7.5, 0.887, 0.83), ('IE3', 7.5, 0.904, 0.83), ('IE3', 37.0, 0.939, 0.86), ('IE3', 110.0, 0.954, 0.87)]
XPU = 0.07   # stator and rotor leakage reactance, per unit of V/I (assumption)

class Motor:
    def __init__(self, cls, P, eta, pf, p=2):
        f = FITS[cls]; sh = {k: (a2*P*P + a1*P + a0)/100 for k, (a2, a1, a0) in f.items()}; sh['ll'] = 1 - sum(sh.values())
        self.name = '%s %.1f kW' % (cls, P); self.P2 = P*1e3; self.sh = sh; self.p = p
        P1 = self.P2/eta; L = P1 - self.P2
        self.Pfe0, self.Pfw0, self.PLL0 = sh['fe']*L, sh['fw']*L, sh['ll']*L
        Pcu1, Pcu2 = sh['cu1']*L, sh['cu2']*L
        I = (P1 - self.Pfe0)/(3*V*pf); I1 = I*cmath.exp(-1j*math.acos(pf))
        self.Rs = Pcu1/(3*I*I); Zb = V/I; self.Xs = self.Xr = XPU*Zb
        Pag = P1 - Pcu1 - self.Pfe0; self.s0 = Pcu2/Pag
        E = V - I1*(self.Rs + 1j*self.Xs); self.E0 = abs(E)
        lo, hi = 0.5*Zb, 50*Zb                     # Xm such that the rotor branch reactance equals Xr
        for _ in range(200):
            Xm = 0.5*(lo+hi); I2 = I1 - E/(1j*Xm); z = E/I2
            if z.imag > self.Xr: hi = Xm
            else: lo = Xm
        self.Xm = Xm; self.Rr = self.s0*(E/I2).real; self.I20 = abs(I2); self.n0 = (1-self.s0)*60*F/p
    def solve(self, s, Rs, kfe=1.0):
        Zr = self.Rr/s + 1j*self.Xr; Zp = 1/(1/(1j*self.Xm) + 1/Zr)
        I1 = V/(Rs + 1j*self.Xs + Zp); E = I1*Zp; I2 = E/Zr
        Pfe = kfe*self.Pfe0*(abs(E)/self.E0)**2
        Pag = 3*abs(I2)**2*self.Rr/s; n = (1-s)*60*F/self.p
        Pfw = self.Pfw0*(n/self.n0)**2; PLL = self.PLL0*(abs(I2)/self.I20)**2
        P2 = Pag*(1-s) - Pfw - PLL; P1 = 3*(V*I1.conjugate()).real + Pfe
        return dict(P2=P2, P1=P1, eta=P2/P1, I1=abs(I1), pf=math.cos(cmath.phase(I1)), P_cu1=3*abs(I1)**2*Rs,
                    P_cu2=s*Pag, P_fe=Pfe, s=s)
    def at(self, P2, Rs=None, kfe=1.0):
        Rs = self.Rs if Rs is None else Rs; lo, hi = 1e-5, 0.3
        for _ in range(80):
            m = 0.5*(lo+hi)
            if self.solve(m, Rs, kfe)['P2'] < P2: lo = m
            else: hi = m
        return self.solve(0.5*(lo+hi), Rs, kfe)

def a_req(m, dk, mlt=1.0, load=1.0):
    P2 = load*m.P2; eta0 = m.at(P2)['eta']; lo, hi = 0.8, 4.0
    if m.at(P2, m.Rs*mlt/hi, 1+dk)['eta'] < eta0: return float('nan')
    for _ in range(50):
        a = 0.5*(lo+hi)
        if m.at(P2, m.Rs*mlt/a, 1+dk)['eta'] < eta0: lo = a
        else: hi = a
    return 0.5*(lo+hi)

def dk_max(m, a, load=1.0):
    P2 = load*m.P2; eta0 = m.at(P2)['eta']; lo, hi = 0.0, 3.0
    for _ in range(50):
        d = 0.5*(lo+hi)
        if m.at(P2, m.Rs/a, 1+d)['eta'] < eta0: hi = d
        else: lo = d
    return 0.5*(lo+hi)

def quant(m, load=1.0):
    b = m.at(load*m.P2)
    return dict(r=b['P_fe']/b['P_cu1'], r2=b['P_cu2']/b['P_cu1'], delta=b['I1']*m.Rs*b['pf']/V, eta0=b['eta'], s=b['s'], I1=b['I1'])

if __name__ == '__main__':
    out = []
    for cls, P, eta, pf in MACH:
        m = Motor(cls, P, eta, pf); q = quant(m)
        row = dict(name=m.name, eta0=q['eta0'], s=q['s'], r=q['r'], r2=q['r2'], delta=q['delta'], shares=m.sh)
        for dk in (0.2, 0.5):
            for mlt in (1.0, 1.05):
                am = a_req(m, dk, mlt); a1 = rule_v1(dk, q['r'], mlt); a2 = rule_v2(dk, q['r'], q['r2'], q['delta'], mlt)
                row['dk%.1f_m%.2f' % (dk, mlt)] = dict(mod=am, v1=a1, v2=a2)
        row['dkmax'] = {str(a): dk_max(m, a) for a in (1.05, 1.10, 1.20)}
        for load in (0.5, 0.75):
            ql = quant(m, load); am = a_req(m, 0.2, 1.0, load)
            row['load%.2f' % load] = dict(r=ql['r'], mod=am, v1=rule_v1(0.2, ql['r'], 1.0), v2=rule_v2(0.2, ql['r'], ql['r2'], ql['delta'], 1.0))
        out.append(row)
        print('%-12s eta0 %.2f%% s %.4f  r %.3f r2 %.3f delta %.3f | ' % (m.name, 100*q['eta0'], q['s'], q['r'], q['r2'], q['delta']) +
              ' '.join('%s mod %.3f v1 %.3f v2 %.3f;' % (k, v['mod'], v['v1'], v['v2']) for k, v in row.items() if k.startswith('dk0')))
        print('     dk_max @ +5/10/20%% copper: %s' % ', '.join('%.1f%%' % (100*v) for v in row['dkmax'].values()),
              '| 50%% load r %.2f mod %.3f v2 %.3f | 75%% load r %.2f mod %.3f v2 %.3f' % (
              row['load0.50']['r'], row['load0.50']['mod'], row['load0.50']['v2'], row['load0.75']['r'], row['load0.75']['mod'], row['load0.75']['v2']))
    json.dump(out, open('ecgen.json', 'w'), indent=1)
