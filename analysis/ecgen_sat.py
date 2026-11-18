"""IE2/IE3 motors of ecgen.py with a saturating magnetizing branch, I_m = I_m0 (E/E0)^mu (mu = 1.6 as the
reference motor at rated load, and 2.5): does the compensation need or the criterion accuracy change?"""
import math, cmath, json
import ecgen as E
from rule3 import rule_v3
from rule2 import rule_v1
class SatMotor(E.Motor):
    mu = 1.0
    def solve(self, s, Rs, kfe=1.0):
        Zr = self.Rr/s + 1j*self.Xr; Xm = self.Xm
        for _ in range(200):
            Zp = 1/(1/(1j*Xm) + 1/Zr); I1 = E.V/(Rs + 1j*self.Xs + Zp); Em = abs(I1*Zp)
            Xn = self.Xm*(Em/self.E0)**(1 - self.mu)              # I_m = (E/Xm0)(E/E0)^(mu-1)
            if abs(Xn - Xm) < 1e-9*self.Xm: break
            Xm = 0.5*(Xm + Xn)
        self._Xm_op = Xm
        Em_ = I1*Zp; I2 = Em_/Zr
        Pfe = kfe*self.Pfe0*(abs(Em_)/self.E0)**2
        Pag = 3*abs(I2)**2*self.Rr/s; n = (1-s)*60*E.F/self.p
        Pfw = self.Pfw0*(n/self.n0)**2; PLL = self.PLL0*(abs(I2)/self.I20)**2
        P2 = Pag*(1-s) - Pfw - PLL; P1 = 3*(E.V*I1.conjugate()).real + Pfe
        return dict(P2=P2, P1=P1, eta=P2/P1, I1=abs(I1), pf=math.cos(cmath.phase(I1)), P_cu1=3*abs(I1)**2*Rs, P_cu2=s*Pag, P_fe=Pfe, s=s)
out = []
for mu in (1.0, 1.6, 2.5):
    SatMotor.mu = mu
    for cls, P, eta, pf in E.MACH:
        m = SatMotor(cls, P, eta, pf); q = E.quant(m); b = m.at(m.P2)
        row = dict(name=m.name, mu=mu, eta0=q['eta0'], r=q['r'])
        for dk in (0.2, 0.5):
            am = E.a_req(m, dk); v1 = rule_v1(dk, q['r'], 1.0); v3 = rule_v3(dk, q['r'], q['r2'], q['delta'], 1-b['pf']**2, mu)
            row['dk%.1f' % dk] = dict(mod=am, v1=v1, v3=v3, e1=100*(v1-am)/am, e3=100*(v3-am)/am)
        row['dkmax105'] = E.dk_max(m, 1.05)
        out.append(row)
        print('mu %.1f %-12s A_req +20%%: %.3f (v1 %+.1f%%, v3 %+.1f%%)  +50%%: %.3f (v1 %+.1f%%, v3 %+.1f%%)  dk_max(+5%% Cu) %.1f%%' % (
            mu, m.name, row['dk0.2']['mod'], row['dk0.2']['e1'], row['dk0.2']['e3'], row['dk0.5']['mod'], row['dk0.5']['e1'], row['dk0.5']['e3'], 100*row['dkmax105']))
json.dump(out, open('ecgen_sat.json', 'w'), indent=1)
