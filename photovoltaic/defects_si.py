import photovoltaic as pv
import numpy as np
import matplotlib.pyplot as plt



def lifetime_SRH(n, p, Et, τ_n, τ_p, ni_eff=8.5e9, T=298.15):
    q = pv.q
    k = pv.k
    """Return the shockley read hall recombination cm-3
    given Et (eV) trap level from intrinsic"""
    print(k*T/q)
    n1 = ni_eff * np.exp(k*T * Et /q)
    p1 = ni_eff * np.exp(-k * T * Et /q)
    print(n1)
    print(p1)
    U_SRH = (τ_p * (n + n1) + τ_n * (p + p1))/(n * p - ni_eff ** 2)
    return U_SRH

ND = 1e15
n0, p0 = pv.semi.equilibrium_carrier(ND)

delta_n = np.logspace(13,16)

Nt = 1e11

delta_n = 1e14
n = n0 + delta_n
p = p0 + delta_n
Ei = 0.56
defect = ('Au_a',0.550,1.4E-16,7.6E-15)
defect = ('Fe_a',-0.185,1.3E-14,7.0E-17)
Et = defect[1]
sigma_n = defect[2]
sigma_p = defect[3]
vth = 10966214.78
tau_n = 1/(vth*Nt*sigma_n)
tau_p = 1/(vth*Nt*sigma_p)
print('Et from Ei',Et)
print('Et from Ev',Et+Ei)
print('Et from Ec',Ei-Et)
print('tau_n',tau_n)
print('tau_p',tau_p)
tau = lifetime_SRH(n,p,Et,tau_n, tau_p)
print('lifetime ms',tau*1e3)


