import photovoltaic as pv
import numpy as np
import matplotlib.pyplot as plt
Sn = 10
Sp = 10
JL = 0.036
W = 0.01
q = pv.q
fig, ax = plt.subplots()
dopings = [1e13, 3e13, 1e14, 3e14, 1e15, 3e15, 1e16]
for doping in dopings:
    n0, p0 = pv.carrier_conc(doping)
    #print(n0,p0)
    Δn = np.logspace(10,17,100)
    n = n0 + Δn
    p = p0 + Δn
    U_surface = pv.recombination_surface(n,p,Sn,Sp)
    U_auger = pv.Uauger(n0, p0, Δn, 8.5e9)
    U_SRH = pv.recombination_SRH(n, p, 0.0, 5e-3, 5e-3)

    #print(Usurface)
    JSRH = q * W * U_SRH
    J_auger = q * W * U_auger
    J_surface = q * U_surface
    Jrec = JSRH + J_auger + J_surface
    V = pv.impliedV(Δn, doping)
    J = JL - Jrec
    V[0] = 0
    J[0] = JL
    Voc, Isc, FF, Vmp, Imp = pv.cell_params(V,J)
    print('{:.2e}'.format(doping),Voc, Isc, FF, Vmp, Imp, Vmp*Imp)
    #ax.axis(V,J, label='{:.2e}'.format(doping))
    ax.axes()

plt.xlim(0,0.8)
plt.ylim(0,0.04)
plt.legend()
plt.xlabel('voltage (V)')
plt.ylabel('current density')
plt.title('Sn= '+str(Sn)+', Sp= '+str(Sp))
plt.show()

plt.figure('surface recombination')
dopings = [1e13, 1e14, 1e15, 1e16]
for doping in dopings:
    n0, p0 = pv.carrier_conc(doping)
    #print(n0,p0)
    Δn = np.logspace(10,17)
    n = n0 + Δn
    p = p0 + Δn
    U_surface = pv.recombination_surface(n,p,Sn,Sp)
    #print(Usurface)
    Seff = U_surface/Δn
    plt.plot(Δn,Seff, label='{:.2e}'.format(doping))
plt.semilogx()
plt.legend()
plt.xlabel('Δn (cm-3)')
plt.ylabel('surface recombination')
plt.title('Sn= '+str(Sn)+', Sp= '+str(Sp))
plt.show()
