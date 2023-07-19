import matplotlib.pyplot as plt
import numpy as np
import cantera as ct

adiabaticFlame = np.genfromtxt("result_phi_1.00_fuel_H2_mdot_2.00.txt",unpack=True)
#adiabaticFlame = np.genfromtxt("result_phi_1.00_fuel_H2_mdot_10.00.txt",unpack=True)
gas = ct.Solution("data/H2_K_10_12_0_OD.yaml")
gas.set_equivalence_ratio(1.0,"H2:1","O2:0.21,N2:0.79")
gas.TP = 300., 101325.;
gas.equilibrate("HP")
Tad = gas.T;

plt.plot([0,np.max(adiabaticFlame[0])],[Tad,Tad],'--k',label="Tad")
plt.plot(adiabaticFlame[0],adiabaticFlame[1],'-r',label="Tg")
plt.plot(adiabaticFlame[0],adiabaticFlame[2],'-b',label="Ts")
plt.xlabel("x (m)")
plt.ylabel("T (K)")
plt.legend()
plt.savefig("T.png")
plt.close("all")

plt.plot(adiabaticFlame[0],adiabaticFlame[4],'-k',label="Ts")
plt.xlabel("x (m)")
plt.ylabel("porosity")
plt.savefig("porosity.png")



