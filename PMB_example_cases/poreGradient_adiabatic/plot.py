import matplotlib.pyplot as plt
import numpy as np
import cantera as ct

adiabaticFlame = np.genfromtxt("result_phi_1.30_fuel_NH3_mdot_0.01.txt",unpack=True)
gas = ct.Solution("data/Arunthanayothin.yaml")
gas.set_equivalence_ratio(1.3,"NH3:1","O2:0.21,N2:0.79")
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



