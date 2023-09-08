import matplotlib.pyplot as plt
import numpy as np


new = np.genfromtxt("result_phi_1.30_fuel_NH3:0.7,H2:0.3_mdot_0.40.txt",unpack=True)

plt.plot(new[0],new[1],'-r',label="Tg")
plt.plot(new[0],new[2],'--k',label="Ts")
plt.xlabel("x (m)")
plt.ylabel("T (K)")
plt.legend()
plt.savefig("T.png")




