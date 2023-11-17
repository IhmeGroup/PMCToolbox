import matplotlib.pyplot as plt
import numpy as np
#x u rho mdot phi T Ts N2 H2 H O2 O H2O OH H2O2 HO2 NO N2O NO2 HNO HNO2 HONO HONO2 N2H2 H2NN NH2OH HNOH NH3 N2H4 N NO3 NH NNH NH2 H2NO N2H3 

new = np.genfromtxt("table.txt",unpack=True)

plt.plot(new[0],new[5],'-r',label="Tg")
plt.plot(new[0],new[6],'--k',label="Ts")
plt.xlabel("x (m)")
plt.ylabel("T (K)")
plt.legend()
plt.savefig("T.png")
plt.close("all")

plt.plot(new[0],new[4],'-r',label="phi")
plt.xlabel("x (m)")
plt.ylabel("phi")
plt.legend()
plt.savefig("phi.png")
plt.close("all")

plt.plot(new[0],new[3],'-r',label="mdot")
plt.xlabel("x (m)")
plt.ylabel("mdot (kg/m2/s)")
plt.legend()
plt.savefig("mdot.png")
plt.close("all")







