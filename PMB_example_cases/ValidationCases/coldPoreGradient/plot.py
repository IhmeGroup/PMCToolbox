import matplotlib.pyplot as plt
import numpy as np

pmb = np.genfromtxt("./result_phi_0.00_fuel_NH3_mdot_0.08.txt",names=True)

xref = 0
xpmb = 0
plt.plot(pmb["x"]-xpmb,pmb["u"],'--k',label="u (m/s)")
plt.plot(pmb["x"]-xpmb,pmb["porosity"],':g',label="porosity")
plt.plot(pmb["x"]-xpmb,pmb["u"]*pmb["porosity"],'-r',label="u*porosity (m/s)")
plt.ylabel("quantitiy")
plt.xlabel("x (m)")
plt.legend()
plt.savefig("u.png")
plt.close("all")

