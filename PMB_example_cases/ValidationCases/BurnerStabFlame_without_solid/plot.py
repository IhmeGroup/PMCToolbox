import matplotlib.pyplot as plt
import numpy as np

ref = np.genfromtxt("./reference_data.txt",names=True)
pmb = np.genfromtxt("./result_phi_1.00_fuel_NH3_mdot_0.02.txt",names=True)

xref = 0
xpmb = 0
plt.plot(ref["x"]-xref,ref["T"],'-r',label="Cantera")
plt.plot(pmb["x"]-xpmb,pmb["Tg"],'--k',label="CanteraPMB")
plt.ylabel("T (K)")
plt.xlabel("x (m)")
plt.legend()
plt.xlim([0,0.02])
plt.savefig("T.png")
plt.close("all")

plt.plot(ref["x"]-xref,ref["u"],'-r',label="Cantera")
plt.plot(pmb["x"]-xpmb,pmb["u"],'--k',label="CanteraPMB")
plt.ylabel("u (m/s)")
plt.xlabel("x (m)")
plt.legend()
plt.savefig("u.png")
plt.close("all")

