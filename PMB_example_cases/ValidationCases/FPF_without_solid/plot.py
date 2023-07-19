import matplotlib.pyplot as plt
import numpy as np

ref = np.genfromtxt("./reference_data.txt",names=True)
pmb = np.genfromtxt("./result_phi_1.00_fuel_NH3_mdot_0.30.txt",names=True)

xref = ref["x"][np.argmax(ref["H"])]
xpmb = pmb["x"][np.argmax(pmb["H"])]
plt.plot(ref["x"]-xref,ref["T"],'-r',label="Cantera")
plt.plot(pmb["x"]-xpmb,pmb["Tg"],'--k',label="CanteraPMB")
plt.ylabel("T (K)")
plt.xlabel("x (m)")
plt.xlim([-0.01,0.01])
plt.legend()
plt.savefig("T.png")
plt.close("all")

plt.plot(ref["x"]-xref,ref["u"],'-r',label="Cantera")
plt.plot(pmb["x"]-xpmb,pmb["u"],'--k',label="CanteraPMB")
plt.ylabel("u (m/s)")
plt.xlabel("x (m)")
plt.legend()
plt.savefig("u.png")
plt.close("all")

