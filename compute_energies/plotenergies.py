import numpy as np
import matplotlib.pyplot as plt

energies=np.zeros(30)
energies2=np.zeros(30)
harray=np.zeros(30)
for h_i in xrange(30):
  h=+0.01*h_i+0.01
  print h
  harray[h_i]=h
  en=np.loadtxt("energies{}.txt".format(h))
  energies[h_i]=np.min(en[12:])


energies_exact=np.loadtxt("gs_energies_4states.txt")[:30,0]
energies_exact1=np.loadtxt("gs_energies_hx_hz.txt")[:30,1]
energies_exact2=np.loadtxt("gs_energies_hx_hz.txt")[:30,2]
energies_exact3=np.loadtxt("gs_energies_hx_hz.txt")[:30,3]

#energies=np.hstack((energies,energies2))
print energies
#harray=np.hstack((harray,harray2))
print np.shape(harray)
print np.shape(energies_exact)
plt.figure()
#plt.semilogy(harray,energies-energies_exact1)
plt.plot(harray,energies_exact)
plt.plot(harray,energies_exact1)
plt.plot(harray,energies_exact2)
plt.plot(harray,energies_exact3)
plt.plot(harray,energies,'o')
plt.savefig("energies_L3.png")
plt.show()

plt.figure()
plt.semilogy(harray,energies-energies_exact)
plt.savefig("energy_diff_L3.png")
plt.show()



