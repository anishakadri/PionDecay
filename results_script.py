import Particles as pa
import chamber as ch
import numpy as np
from matplotlib import pyplot as plt
import pylab as pl

"""
A program to collect data from the chamber.
"""
c = ch.Chamber()

for i in [800, 4900, 9000]:
    c.muelectron_scatter(i,1000)
    c.pielectron_scatter(i, 1000)
    c.muon_scatter(i,1000)
    c.pion_scatter(i,1000)
    c.electron_scatter(400,1000)

E= np.linspace(500,10000,100)
ratio=np.array([])

for i in E:
   ratio = np.append(ratio, c.escape_ratio(i))
fig = plt.figure()  
pl.plot(E,ratio,'b-')
plt.title("Ratio of muons leaving chamber at different energies")
pl.xlabel("Energy of Pion source/ MeV")
pl.ylabel("Escape ratio")
pl.show()
