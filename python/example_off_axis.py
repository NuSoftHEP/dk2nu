import pydk2nu as dk

import sys
import time

from matplotlib import pyplot as plt
import numpy as np
from math import pi

fr = dk.dk2nuFileReader(sys.argv[1:])

POT = np.sum([m.pots for m in fr.metas()])

print(f"{POT} Protons on target") 

ND_Z_cm = 57400
offaxis_positions_m = [0, 5, 10, 20, 25]

nus = fr.decay_all_through_points([ [x*100, 0, ND_Z_cm] for x in offaxis_positions_m ])

numus = nus[nus[:,0] == 14]

for i, oam in enumerate(offaxis_positions_m):
  
  offaxis_label = f"{oam} m Off Axis" if oam > 0 else "On Axis"
  plt.hist(numus[:,1 + i*2], bins=np.linspace(0,10,200), 
           weights=numus[:,2 + i*2]/(POT*10.0/200.0), 
           label=f"ND Flux ({offaxis_label})")

plt.legend()
plt.ylabel(r"flux /cm${}^2$ /GeV /POT")
plt.xlabel(r"$E_\nu$ (GeV)")
plt.savefig("flux_ND.pdf", bbox_inches='tight')