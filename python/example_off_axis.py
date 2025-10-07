import pydk2nu as dk

import sys

from matplotlib import pyplot as plt
import numpy as np
from math import pi

fr = dk.dk2nuFileReader(sys.argv[1:])

POT = np.sum([m.pots for m in fr.metas()])

print(f"{POT} Protons on target") 

ND_Z_cm = 57400

for offaxis_m in [0, 5, 10, 20, 25]:
  nus = np.array([(dk.decay.ntype,*dk.decay_through_point([offaxis_m*100, 0, ND_Z_cm])) for dk in fr.decays()])
  numus = nus[nus[:,0] == 14]

  offaxis_label = f"{offaxis_m} m Off Axis" if offaxis_m > 0 else "On Axis"
  plt.hist(numus[:,1], bins=np.linspace(0,10,50), weights=numus[:,2]/(POT*10.0/50.0), label=f"ND Flux ({offaxis_label})")

plt.legend()
plt.ylabel(r"flux /cm${}^2$ /GeV /POT")
plt.xlabel(r"$E_\nu$ (GeV)")
plt.savefig("flux_ND.pdf", bbox_inches='tight')