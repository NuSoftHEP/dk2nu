import pydk2nu as dk

import sys

from matplotlib import pyplot as plt
import numpy as np
from math import pi

fr = dk.dk2nuFileReader(sys.argv[1:])

POT = np.sum([m.pots for m in fr.metas()])

print(f"{POT} Protons on target") 

ND_Z_cm = 57400

lswheel = ["solid", "dashed", "dotted"]
cwheel = ['#0077bb', '#ee3377', '#009988', '#ee7733']
for i,offaxis_m in enumerate([0, 5, 10]):
  nus = np.array(
    [(dk.decay.ntype, *dk.decay_through_point([offaxis_m*100, 0, ND_Z_cm]), dk.ancestor[-2].pdg) \
      for dk in fr.decays()] 
      )
  numus = nus[nus[:,0] == 14]
  for j,ppdg in enumerate(np.nditer(np.unique(numus[:,3]))):
    numus_wparent = numus[numus[:,3] == ppdg]

    bin_vals, bins = np.histogram(numus_wparent[:,1], bins=np.linspace(0,10,50), weights=numus_wparent[:,2])

    bcs = (bins[1:] + bins[:-1])/2.0
    bws = (bins[1:] - bins[:-1])

    offaxis_label = f"{offaxis_m} m Off Axis" if offaxis_m > 0 else "On Axis"
    plt.plot(bcs, bin_vals/(POT * bws), c=cwheel[j], linestyle=lswheel[i], label=f"ND Flux (parent={ppdg}), {offaxis_label}")

plt.legend(bbox_to_anchor=(1, 0, 0.75, 1))
plt.yscale('log')
plt.ylabel(r"flux /cm${}^2$ /GeV /POT")
plt.xlabel(r"$E_\nu$ (GeV)")
plt.savefig("flux_ND_byparent.pdf", bbox_inches='tight')