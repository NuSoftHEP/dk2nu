# pydk2nu

## Building

### Requirements

Python3, CMake >3.14

Add `-DDK2NU_PYTHON_ENABLED=ON` to your cmake invocation for the main dk2nu 
project to enable thepython bindings.

## Usage

The main envisioned use case for these bindings is calculating neutrino rays 
and flux predictions from input dk2nu files. We do not actively support 
workflows looking to modify and re-persist dk2nu files.

The main dk2nu classes have python wrappers that expose their data members to 
python scripts.

A helper class is provided for reading lists of dk2nu files and iterating over
entries, rather than relying on pyroot or uproot to get instances out of trees 
on the python side. This makes the bindings a little less flexible, but a lot 
simpler and more self contained.

To read the total POT and vector of decay products from 2 input dk2nu files, 
you might use something like:

```python
import pydk2nu as dk

rdr = dk.dk2nuFileReader(["./dk2nu.file1.root", "./dk2nu.file2.root"])

POT = sum([m.pots for m in rdr.metas()])

for d in dk.decays():
  # do something with each decay

```

A common task is calculating the weight for all decays in an input chain
to send a neutrino through a small flux window about a point. The below
snippet gives an example of how to do so:

```python
import sys
import numpy as np
import pydk2nu as dk

rdr = dk.dk2nuFileReader(sys.argv[1:])

POT = sum([m.pots for m in rdr.metas()])

point_cm = [0, 0, 57400]

nus = np.array([(dk.decay.ntype, *dk.decay_through_point(point_cm)) for dk in fr.decays()])
# select only the muon neutrinos
numus = nus[nus[:,0] == 14]
```

In fact, this task is so common, we provide an interface to perform the
tight loop in C++ and only return the resulting numpy array.

```python
import sys
import numpy as np
import pydk2nu as dk

rdr = dk.dk2nuFileReader(sys.argv[1:])

POT = sum([m.pots for m in rdr.metas()])

point_cm = [0, 0, 57400]
nus = rdr.decay_all_through_point(point_cm)
# select only the muon neutrinos
numus = nus[nus[:,0] == 14]
```

This approach has been benchmarked to be approximately 30% faster on identical
inputs, so if the greater degree of flexibility afforded by the python looping
approach is useful, you don't sacrifice much performance for doing so. The vast
majority of the computation time is spent on ROOT I/O, and then on the 
python/C++ interface, and then neglibly on the weight and energy calculations.

When the flux through multiple points is needed then it is significantly more
performant to only loop over the dk2nu entries once. We provide a single 
interface for calculating the neutrino ray energies and weights through an 
arbitrary number of flux windows.

```python
import sys
import numpy as np
import pydk2nu as dk

rdr = dk.dk2nuFileReader(sys.argv[1:])

POT = sum([m.pots for m in rdr.metas()])

ND_Z_cm = 57400
offaxis_positions_m = [0, 5, 10, 20, 25]

nus = fr.decay_all_through_points([ [x*100, 0, ND_Z_cm] for x in offaxis_positions_m ])
# select only the muon neutrinos
numus = nus[nus[:,0] == 14]

window1_energies = numus[:,1]
window1_weights = numus[:,2]

window2_energies = numus[:,3]
window2_weights = numus[:,4]
```

Two complete example scripts for making Near Detector flux predictions can 
be found below:
* [example_off_axis.py](./example_off_axis.py)
* [example_decay_parent.py](./example_decay_parent.py)
