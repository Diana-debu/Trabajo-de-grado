import yaml
import numpy as np
import sys
# Cambia las líneas de "scripts" por estas:
from makegraphitics.read_lammpsdata import ReadLammpsData, Writer

sim = ReadLammpsData(sys.argv[1])


name = "out"
output = Writer(sim, name)
output.write_xyz(name + ".xyz")
output.write_lammps(name + ".data")
