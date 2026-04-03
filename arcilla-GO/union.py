import os
from ase.io import read, write
import numpy as np

# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Use paths relative to the script directory
atoms = read(os.path.join(script_dir, 'dump.unitcell'), format='lammps-data', atom_style='full')

# ... rest of the code ...

write(os.path.join(script_dir, 'POSCAR_arcilla'), atoms, format='vasp')