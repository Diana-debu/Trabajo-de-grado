import os
from ase.io import read, write

# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

atoms = read(os.path.join(script_dir, 'go_cvff.data'), format='lammps-data', atom_style='full')

type_to_symbol = {
    1: 'C',
    2: 'O',
    3: 'H',
    4: 'C'
}

symbols = [type_to_symbol[t] for t in atoms.arrays['type']]
atoms.set_chemical_symbols(symbols)

write(os.path.join(script_dir, 'POSCAR_OG'), atoms, format='vasp')
print("POSCAR_OG generado con", len(atoms), "átomos")