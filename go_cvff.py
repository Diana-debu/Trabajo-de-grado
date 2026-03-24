import numpy as np
import math
import os
import makegraphitics as mg

# --- CONFIGURACIÓN ---
FILENAME_DATA = "go_cvff.data"
FILENAME_IN = "run_go_cvff.in"

# --- CLASES ---
class Atom:
    def __init__(self, id, type, x, y, z):
        self.id = id
        self.type = type # 1: C_sp2, 2: O, 3: H, 4: C_sp3 (oxidado)
        self.q = 0.0     
        self.x = x
        self.y = y
        self.z = z
        self.neighbors =[]

def calc_dist(a, b):
    return math.sqrt((a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2)

# =====================================================================
# PASO 1: GENERACIÓN CON MAKE-GRAPHITICS
# =====================================================================
print("1. Generando molécula de GO con make-graphitics...")

flake_radius = 10  # Radio del flake 
layout = [1, 1, 1] 

motif = mg.molecules.Hexagon_Graphene(flake_radius)
flake = mg.Crystal(motif, layout)

oxidiser = mg.reactors.Oxidiser(
    ratio=2.5, video_xyz=20, new_island_freq=1e14, method="rf"
)
flake = oxidiser.react(flake)
mg.Parameterise(flake)

# Guardar un XYZ temporal para extraer las coordenadas y elementos fácilmente
temp_name = "temp_go"
output = mg.Writer(flake, temp_name)
output.write_xyz(temp_name + ".xyz")

print("   Leyendo coordenadas generadas...")
atoms =[]
with open(temp_name + ".xyz", 'r') as f:
    lines = f.readlines()
    num_atoms = int(lines[0].strip())
    # lines[1] es el comentario del XYZ
    for i in range(2, 2 + num_atoms):
        parts = lines[i].split()
        elem = parts[0]
        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        
        # Mapear el símbolo químico a nuestro Tipo de Átomo CVFF
        if elem == 'C': atype = 1
        elif elem == 'O': atype = 2
        elif elem == 'H': atype = 3
        else: continue
        
        atoms.append(Atom(id=i-1, type=atype, x=x, y=y, z=z))

# Eliminar archivo temporal
#if os.path.exists(temp_name + ".xyz"):
#    os.remove(temp_name + ".xyz")
#if os.path.exists(temp_name + ".data"): # make-graphitics a veces saca data tmb
#    os.remove(temp_name + ".data")

print("   Reconstruyendo red de enlaces topológica por distancia...")
# Conectar átomos que estén a distancia de enlace (< 1.75 Å)
for i in range(len(atoms)):
    for j in range(i+1, len(atoms)):
        if calc_dist(atoms[i], atoms[j]) < 1.75:
            atoms[i].neighbors.append(atoms[j])
            atoms[j].neighbors.append(atoms[i])

# =====================================================================
# PASO 2: IDENTIFICACIÓN DEL TIPO 4 Y ASIGNACIÓN DE CARGAS (CVFF)
# =====================================================================
print("2. Identificando carbonos oxidados (Tipo 4) y balanceando cargas...")

CHARGE_O = -0.40
CHARGE_H = +0.20
CHARGE_C_OX = +0.20
CHARGE_C_SP2 = 0.00

for a in atoms:
    if a.type == 2: a.q = CHARGE_O
    elif a.type == 3: a.q = CHARGE_H
    elif a.type == 1:
        has_oxygen = any(n.type == 2 for n in a.neighbors)
        if has_oxygen:
            a.type = 4
            a.q = CHARGE_C_OX
        else:
            a.q = CHARGE_C_SP2

# Neutralidad exacta de carga
total_q = sum(a.q for a in atoms)
c_sp2 = [a for a in atoms if a.type == 1]
if c_sp2:
    correction = total_q / len(c_sp2)
    for a in c_sp2: a.q -= correction
print(f"   Carga neta del sistema ajustada a: {sum(a.q for a in atoms):.5f}")

# =====================================================================
# PASO 3: CONSTRUCCIÓN TOPOLÓGICA EXPLÍCITA
# =====================================================================
print("3. Generando listas topológicas (Bonds, Angles, Dihedrals, Impropers)...")

bond_set, bonds = set(),[]
for at in atoms:
    for n in at.neighbors:
        pair = tuple(sorted((at.id, n.id)))
        if pair not in bond_set:
            bond_set.add(pair)
            ts = set([at.type, n.type])
            if ts.issubset({1, 4}): btype = 1
            elif ts == {2, 4} or ts == {1, 2}: btype = 2
            elif ts == {2, 3}: btype = 3
            else: btype = 1
            bonds.append([len(bonds)+1, btype, pair[0], pair[1]])

angles =[]
for b in atoms:
    neis = b.neighbors
    if len(neis) < 2: continue
    for i in range(len(neis)):
        for j in range(i+1, len(neis)):
            a, c = neis[i], neis[j]
            center = b.type
            sides = set([a.type, c.type])
            
            if center in [1, 4]:
                if sides.issubset({1, 4}): atype = 1
                elif 2 in sides and sides.intersection({1, 4}): 
                    if any(n.id == c.id for n in a.neighbors): atype = 6 # Epóxido
                    else: atype = 2 # Hidroxilo
                elif sides == {2}: atype = 3
                else: atype = 1
            elif center == 2:
                if sides.issubset({1, 4}): atype = 4
                elif 3 in sides: atype = 5
                else: atype = 1
            else: atype = 1
            angles.append([len(angles)+1, atype, a.id, b.id, c.id])

dihedrals =[]
for b_data in bonds:
    b, c = next(x for x in atoms if x.id == b_data[2]), next(x for x in atoms if x.id == b_data[3])
    neis_b =[n for n in b.neighbors if n.id != c.id]
    neis_c =[n for n in c.neighbors if n.id != b.id]
    for n_b in neis_b:
        for n_c in neis_c:
            if n_b.id != n_c.id:
                dihedrals.append([len(dihedrals)+1, 1, n_b.id, b.id, c.id, n_c.id])

impropers =[]
for at in atoms:
    if at.type == 1 and len(at.neighbors) == 3:
        n = at.neighbors
        impropers.append([len(impropers)+1, 1, at.id, n[0].id, n[1].id, n[2].id])

# =====================================================================
# PASO 4: ESCRITURA DEL ARCHIVO DATA
# =====================================================================
print(f"4. Escribiendo archivo de datos: {FILENAME_DATA}")

# IMPORTANTE: Al ser un "flake", añadimos espacio vacío (padding) para evitar interacciones periódicas
PADDING = 15.0
xlo, xhi = min(a.x for a in atoms) - PADDING, max(a.x for a in atoms) + PADDING
ylo, yhi = min(a.y for a in atoms) - PADDING, max(a.y for a in atoms) + PADDING
zlo, zhi = -20.0, 20.0

with open(FILENAME_DATA, 'w') as f:
    f.write("Graphene Oxide Flake with CVFF Force Field\n\n")
    f.write(f"{len(atoms)} atoms\n{len(bonds)} bonds\n{len(angles)} angles\n")
    f.write(f"{len(dihedrals)} dihedrals\n{len(impropers)} impropers\n\n")
    f.write("4 atom types\n3 bond types\n6 angle types\n1 dihedral types\n1 improper types\n\n")
    f.write(f"{xlo:.4f} {xhi:.4f} xlo xhi\n{ylo:.4f} {yhi:.4f} ylo yhi\n{zlo:.4f} {zhi:.4f} zlo zhi\n\n")
    
    f.write("Masses\n\n1 12.011\n2 15.999\n3 1.008\n4 12.011\n\n")
    f.write("Atoms\n\n")
    for a in atoms: f.write(f"{a.id} 1 {a.type} {a.q:.4f} {a.x:.5f} {a.y:.5f} {a.z:.5f}\n")
    f.write("\nBonds\n\n")
    for b in bonds: f.write(f"{b[0]} {b[1]} {b[2]} {b[3]}\n")
    f.write("\nAngles\n\n")
    for a in angles: f.write(f"{a[0]} {a[1]} {a[2]} {a[3]} {a[4]}\n")
    f.write("\nDihedrals\n\n")
    for d in dihedrals: f.write(f"{d[0]} {d[1]} {d[2]} {d[3]} {d[4]} {d[5]}\n")
    f.write("\nImpropers\n\n")
    for i in impropers: f.write(f"{i[0]} {i[1]} {i[2]} {i[3]} {i[4]} {i[5]}\n")

# =====================================================================
# PASO 5: ESCRITURA DEL SCRIPT LAMMPS (CVFF)
# =====================================================================
print(f"5. Escribiendo input para LAMMPS: {FILENAME_IN}")
lammps_script = """# Graphene Oxide FLAKE Relaxation using CVFF

##---------------INITIALIZATION-------------------------------
units          real
dimension      3
boundary       p p p
atom_style     full

##---------------FORCE FIELDS (CVFF)--------------------------
bond_style      harmonic
angle_style     harmonic
dihedral_style  charmm
improper_style  cvff

pair_style      lj/cut/coul/long 10.0
pair_modify     mix arithmetic
kspace_style    pppm 1.0e-4
special_bonds   lj/coul 0.0 0.0 1.0 

##---------------READ DATA------------------------------------
read_data       go_cvff.data

##---------------COEFFICIENTS (CVFF)--------------------------
pair_coeff 1 1 0.0700 3.5500 # C sp2
pair_coeff 2 2 0.1500 3.0000 # O
pair_coeff 3 3 0.0250 1.0000 # H 
pair_coeff 4 4 0.0700 3.5500 # C sp3

bond_coeff 1 350.0 1.400  # C-C
bond_coeff 2 320.0 1.420  # C-O
bond_coeff 3 450.0 0.960  # O-H

angle_coeff 1 40.0 120.0  # C-C-C
angle_coeff 2 50.0 109.5  # C-C-O
angle_coeff 3 40.0 109.5  # O-C-O
angle_coeff 4 60.0 60.0   # C-O-C 
angle_coeff 5 35.0 108.0  # C-O-H
angle_coeff 6 50.0 60.0   # C-C-O (Epoxido anillo)

dihedral_coeff 1 2.5 2 180 0.0
improper_coeff 1 2.0 -1 2  

##---------------MINIMIZATION---------------------------------
thermo          10
min_modify      dmax 0.05
minimize        1.0e-6 1.0e-8 5000 10000

##---------------DYNAMICS-------------------------------------
reset_timestep  0
timestep        0.5

dump            1 all custom 1000 go_flake_relaxed.lammpstrj id type q x y z

# IMPORTANTE: Como es un flake (no infinito), relajamos en NVT (Volumen Constante)
# Si usamos NPT, la caja intentará encogerse sobre el vacío y colapsará.
fix             1 all nvt temp 300 300 50
thermo          100
thermo_style    custom step temp pe ke press
run             50000

print "Simulación del Flake Completada con CVFF"
"""

with open(FILENAME_IN, 'w') as f:
    f.write(lammps_script)

print("¡Proceso Finalizado con Éxito!")