import numpy as np
import random

# --- CONFIGURACIÓN ---
FILENAME_DATA = "graphene_oxide_full.data"
FILENAME_IN = "go_full.in.txt"

# Tamaño de la lámina
NX = 8
NY = 8
CC_DIST = 1.42

# Química
COVERAGE = 0.20  # 20% de oxidación
RATIO_EPOXY = 0.6 # 60% epóxidos, 40% hidroxilos

# --- PARÁMETROS DEL CAMPO DE FUERZA (Harmonic/Charmm) ---
# Estos valores son aproximaciones típicas para grafeno/GO en campos armónicos
MASSES = {1: 12.011, 2: 15.999, 3: 1.008} # C, O, H

# --- CLASES AUXILIARES ---
class Atom:
    def __init__(self, id, type, x, y, z):
        self.id = id
        self.type = type # 1:C, 2:O, 3:H
        self.x = x
        self.y = y
        self.z = z
        self.neighbors = [] # Lista de átomos conectados

# --- FUNCIONES DE TOPOLOGÍA ---
def get_bond_type(t1, t2):
    ts = sorted([t1, t2])
    if ts == [1, 1]: return 1 # C-C
    if ts == [1, 2]: return 2 # C-O
    if ts == [2, 3]: return 3 # O-H
    return 1

def get_angle_type(t1, t2, t3):
    # t2 es el átomo central
    center = t2
    sides = sorted([t1, t3])
    if center == 1: # C central
        if sides == [1, 1]: return 1 # C-C-C (120 deg)
        if sides == [1, 2]: return 2 # C-C-O 
        if sides == [2, 2]: return 3 # O-C-O
    if center == 2: # O central
        if sides == [1, 1]: return 4 # C-O-C (Epóxido, ~60 deg)
        if sides == [1, 3]: return 5 # C-O-H (Hidroxilo, ~109 deg)
    return 1

# --- GENERACIÓN DE ESTRUCTURA ---
atoms = []
bonds = []
angles = []
dihedrals = []
impropers = []

print("1. Generando red de carbonos...")
atom_id = 1
# Vectores base grafeno
a1 = np.array([3/2 * CC_DIST, np.sqrt(3)/2 * CC_DIST, 0])
a2 = np.array([3/2 * CC_DIST, -np.sqrt(3)/2 * CC_DIST, 0])

carbons = []
for i in range(NX):
    for j in range(NY):
        # Átomo A
        pos_a = i*a1 + j*a2
        a = Atom(atom_id, 1, pos_a[0], pos_a[1], 0.0)
        atoms.append(a); carbons.append(a); atom_id += 1
        
        # Átomo B (desplazado CC_DIST en X)
        pos_b = pos_a + np.array([CC_DIST, 0, 0])
        b = Atom(atom_id, 1, pos_b[0], pos_b[1], 0.0)
        atoms.append(b); carbons.append(b); atom_id += 1

# Ajustar caja
xs = [a.x for a in atoms]; ys = [a.y for a in atoms]
xlo, xhi = min(xs)-2, max(xs)+2
ylo, yhi = min(ys)-2, max(ys)+2
zlo, zhi = -10, 10

# Construir lista de vecinos inicial para carbonos (distancia < 1.55)
print("   Conectando carbonos...")
for i in range(len(carbons)):
    for j in range(i+1, len(carbons)):
        d2 = (carbons[i].x - carbons[j].x)**2 + (carbons[i].y - carbons[j].y)**2
        if d2 < 1.55**2:
            carbons[i].neighbors.append(carbons[j])
            carbons[j].neighbors.append(carbons[i])

print("2. Añadiendo grupos funcionales...")
num_sites = int(len(carbons) * COVERAGE)
target_idxs = random.sample(range(len(carbons)), num_sites)

for idx in target_idxs:
    c = carbons[idx]
    
    # Determinar tipo: Epóxido o Hidroxilo
    if random.random() < RATIO_EPOXY:
        # EPÓXIDO (Puente C-O-C)
        if len(c.neighbors) > 0:
            # Elegir un vecino C para hacer el puente
            c2 = c.neighbors[0] 
            
            # Posición O
            ox = (c.x + c2.x)/2
            oy = (c.y + c2.y)/2
            oz = 1.2 if random.random() > 0.5 else -1.2
            
            o_atom = Atom(atom_id, 2, ox, oy, oz)
            atom_id += 1
            atoms.append(o_atom)
            
            # Conectar topología (sin crear Bond object aun, solo grafo)
            c.neighbors.append(o_atom)
            o_atom.neighbors.append(c)
            c2.neighbors.append(o_atom)
            o_atom.neighbors.append(c2)
            
            # Deformar C ligeramente fuera del plano
            dir_z = 1 if oz > 0 else -1
            c.z += 0.3 * dir_z
            c2.z += 0.3 * dir_z
            
    else:
        # HIDROXILO (C-O-H)
        dir_z = 1 if random.random() > 0.5 else -1
        ox, oy, oz = c.x, c.y, c.z + (1.45 * dir_z)
        hx, hy, hz = ox + 0.8, oy, oz + (0.6 * dir_z)
        
        o_atom = Atom(atom_id, 2, ox, oy, oz)
        atom_id += 1
        h_atom = Atom(atom_id, 3, hx, hy, hz)
        atom_id += 1
        
        atoms.extend([o_atom, h_atom])
        
        # Conectar
        c.neighbors.append(o_atom)
        o_atom.neighbors.append(c)
        o_atom.neighbors.append(h_atom)
        h_atom.neighbors.append(o_atom)
        
        c.z += 0.2 * dir_z

print("3. Generando listas topológicas explícitas...")

# --- GENERAR BONDS ---
bond_set = set() # Para evitar duplicados 1-2 y 2-1
bond_list_final = []
b_id = 1

for at in atoms:
    for n in at.neighbors:
        pair = tuple(sorted((at.id, n.id)))
        if pair not in bond_set:
            bond_set.add(pair)
            btype = get_bond_type(at.type, n.type)
            bond_list_final.append([b_id, btype, pair[0], pair[1]])
            b_id += 1

# --- GENERAR ANGLES (A-B-C) ---
# Iterar sobre átomo central B
angle_list_final = []
ang_id = 1

for b in atoms:
    neis = b.neighbors
    if len(neis) < 2: continue
    
    # Combinaciones de 2 vecinos
    for i in range(len(neis)):
        for j in range(i+1, len(neis)):
            a = neis[i]
            c = neis[j]
            
            # Evitar triangulos cerrados pequeños en definición de ángulos si es necesario
            # Pero en MD clásico definimos todos.
            atype = get_angle_type(a.type, b.type, c.type)
            angle_list_final.append([ang_id, atype, a.id, b.id, c.id])
            ang_id += 1

# --- GENERAR DIHEDRALS (A-B-C-D) ---
# Iterar sobre enlace central B-C
dihedral_list_final = []
dih_id = 1

for bond in bond_list_final:
    # ids de bond: bond[2], bond[3]
    id_b, id_c = bond[2], bond[3]
    
    # Objetos atomos
    atom_b = next(x for x in atoms if x.id == id_b)
    atom_c = next(x for x in atoms if x.id == id_c)
    
    # Vecinos de B (excluyendo C)
    neis_a = [n for n in atom_b.neighbors if n.id != id_c]
    # Vecinos de C (excluyendo B)
    neis_d = [n for n in atom_c.neighbors if n.id != id_b]
    
    for atom_a in neis_a:
        for atom_d in neis_d:
            if atom_a.id == atom_d.id: continue # Es un anillo de 3, no diedro propio
            
            # Tipo simple basado en átomos centrales (para simplificar)
            dtype = 1 
            dihedral_list_final.append([dih_id, dtype, atom_a.id, id_b, id_c, atom_d.id])
            dih_id += 1

# --- GENERAR IMPROPERS (I-J-K-L) ---
# Se usan para mantener planaridad en carbonos sp2 (3 vecinos)
# I es el centro, J, K, L son vecinos
improper_list_final = []
imp_id = 1

for at in atoms:
    if at.type == 1 and len(at.neighbors) == 3: # Carbono con 3 vecinos
        n = at.neighbors
        # Definición típica LAMMPS: Centro, Vecino1, Vecino2, Vecino3
        improper_list_final.append([imp_id, 1, at.id, n[0].id, n[1].id, n[2].id])
        imp_id += 1

print(f"   Estadísticas: {len(atoms)} Atoms, {len(bond_list_final)} Bonds, {len(angle_list_final)} Angles, {len(dihedral_list_final)} Dihedrals, {len(improper_list_final)} Impropers")

# --- ESCRIBIR ARCHIVO DATA ---
print(f"4. Escribiendo {FILENAME_DATA}...")
with open(FILENAME_DATA, 'w') as f:
    f.write("# Graphene Oxide FULL Topology\n\n")
    
    f.write(f"{len(atoms)} atoms\n")
    f.write(f"{len(bond_list_final)} bonds\n")
    f.write(f"{len(angle_list_final)} angles\n")
    f.write(f"{len(dihedral_list_final)} dihedrals\n")
    f.write(f"{len(improper_list_final)} impropers\n\n")
    
    f.write("3 atom types\n3 bond types\n5 angle types\n1 dihedral types\n1 improper types\n\n")
    
    f.write(f"{xlo:.4f} {xhi:.4f} xlo xhi\n")
    f.write(f"{ylo:.4f} {yhi:.4f} ylo yhi\n")
    f.write(f"{zlo:.4f} {zhi:.4f} zlo zhi\n\n")
    
    f.write("Masses\n\n1 12.011\n2 15.999\n3 1.008\n\n")
    
    # Coeficientes (Ejemplo genérico, deben validarse para publicación)
    f.write("Pair Coeffs\n\n")
    f.write("1 0.070 3.55\n") # C
    f.write("2 0.152 3.15\n") # O
    f.write("3 0.000 0.00\n") # H (LJ cero a menudo si hay carga, o pequeño)
    
    f.write("\nBond Coeffs\n\n")
    f.write("1 469.0 1.40\n") # C-C
    f.write("2 350.0 1.43\n") # C-O
    f.write("3 450.0 0.96\n") # O-H
    
    f.write("\nAngle Coeffs\n\n")
    f.write("1 63.0 120.0\n")  # C-C-C
    f.write("2 50.0 109.5\n")  # C-C-O
    f.write("3 70.0 110.0\n")  # O-C-O
    f.write("4 60.0 61.0\n")   # C-O-C (Epóxido muy tenso)
    f.write("5 55.0 108.0\n")  # C-O-H
    
    f.write("\nDihedral Coeffs\n\n")
    # K, n, d, w (Charmm style: E = K * (1 + cos(n*phi - d)))
    f.write("1 2.5 2 180 0.0\n") # Genérico
    
    f.write("\nImproper Coeffs\n\n")
    # K, d, n (Cvff style: E = K * (1 + d * cos(n * phi))) 
    # O Harmonic: K * (chi - chi0)^2. LAMMPS 'cvff' usa K(1+d cos(n phi))
    f.write("1 2.2 -1 2\n") # Mantiene planaridad
    
    f.write("\nAtoms\n\n")
    for a in atoms:
        f.write(f"{a.id} 1 {a.type} 0.0 {a.x:.5f} {a.y:.5f} {a.z:.5f}\n")
        
    f.write("\nBonds\n\n")
    for b in bond_list_final:
        f.write(f"{b[0]} {b[1]} {b[2]} {b[3]}\n")

    f.write("\nAngles\n\n")
    for a in angle_list_final:
        f.write(f"{a[0]} {a[1]} {a[2]} {a[3]} {a[4]}\n")
        
    f.write("\nDihedrals\n\n")
    for d in dihedral_list_final:
        f.write(f"{d[0]} {d[1]} {d[2]} {d[3]} {d[4]} {d[5]}\n")
        
    f.write("\nImpropers\n\n")
    for i in improper_list_final:
        f.write(f"{i[0]} {i[1]} {i[2]} {i[3]} {i[4]} {i[5]}\n")

# --- ESCRIBIR INPUT LAMMPS ---
print(f"5. Escribiendo {FILENAME_IN}...")
lammps_script = """# Uniaxial Tensile Test of Graphene Oxide (Full Topology)

##---------------INITIALIZATION-------------------------------
units          real
dimension      3
boundary       p p p
atom_style     full
newton         on

##---------------FORCE FIELDS---------------------------------
bond_style      harmonic
angle_style     harmonic
dihedral_style  charmm
improper_style  cvff
pair_style      lj/cut 10.0

##---------------ATOM DEFINITION------------------------------
read_data       graphene_oxide_full.data

##---------------SETTINGS-------------------------------------
timestep        0.5
variable        ts equal 0.5

##---------------COMPUTES-------------------------------------
compute         1 all stress/atom NULL
compute         2 all reduce sum c_1[1] c_1[2]

variable        Lx equal lx
variable        Ly equal ly
variable        Lz equal lz
variable        Vol equal vol
variable        thickn equal 3.4

##---------------RELAXATION--------------------------------------
# Crucial para GO generado proceduralmente: minimizar tensiones de epoxidos
minimize 1.0e-4 1.0e-6 1000 10000

reset_timestep  0
fix             1 all npt temp 300 300 50 x 1 1 1000 y 1 1 1000
thermo          1000
run             10000

##---------------DEFORMATION--------------------------------------
unfix           1
reset_timestep  0

fix             1 all npt temp 300 300 50 x 1 1 1000 z 1 1 1000
fix             2 all ave/time 1 100 100 c_2[1] c_2[2]
fix             3 all ave/time 1 100 100 v_Lx v_Ly v_Lz v_Vol

variable        srate equal 1.0e9
variable        srate1 equal "v_srate / 1.0e15"

fix             4 all deform 1 y erate ${srate1} units box remap x

##---------------THERMO-OUTPUTS--------------------------------------
variable        CorVol equal f_3[4]*v_thickn/(f_3[3])
variable        ConvoFac equal 1/1.01325e4
variable        sigmaxx equal f_2[1]*v_ConvoFac/v_CorVol
variable        sigmayy equal f_2[2]*v_ConvoFac/v_CorVol
variable        StrainPerTs equal v_srate1*v_ts
variable        strain equal v_StrainPerTs*step

thermo          100
thermo_style    custom step temp v_strain v_sigmaxx v_sigmayy pe ke lx ly vol

# Dump con tensiones para ver en OVITO
dump            1 all custom 1000 tensile_go_full.lammpstrj id type x y z c_1[1] c_1[2]

run             100000
"""

with open(FILENAME_IN, 'w') as f:
    f.write(lammps_script)

print("¡Proceso completado con topología completa!")