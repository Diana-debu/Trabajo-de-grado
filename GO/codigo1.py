import makegraphitics as mg

# 1. Definir la base (Rectangular Graphene)
# Dimensiones en Angstroms (50x50)
motif = mg.molecules.Rectangle_Graphene(50, 50)

# 2. Crear el cristal (lámina única)
# [1, 1, 1] indica que no se repite la unidad en ninguna dirección
flake = mg.Crystal(motif, [1, 1, 1])

# 3. Configurar el Oxidador
# ratio=2.5 es la relación C/O (un valor común para GO)
# method="rf" usa Random Forest para una distribución realista de grupos funcionales
oxidiser = mg.reactors.Oxidiser(
    ratio=2.5, 
    video_xyz=20, 
    new_island_freq=1e14, 
    method="rf"
)

# 4. Ejecutar la oxidación
print("Oxidando la estructura... esto puede tardar unos segundos.")
flake = oxidiser.react(flake)

# 5. Parametrizar (Asigna tipos de átomos para campos de fuerza como OPLS)
mg.Parameterise(flake, flake.vdw_defs)

# 6. Guardar archivos
name = "GO_structure"
output = mg.Writer(flake, name)

# Escribir XYZ (Ideal para visualización rápida en OVITO)
output.write_xyz(name + ".xyz")

# Escribir LAMMPS data (Ideal si quieres ver enlaces y tipos de átomos en OVITO)
output.write_lammps(name + ".data")

print(f"Archivos generados: {name}.xyz y {name}.data")