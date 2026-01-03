from rdkit import Chem  # type: ignore
from rdkit.Chem import AllChem  # type: ignore

def element_color(el):
    colors = {
        'H': (1, 1, 1, 1),
        'C': (0.2, 0.2, 0.2, 1),
        'O': (1, 0, 0, 1),
        'N': (0, 0, 1, 1),
        'S': (1, 1, 0, 1),
        'Cl': (0, 1, 0, 1),
        'F': (0.5, 1, 0.5, 1)
    }
    return colors.get(el, (0.5, 0.5, 0.5, 1))

def generate_blender_py(xyz_path, name):
    import os

    with open(xyz_path, 'r') as f:
        lines = f.readlines()[2:]

    atoms = [line.strip().split() for line in lines]
    coordinates = [(float(x), float(y), float(z)) for (_, x, y, z) in atoms]
    elements = [el for (el, _, _, _) in atoms]

    # Known SMILES dict
    smiles_dict = {
        "methane": "C", "ethane": "CC", "propane": "CCC", "butane": "CCCC",
        "methanol": "CO", "cholesterol": "C[C@H](CCC(=O)O)C1CCC2C3CCC4=CC(=O)CCC4(C)C3CCC12C",
            "methane": "C",
    "ethane": "CC",
    "propane": "CCC",
    "butane": "CCCC",
    "pentane": "CCCCC",
    "hexane": "CCCCCC",
    "ethene": "C=C",
    "ethyne": "C#C",
    "benzene": "c1ccccc1",
    "water":"H2O", 
    # Alcohols
    "methanol": "CO",
    "ethanol": "CCO",
    "propanol": "CCCO",
    "butanol": "CCCCO",
    "isopropanol": "CC(O)C",
    "glycerol": "C(C(CO)O)O",

    # Acids
    "formic acid": "C(=O)O",
    "acetic acid": "CC(=O)O",
    "propanoic acid": "CCC(=O)O",
    "butanoic acid": "CCCC(=O)O",
    "benzoic acid": "c1ccc(cc1)C(=O)O",

    # Amines
    "methylamine": "CN",
    "ethylamine": "CCN",
    "aniline": "c1ccccc1N",

    # Aldehydes and Ketones
    "formaldehyde": "C=O",
    "acetaldehyde": "CC=O",
    "acetone": "CC(=O)C",

    # Carbohydrates
    "glucose": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    "fructose": "C(C1C(C(C(C(O1)(CO))O)O)O)=O",
    "sucrose": "OCC1OC(O)C(O)C(O)C1OC2(C(C(C(C(O2)CO)O)O)O)",

    # Amino acids
    "glycine": "NCC(=O)O",
    "alanine": "CC(C(=O)O)N",
    "valine": "CC(C)C(C(=O)O)N",
    "leucine": "CC(C)CC(C(=O)O)N",
    "phenylalanine": "C1=CC=C(C=C1)CC(C(=O)O)N",
    "tryptophan": "C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N",

    # Inorganic molecules
    "ammonia": "N",
    "carbon dioxide": "O=C=O",
    "nitric acid": "O=N(=O)O",
    "sulfuric acid": "OS(=O)(=O)O",
    "phosphoric acid": "OP(=O)(O)O",
    "hydrogen peroxide": "OO",
    "ozone": "O=O[O]",

    # Nucleotides / Bases
    "adenine": "C1=NC2=C(N1)N=CN2",
    "guanine": "C1=NC2=C(N1)C(=O)N=CN2",
    "cytosine": "C1=CN=CN1",
    "uracil": "C1=CC(=O)NC(=O)N1",
    "thymine": "CC1=CN(C(=O)NC1=O)C",

    # Other biological
    "cholesterol": "C[C@H](CCC(=O)O)C1CCC2C3CCC4=CC(=O)CCC4(C)C3CCC12C",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "nicotine": "CN1CCCC1C2=CN=CC=C2",

    # Gases
    "oxygen": "O=O",
    "hydrogen": "[H][H]",
    "nitrogen": "N#N",
    "chlorine": "ClCl",
    "fluorine": "F[F]",

    # Extra important organics
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "paracetamol": "CC(=O)NC1=CC=C(C=C1)O",
    "acetylsalicylic acid": "CC(=O)Oc1ccccc1C(=O)O"
        # add more as needed...
    }

    smiles = smiles_dict[name]
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

    # Get bond info
    bonds = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        order = bond.GetBondTypeAsDouble()
        bonds.append((i, j, order))

    conf = mol.GetConformer()
    rdkit_coords = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]

    py_path = f"molecules/{name}_blender.py"
    os.makedirs("molecules", exist_ok=True)

    with open(py_path, 'w') as f:
        f.write("import bpy\nimport mathutils\nimport math\n\n")
        f.write("bpy.ops.object.select_all(action='SELECT')\n")
        f.write("bpy.ops.object.delete(use_global=False)\n\n")

        f.write("# === ATOMS ===\n")
        for i, (el, (x, y, z)) in enumerate(zip(elements, coordinates)):
            color = element_color(el)
            offset_x = 0

            f.write(f"# Atom {i}: {el}\n")
            f.write(f"bpy.ops.mesh.primitive_uv_sphere_add(radius=0.3, location=({x + offset_x}, {y}, {z}))\n")
            f.write(f"atom = bpy.context.object\n")
            f.write(f"mat = bpy.data.materials.new(name='Mat_{i}')\n")
            f.write(f"mat.diffuse_color = {color}\n")
            f.write("atom.data.materials.append(mat)\n")

            f.write(f"atom.keyframe_insert(data_path='location', frame=1)\n")
            f.write(f"atom.location = ({x}, {y}, {z})\n")
            f.write(f"atom.keyframe_insert(data_path='location', frame=50)\n\n")

            f.write(f"bpy.ops.object.text_add(location=({x + 0.4}, {y}, {z}))\n")
            f.write("text_obj = bpy.context.object\n")
            f.write(f"text_obj.data.body = '{el}'\n")
            f.write(f"text_obj.scale = (0.3, 0.3, 0.3)\n\n")

        f.write("# === BONDS ===\n")
        for i, (start_idx, end_idx, order) in enumerate(bonds):
            start = coordinates[start_idx]
            end = coordinates[end_idx]
            dx, dy, dz = end[0] - start[0], end[1] - start[1], end[2] - start[2]
            bond_count = int(order)
            shift = 0.15

            for b in range(bond_count):
                offset = (b - (bond_count - 1)/2) * shift
                f.write(f"# Bond {i} - {b+1}/{bond_count}\n")
                f.write(f"start = mathutils.Vector(({start[0]}, {start[1]}, {start[2]}))\n")
                f.write(f"end = mathutils.Vector(({end[0]}, {end[1]}, {end[2]}))\n")
                f.write("mid = (start + end) / 2\n")
                f.write("vec = end - start\n")
                f.write("length = vec.length\n")
                f.write("bpy.ops.mesh.primitive_cylinder_add(radius=0.07, depth=length, location=mid)\n")
                f.write("cyl = bpy.context.object\n")
                f.write("cyl.rotation_mode = 'QUATERNION'\n")
                f.write("cyl.rotation_quaternion = vec.to_track_quat('Z', 'Y')\n")
                if bond_count > 1:
                    f.write(f"cyl.location += vec.cross(mathutils.Vector((0, 0, 1))).normalized() * {offset}\n")
                f.write("\n")

        f.write(f"bpy.ops.export_scene.gltf(filepath='outputs/{name}.glb', export_animations=True)\n")

    print("✅ Blender script finished")
    print(f"✅ Checking for file: outputs/{name}.glb")

    return py_path
