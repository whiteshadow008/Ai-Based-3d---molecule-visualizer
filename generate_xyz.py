from rdkit import Chem  # type: ignore
from rdkit.Chem import AllChem  # type: ignore
import os

def generate_xyz(smiles: str, name: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol) != 0:
        raise RuntimeError("Embedding failed for molecule: " + name)

    if AllChem.UFFOptimizeMolecule(mol) != 0:
        print("⚠️ Warning: Optimization not fully converged.")

    conf = mol.GetConformer()
    atoms = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        atoms.append(f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}")

    os.makedirs("molecules", exist_ok=True)
    xyz_path = f"molecules/{name}.xyz"
    with open(xyz_path, "w") as f:
        f.write(f"{len(atoms)}\n{name}\n" + "\n".join(atoms))

    print(f"✅ XYZ written: {xyz_path}")
    return xyz_path
