from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdDetermineBonds
from rdkit.Geometry import Point3D

COMPLEXES = ["181l", "182l", "183l", "184l", "185l"]
DOCKING_DIR = Path("docking")
DIFFDOCK_DIR = Path("predictions/diffdock")
BOLTZ_DIR = Path("predictions/boltz")
POSES_DIR = Path("poses")

BOLTZ_LIGAND_CHAIN = "L"

AUTODOCK_ELEMENTS = {"A": "C", "OA": "O", "NA": "N", "NS": "N", "SA": "S", "HD": "H", "HS": "H"}


def native_template(pdb: str) -> Chem.Mol:
    return Chem.MolFromMolFile(str(DOCKING_DIR / pdb / "native_ligand.sdf"))


def build_ligand(atoms: list[tuple[str, np.ndarray]], template: Chem.Mol) -> Chem.Mol:
    heavy = [(element, xyz) for element, xyz in atoms if element != "H"]
    mol = Chem.RWMol()
    conformer = Chem.Conformer(len(heavy))
    for index, (element, xyz) in enumerate(heavy):
        mol.AddAtom(Chem.Atom(element))
        conformer.SetAtomPosition(index, Point3D(*map(float, xyz)))
    mol.AddConformer(conformer)
    rdDetermineBonds.DetermineConnectivity(mol)
    ligand = AllChem.AssignBondOrdersFromTemplate(template, mol.GetMol())
    for atom in ligand.GetAtoms():
        atom.SetNoImplicit(False)
        atom.SetNumRadicalElectrons(0)
    Chem.SanitizeMol(ligand)
    return ligand


def save_pose(mol: Chem.Mol, method: str, pdb: str) -> None:
    mol = Chem.AddHs(Chem.RemoveHs(mol), addCoords=True)
    Chem.MolToMolFile(mol, str(POSES_DIR / method / f"{pdb}.sdf"))


def pdb_coords(line: str) -> np.ndarray:
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])


def vina_pose(pdb: str) -> None:
    atoms = []
    for line in (DOCKING_DIR / pdb / "output.pdbqt").read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM", "HETATM")):
            atom_type = line[77:79].strip()
            element = AUTODOCK_ELEMENTS.get(atom_type, atom_type.capitalize())
            atoms.append((element, pdb_coords(line)))
    save_pose(build_ligand(atoms, native_template(pdb)), "vina", pdb)


def diffdock_pose(pdb: str) -> None:
    path = DIFFDOCK_DIR / pdb / "rank1.sdf"
    mol = Chem.MolFromMolFile(str(path))
    if mol is None:
        raise RuntimeError(f"RDKit не смог прочитать {path}")
    save_pose(mol, "diffdock", pdb)


def kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mobile_center = mobile.mean(axis=0)
    target_center = target.mean(axis=0)
    u, _, vt = np.linalg.svd((mobile - mobile_center).T @ (target - target_center))
    sign = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ np.diag([1, 1, sign]) @ u.T
    return rotation, target_center - rotation @ mobile_center


def ca_coords(lines: list[str]) -> np.ndarray:
    return np.array([pdb_coords(line) for line in lines if line.startswith("ATOM") and line[12:16].strip() == "CA"])


def boltz_pose(pdb: str) -> None:
    path = next(BOLTZ_DIR.rglob(f"{pdb}_model_0.pdb"), None)
    if path is None:
        raise RuntimeError(f"не найден {pdb}_model_0.pdb в {BOLTZ_DIR}")

    lines = path.read_text().splitlines()
    protein = [line for line in lines if line.startswith("ATOM") and line[21] != BOLTZ_LIGAND_CHAIN]
    ligand = [line for line in lines if line.startswith(("ATOM", "HETATM")) and line[21] == BOLTZ_LIGAND_CHAIN]

    predicted_ca = ca_coords(protein)
    crystal_ca = ca_coords((DOCKING_DIR / pdb / "receptor.pdb").read_text().splitlines())
    if len(predicted_ca) != len(crystal_ca):
        raise RuntimeError(f"разное число CA: {len(predicted_ca)} и {len(crystal_ca)}")

    rotation, shift = kabsch(predicted_ca, crystal_ca)
    aligned = predicted_ca @ rotation.T + shift
    ca_rmsd = np.sqrt(((aligned - crystal_ca) ** 2).sum(axis=1).mean())

    def move(line: str) -> np.ndarray:
        return rotation @ pdb_coords(line) + shift

    protein_lines = [line[:30] + "{:8.3f}{:8.3f}{:8.3f}".format(*move(line)) + line[54:] for line in protein]
    (POSES_DIR / "boltz" / f"{pdb}_protein.pdb").write_text("\n".join(protein_lines) + "\nEND\n")

    atoms = [(line[76:78].strip().capitalize() or line[12:16].strip()[0], move(line)) for line in ligand]
    save_pose(build_ligand(atoms, native_template(pdb)), "boltz", pdb)
    print(f"    CA RMSD белка после выравнивания: {ca_rmsd:.2f} Å")


METHODS = {"vina": vina_pose, "diffdock": diffdock_pose, "boltz": boltz_pose}


def main() -> None:
    RDLogger.DisableLog("rdApp.warning")
    for method, collect in METHODS.items():
        (POSES_DIR / method).mkdir(parents=True, exist_ok=True)
        for pdb in COMPLEXES:
            print(f"{method} {pdb}")
            try:
                collect(pdb)
            except (OSError, RuntimeError, ValueError) as error:
                print(f"    пропущен: {error}")


if __name__ == "__main__":
    main()
