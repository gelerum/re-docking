import csv
from pathlib import Path

from rdkit import Chem

COMPLEXES = ["181l", "182l", "183l", "184l", "185l"]
DOCKING_DIR = Path("docking")
INPUTS_DIR = Path("inputs")

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def ligand_smiles(sdf: Path) -> str:
    mol = Chem.MolFromMolFile(str(sdf))
    if mol is None:
        raise RuntimeError(f"RDKit не смог прочитать {sdf}")
    return Chem.MolToSmiles(mol)


def chain_sequences(pdb: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            chain = line[21]
            sequences[chain] = sequences.get(chain, "") + THREE_TO_ONE.get(line[17:20], "X")
    if not sequences:
        raise RuntimeError(f"в {pdb} нет атомов CA")
    return sequences


def boltz_yaml(sequences: dict[str, str], smiles: str) -> str:
    lines = ["version: 1", "sequences:"]
    for chain, sequence in sequences.items():
        lines += ["  - protein:", f"      id: {chain}", f"      sequence: {sequence}"]
    lines += ["  - ligand:", "      id: L", f"      smiles: '{smiles}'"]
    return "\n".join(lines) + "\n"


def main() -> None:
    boltz_dir = INPUTS_DIR / "boltz"
    boltz_dir.mkdir(parents=True, exist_ok=True)

    diffdock_rows = []
    for pdb in COMPLEXES:
        pdb_dir = DOCKING_DIR / pdb
        smiles = ligand_smiles(pdb_dir / "native_ligand.sdf")
        sequences = chain_sequences(pdb_dir / "receptor.pdb")

        diffdock_rows.append([pdb, (pdb_dir / "receptor.pdb").resolve(), smiles, ""])
        (boltz_dir / f"{pdb}.yaml").write_text(boltz_yaml(sequences, smiles))

        lengths = ", ".join(f"{chain}: {len(seq)}" for chain, seq in sequences.items())
        print(f"{pdb}: {smiles}, цепи {lengths}")

    with (INPUTS_DIR / "diffdock.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["complex_name", "protein_path", "ligand_description", "protein_sequence"])
        writer.writerows(diffdock_rows)

    print(f"Входные файлы сохранены в {INPUTS_DIR}")


if __name__ == "__main__":
    main()
