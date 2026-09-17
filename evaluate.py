import csv
from pathlib import Path

import numpy as np
from posebusters import PoseBusters
from rdkit import Chem
from spyrmsd import molecule, rmsd

COMPLEXES = ["181l", "182l", "183l", "184l", "185l"]
METHODS = {"vina": "AutoDock Vina", "diffdock": "DiffDock", "boltz": "Boltz-1"}
DOCKING_DIR = Path("docking")
POSES_DIR = Path("poses")
RESULTS_DIR = Path("results")
RMSD_THRESHOLD = 2.0


def heavy_atom_rmsd(pose: Path, native: Path) -> float:
    predicted = molecule.Molecule.from_rdkit(Chem.MolFromMolFile(str(pose)))
    reference = molecule.Molecule.from_rdkit(Chem.MolFromMolFile(str(native)))
    return float(rmsd.symmrmsd(
        reference.coordinates,
        predicted.coordinates,
        reference.atomicnums,
        predicted.atomicnums,
        reference.adjacency_matrix,
        predicted.adjacency_matrix,
    ))


def failed_checks(buster: PoseBusters, pose: Path, native: Path, protein: Path) -> list[str]:
    report = buster.bust(mol_pred=pose, mol_true=native, mol_cond=protein).iloc[0]
    return [
        name
        for name, passed in report.items()
        if isinstance(passed, (bool, np.bool_)) and "rmsd" not in name.lower() and not passed
    ]


def evaluate(buster: PoseBusters, method: str, pdb: str) -> dict:
    pose = POSES_DIR / method / f"{pdb}.sdf"
    native = DOCKING_DIR / pdb / "native_ligand.sdf"
    protein = POSES_DIR / method / f"{pdb}_protein.pdb"
    if not protein.exists():
        protein = DOCKING_DIR / pdb / "receptor.pdb"

    if not pose.exists():
        return {"rmsd": None, "failed": ["нет позы"]}
    return {
        "rmsd": heavy_atom_rmsd(pose, native),
        "failed": failed_checks(buster, pose, native, protein),
    }


def percent(count: int, total: int) -> str:
    return f"{count}/{total} ({count / total * 100:.0f}%)"


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    buster = PoseBusters(config="redock")

    details = []
    summary = ["| Метод | RMSD < 2 Å | PB-valid | RMSD < 2 Å и PB-valid |", "|---|---|---|---|"]
    for method, title in METHODS.items():
        accurate = valid = both = 0
        for pdb in COMPLEXES:
            result = evaluate(buster, method, pdb)
            is_accurate = result["rmsd"] is not None and result["rmsd"] < RMSD_THRESHOLD
            is_valid = not result["failed"]
            accurate += is_accurate
            valid += is_valid
            both += is_accurate and is_valid

            rmsd_text = "—" if result["rmsd"] is None else f"{result['rmsd']:.2f}"
            details.append([title, pdb, rmsd_text, "Yes" if is_valid else "No", ", ".join(result["failed"])])
            print(f"{title} {pdb}: RMSD {rmsd_text}, не пройдено: {', '.join(result['failed']) or 'ничего'}")

        total = len(COMPLEXES)
        summary.append(f"| {title} | {percent(accurate, total)} | {percent(valid, total)} | {percent(both, total)} |")

    with (RESULTS_DIR / "details.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Method", "Complex", "RMSD_A", "PB_valid", "Failed_checks"])
        writer.writerows(details)

    table = "\n".join(summary) + "\n"
    (RESULTS_DIR / "summary.md").write_text(table)
    print("\n" + table)


if __name__ == "__main__":
    main()
