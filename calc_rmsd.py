import csv
import shutil
import subprocess
from pathlib import Path

from rdkit import Chem
from spyrmsd import molecule, rmsd

COMPLEXES = ["181l", "182l", "183l", "184l", "185l"]
DOCKING_DIR = Path("docking")
RMSD_THRESHOLD = 2.0


def first_pose(pdbqt: Path) -> str:
    lines = []
    for line in pdbqt.read_text().splitlines(keepends=True):
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("MODEL"):
            lines.append(line)
    return "".join(lines)


def read_sdf(path: Path) -> Chem.Mol:
    mol = Chem.MolFromMolFile(str(path), removeHs=True)
    if mol is None:
        raise RuntimeError(f"RDKit не смог прочитать {path}")
    return mol


def read_docked_pose(pdbqt: Path) -> Chem.Mol:
    result = subprocess.run(
        ["obabel", "-ipdbqt", "-osdf"],
        input=first_pose(pdbqt),
        capture_output=True,
        text=True,
    )
    mol = Chem.MolFromMolBlock(result.stdout, removeHs=True)
    if result.returncode != 0 or mol is None:
        raise RuntimeError(f"не удалось преобразовать {pdbqt}: {result.stderr.strip()}")
    return mol


def symmetric_rmsd(pdb_dir: Path) -> float:
    native = molecule.Molecule.from_rdkit(read_sdf(pdb_dir / "native_ligand.sdf"))
    docked = molecule.Molecule.from_rdkit(read_docked_pose(pdb_dir / "output.pdbqt"))

    if native.natoms != docked.natoms:
        raise RuntimeError(
            f"разное число тяжелых атомов: {native.natoms} и {docked.natoms}"
        )

    return float(rmsd.symmrmsd(
        native.coordinates,
        docked.coordinates,
        native.atomicnums,
        docked.atomicnums,
        native.adjacency_matrix,
        docked.adjacency_matrix,
    ))


def main() -> None:
    if shutil.which("obabel") is None:
        raise SystemExit("Не найдена программа obabel")

    results = []
    for pdb in COMPLEXES:
        try:
            value = symmetric_rmsd(DOCKING_DIR / pdb)
        except (OSError, RuntimeError) as error:
            print(f"{pdb}: {error}")
            continue
        success = value < RMSD_THRESHOLD
        results.append((pdb, value, success))
        print(f"{pdb}: RMSD {value:.2f} Å, {'успех' if success else 'неуспех'}")

    if not results:
        raise SystemExit("Нет результатов")

    successes = sum(success for _, _, success in results)
    share = successes / len(results) * 100
    print(f"\nУспешных поз (RMSD < {RMSD_THRESHOLD} Å): {successes} из {len(results)}, {share:.0f}%")

    table = DOCKING_DIR / "results.csv"
    with table.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Complex", "RMSD_A", "Success"])
        for pdb, value, success in results:
            writer.writerow([pdb, f"{value:.2f}", "Yes" if success else "No"])
    print(f"Таблица сохранена в {table}")


if __name__ == "__main__":
    main()
