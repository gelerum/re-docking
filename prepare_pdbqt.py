import shutil
import subprocess
from pathlib import Path

COMPLEXES = ["181l", "182l", "183l", "184l", "185l"]
DOCKING_DIR = Path("docking")

CONVERSIONS = [
    ("receptor.pdb", "receptor.pdbqt", ["-xr"]),
    ("native_ligand.sdf", "ligand.pdbqt", ["-h"]),
    ("native_ligand.sdf", "native_ligand.pdbqt", []),
]


def obabel(source: Path, target: Path, options: list[str]) -> None:
    result = subprocess.run(
        ["obabel", str(source), "-O", str(target), *options],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0 or not target.exists():
        raise RuntimeError(f"obabel не смог создать {target.name}: {result.stderr.strip()}")


def prepare(pdb_dir: Path) -> None:
    for source, target, options in CONVERSIONS:
        obabel(pdb_dir / source, pdb_dir / target, options)

    if "ROOT" not in (pdb_dir / "ligand.pdbqt").read_text():
        print(f"{pdb_dir.name}: в ligand.pdbqt нет ROOT, Vina может отработать некорректно")


def main() -> None:
    if shutil.which("obabel") is None:
        raise SystemExit("Не найдена программа obabel")

    failed = []
    for pdb in COMPLEXES:
        try:
            prepare(DOCKING_DIR / pdb)
        except (OSError, RuntimeError) as error:
            print(f"{pdb}: {error}")
            failed.append(pdb)
            continue
        print(f"{pdb}: PDBQT-файлы готовы")

    if failed:
        raise SystemExit(f"Не удалось подготовить: {', '.join(failed)}")


if __name__ == "__main__":
    main()
