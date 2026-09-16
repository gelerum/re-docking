import shutil
import subprocess
from pathlib import Path

COMPLEXES = ["181l", "182l", "183l", "184l", "185l"]
DOCKING_DIR = Path("docking")

PADDING = 4.0
EXHAUSTIVENESS = 8
NUM_MODES = 9
ENERGY_RANGE = 3


def ligand_box(pdbqt: Path) -> tuple[list[float], list[float]]:
    coords = [
        (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        for line in pdbqt.read_text().splitlines()
        if line.startswith(("ATOM", "HETATM"))
    ]
    if not coords:
        raise RuntimeError(f"в {pdbqt.name} нет атомов")

    center, size = [], []
    for axis in zip(*coords):
        low, high = min(axis), max(axis)
        center.append((low + high) / 2)
        size.append(high - low + 2 * PADDING)
    return center, size


def write_config(pdb_dir: Path) -> Path:
    center, size = ligand_box(pdb_dir / "native_ligand.pdbqt")
    lines = [
        "receptor = receptor.pdbqt",
        "ligand = ligand.pdbqt",
        "out = output.pdbqt",
        *(f"center_{axis} = {value:.3f}" for axis, value in zip("xyz", center)),
        *(f"size_{axis} = {value:.3f}" for axis, value in zip("xyz", size)),
        f"exhaustiveness = {EXHAUSTIVENESS}",
        f"num_modes = {NUM_MODES}",
        f"energy_range = {ENERGY_RANGE}",
    ]
    config = pdb_dir / "vina_config.txt"
    config.write_text("\n".join(lines) + "\n")
    return config


def best_affinity(vina_output: str) -> str:
    for line in vina_output.splitlines():
        fields = line.split()
        if len(fields) >= 2 and fields[0] == "1":
            return fields[1]
    return "неизвестно"


def dock(pdb_dir: Path) -> str:
    config = write_config(pdb_dir)
    result = subprocess.run(
        ["vina", "--config", config.name],
        cwd=pdb_dir,
        capture_output=True,
        text=True,
    )
    (pdb_dir / "vina_log.txt").write_text(result.stdout + "\n" + result.stderr)
    if result.returncode != 0:
        raise RuntimeError(f"Vina завершилась с ошибкой: {result.stderr.strip()}")
    return best_affinity(result.stdout)


def main() -> None:
    if shutil.which("vina") is None:
        raise SystemExit("Не найдена программа vina")

    failed = []
    for pdb in COMPLEXES:
        try:
            affinity = dock(DOCKING_DIR / pdb)
        except (OSError, RuntimeError) as error:
            print(f"{pdb}: {error}")
            failed.append(pdb)
            continue
        print(f"{pdb}: лучшая поза {affinity} ккал/моль")

    if failed:
        raise SystemExit(f"Докинг не удался: {', '.join(failed)}")


if __name__ == "__main__":
    main()
