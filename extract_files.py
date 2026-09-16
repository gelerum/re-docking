import shutil
from pathlib import Path

COMPLEXES = ["181l", "182l", "183l", "184l", "185l"]
RAW_DIR = Path("raw_data_hiq_sm")
OUT_DIR = Path("docking")


def only_one(items: list[Path], description: str) -> Path:
    if len(items) != 1:
        names = [item.name for item in items]
        raise RuntimeError(f"ожидался один {description}, найдено: {names}")
    return items[0]


def extract(pdb: str) -> None:
    subdirs = [d for d in (RAW_DIR / pdb).iterdir() if d.is_dir()]
    source = only_one(subdirs, "подкаталог")
    protein = only_one(list(source.glob("*_protein_refined.pdb")), "файл белка")
    ligand = only_one(list(source.glob("*_ligand_refined.sdf")), "файл лиганда")

    target = OUT_DIR / pdb
    target.mkdir(parents=True, exist_ok=True)
    shutil.copy2(protein, target / "receptor.pdb")
    shutil.copy2(ligand, target / "native_ligand.sdf")


def main() -> None:
    if not RAW_DIR.is_dir():
        raise SystemExit(f"Не найден каталог {RAW_DIR}, запустите скрипт из корня проекта")

    for pdb in COMPLEXES:
        try:
            extract(pdb)
        except (OSError, RuntimeError) as error:
            print(f"{pdb}: пропущен, {error}")
            continue
        print(f"{pdb}: скопированы receptor.pdb и native_ligand.sdf")


if __name__ == "__main__":
    main()
