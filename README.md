# 0. Notes

[Notes.pdf](https://github.com/user-attachments/files/32350071/looking.li_202609180154_34517.pdf)

# 1. Download HiqBind

```bash
wget -c -O hiqbind.tar.gz "https://ndownloader.figshare.com/files/52379345"
wget -c -O hiqbind_sm_metadata.csv "https://ndownloader.figshare.com/files/52379375"
```

Unpack

```bash
tar -xzvf hiqbind.tar.gz
```

```bash
uv sync
```

```bash
uv run extract_files.py
```

## 2. AutoDock Vina

Нужны `obabel` и `vina` в PATH.

```bash
uv run prepare_pdbqt.py
```

```bash
uv run run_vina.py
```

```bash
uv run calc_rmsd.py
```

## 3. Входные файлы для DiffDock и Boltz-1

```bash
uv run make_inputs.py
```

## 4. DiffDock

```bash
python -m inference --config default_inference_args.yaml \
    --protein_ligand_csv /path/to/re-docking/inputs/diffdock.csv \
    --out_dir /path/to/re-docking/predictions/diffdock
```

Ожидаемый результат: `predictions/diffdock/<pdb>/rank1.sdf`.

## 5. Boltz-1


```bash
boltz predict inputs/boltz --out_dir predictions/boltz \
    --model boltz1 --use_msa_server --output_format pdb
```


Ожидаемый результат: `predictions/boltz/**/<pdb>_model_0.pdb`, белок в цепи `A`, лиганд в цепи `L`.

## 6. Сбор поз

```bash
uv run collect_poses.py
```

## 7. RMSD и PoseBusters

```bash
uv run evaluate.py
```

Выход:
- `results/summary.md`
- `results/details.csv`
