Download HiqBind

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

```bash
uv run prepare_pdbqt.py
```

```bash
uv run run_vina.py
```

```bash
uv run calc_rmsd.py
```


Результаты в docking/results.csv
