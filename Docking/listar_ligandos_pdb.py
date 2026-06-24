from pathlib import Path
import urllib.request
import csv
from collections import defaultdict

pdb_ids = ["1AO6", "1HZH", "2HAV", "3GHG"]

outdir = Path("PDB_originales")
outdir.mkdir(exist_ok=True)

# Cosas que no queremos usar para definir caja
ignorar = {
    "HOH", "WAT", "DOD",
    "NA", "CL", "K", "CA", "MG", "ZN", "MN", "FE", "CU",
    "SO4", "PO4", "ACT", "EDO", "GOL"
}

rows = []

for pdb_id in pdb_ids:
    pdb_file = outdir / f"{pdb_id}.pdb"

    if not pdb_file.exists():
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        print(f"Descargando {pdb_id}...")
        urllib.request.urlretrieve(url, pdb_file)

    grupos = defaultdict(list)

    with open(pdb_file, "r", errors="ignore") as f:
        for line in f:
            if line.startswith("HETATM"):
                resname = line[17:20].strip()
                chain = line[21].strip()
                resseq = line[22:26].strip()
                icode = line[26].strip()

                if resname in ignorar:
                    continue

                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue

                key = (pdb_id, resname, chain, resseq, icode)
                grupos[key].append((x, y, z))

    if not grupos:
        print(f"\n{pdb_id}: no he encontrado ligandos HETATM útiles.")
        continue

    print(f"\n{pdb_id}: posibles ligandos/cofactores:")
    for (pid, resname, chain, resseq, icode), coords in grupos.items():
        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        zs = [c[2] for c in coords]

        cx = sum(xs) / len(xs)
        cy = sum(ys) / len(ys)
        cz = sum(zs) / len(zs)

        row = {
            "pdb": pid,
            "resname": resname,
            "chain": chain,
            "resseq": resseq,
            "icode": icode,
            "n_atoms": len(coords),
            "center_x": round(cx, 3),
            "center_y": round(cy, 3),
            "center_z": round(cz, 3),
        }
        rows.append(row)

        print(
            f"  {resname} cadena {chain} residuo {resseq} "
            f"({len(coords)} átomos) -> "
            f"center_x={cx:.3f}, center_y={cy:.3f}, center_z={cz:.3f}"
        )

with open("hetero_centros.csv", "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "pdb", "resname", "chain", "resseq", "icode",
            "n_atoms", "center_x", "center_y", "center_z"
        ]
    )
    writer.writeheader()
    writer.writerows(rows)

print("\nCreado hetero_centros.csv")