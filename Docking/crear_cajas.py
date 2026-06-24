from pathlib import Path

proteinas = Path("Proteinas")
configs = Path("configs")
configs.mkdir(exist_ok=True)

padding = 8.0

for pdbqt in proteinas.glob("*.pdbqt"):
    xs, ys, zs = [], [], []

    with open(pdbqt, "r", errors="ignore") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    xs.append(x)
                    ys.append(y)
                    zs.append(z)
                except ValueError:
                    pass

    if not xs:
        print(f"No se encontraron átomos en {pdbqt.name}")
        continue

    center_x = (min(xs) + max(xs)) / 2
    center_y = (min(ys) + max(ys)) / 2
    center_z = (min(zs) + max(zs)) / 2

    size_x = (max(xs) - min(xs)) + padding
    size_y = (max(ys) - min(ys)) + padding
    size_z = (max(zs) - min(zs)) + padding

    config_text = f"""center_x = {center_x:.3f}
center_y = {center_y:.3f}
center_z = {center_z:.3f}

size_x = {size_x:.3f}
size_y = {size_y:.3f}
size_z = {size_z:.3f}

exhaustiveness = 8
num_modes = 10
energy_range = 4
cpu = 6
seed = 12345
"""

    out = configs / f"{pdbqt.stem}_box.txt"
    out.write_text(config_text)
    print(f"Creado {out}")