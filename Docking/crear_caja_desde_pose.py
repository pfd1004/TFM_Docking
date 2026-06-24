from pathlib import Path
import sys

if len(sys.argv) != 3:
    print("Uso:")
    print("python crear_caja_desde_pose.py resultados\\proteina_ligando_out.pdbqt configs_refinados\\proteina_box.txt")
    sys.exit(1)

pose_file = Path(sys.argv[1])
out_config = Path(sys.argv[2])
out_config.parent.mkdir(exist_ok=True)

xs, ys, zs = [], [], []
in_model_1 = False

with open(pose_file, "r", errors="ignore") as f:
    for line in f:
        if line.startswith("MODEL"):
            parts = line.split()
            if len(parts) > 1 and parts[1] == "1":
                in_model_1 = True
            elif in_model_1:
                break

        if in_model_1 and line.startswith(("ATOM", "HETATM")):
            try:
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
            except ValueError:
                pass

if not xs:
    print(f"No he encontrado coordenadas en MODEL 1 de {pose_file}")
    sys.exit(1)

center_x = sum(xs) / len(xs)
center_y = sum(ys) / len(ys)
center_z = sum(zs) / len(zs)

config_text = f"""center_x = {center_x:.3f}
center_y = {center_y:.3f}
center_z = {center_z:.3f}

size_x = 28
size_y = 28
size_z = 28

exhaustiveness = 16
num_modes = 10
energy_range = 4
cpu = 6
seed = 12345
"""

out_config.write_text(config_text)
print(f"Creado {out_config}")
print(config_text)