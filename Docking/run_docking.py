from pathlib import Path
import subprocess
import shutil

proteinas = Path("Proteinas_clean")
ligandos = Path("Ligandos_pdbqt")
configs = Path("configs_refinados")
resultados = Path("resultados_final")
logs = Path("logs_final")

resultados.mkdir(exist_ok=True)
logs.mkdir(exist_ok=True)

vina = shutil.which("vina.exe")
if vina is None:
    local_vina = Path("vina.exe")
    if local_vina.exists():
        vina = str(local_vina.resolve())
    else:
        raise FileNotFoundError("No encuentro vina.exe. Pon vina.exe en la carpeta Docking.")

fallidos = []

for prot in proteinas.glob("*.pdbqt"):
    pname = prot.stem
    config = configs / f"{pname}_box.txt"

    if not config.exists():
        print(f"Saltando {pname}: no existe {config}")
        continue

    for lig in ligandos.glob("*.pdbqt"):
        lname = lig.stem

        out = resultados / f"{pname}_{lname}_out.pdbqt"
        log = logs / f"{pname}_{lname}_log.txt"

        # Si ya existe resultado, lo saltamos para no repetir
        if out.exists() and log.exists():
            print(f"Ya existe, saltando: {pname} + {lname}")
            continue

        print(f"Docking final: {pname} + {lname}")

        cmd = [
            vina,
            "--config", str(config),
            "--receptor", str(prot),
            "--ligand", str(lig),
            "--out", str(out),
        ]

        with open(log, "w", errors="ignore") as f:
            result = subprocess.run(
                cmd,
                stdout=f,
                stderr=subprocess.STDOUT,
                check=False
            )

        if result.returncode != 0:
            print(f"FALLÓ: {pname} + {lname}. Mira {log}")
            fallidos.append(f"{pname},{lname},{log}")

if fallidos:
    with open("dockings_fallidos.csv", "w") as f:
        f.write("proteina,ligando,log\n")
        for item in fallidos:
            f.write(item + "\n")

    print("\nAlgunos dockings fallaron. Mira dockings_fallidos.csv")
else:
    print("\nTodos los dockings finales han terminado correctamente.")