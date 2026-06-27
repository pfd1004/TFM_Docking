from pathlib import Path
import re
import pandas as pd

DLG_DIR = Path("datos_sergio/dlg")
OUT = Path("resultados_docking_nanoparticulas.csv")

rows = []

def parse_name(path):
    """
    Espera nombres tipo:
    Au_1ao6.dlg
    Ag_2hav.dlg
    ZnO_1hzh.dlg
    """
    stem = path.stem
    parts = stem.split("_", 1)
    if len(parts) == 2:
        nano, proteina = parts
    else:
        nano, proteina = stem, ""
    return nano, proteina


def extract_energies_from_dlg(path):
    text = path.read_text(errors="ignore")
    lines = text.splitlines()

    energies = []

    # Opción 1: líneas típicas de AutoDock4
    patterns = [
        r"Estimated Free Energy of Binding\s*=\s*(-?\d+(?:\.\d+)?)",
        r"binding energy\s*[=:]\s*(-?\d+(?:\.\d+)?)",
    ]

    for line in lines:
        for pat in patterns:
            m = re.search(pat, line, re.IGNORECASE)
            if m:
                energies.append(float(m.group(1)))

    # Opción 2: bloque de histogram, como dijo Sergio
    if not energies:
        for i, line in enumerate(lines):
            if "histogram" in line.lower():
                block = lines[i:i+250]
                for bline in block:
                    nums = re.findall(r"-?\d+\.\d+", bline)
                    for n in nums:
                        val = float(n)
                        # Filtrado suave: energías de docking suelen ser negativas
                        if val < 0:
                            energies.append(val)
                break

    return energies


for dlg in sorted(DLG_DIR.glob("*.dlg")):
    nano, proteina = parse_name(dlg)
    energies = extract_energies_from_dlg(dlg)

    if energies:
        best = min(energies)
        n_values = len(energies)
        estado = "ok"
    else:
        best = None
        n_values = 0
        estado = "sin_energias"

    rows.append({
        "archivo": dlg.name,
        "nanoparticula": nano,
        "proteina": proteina.lower(),
        "binding_energy_min": best,
        "n_energias_detectadas": n_values,
        "estado": estado,
    })

df = pd.DataFrame(rows)
df.to_csv(OUT, index=False)

print(f"Creado {OUT}")
print(df)