from pathlib import Path
import csv
import re

logs = Path("logs_final")
rows = []

for log in logs.glob("*_log.txt"):
    text = log.read_text(errors="ignore")

    name = log.stem.replace("_log", "")
    parts = name.split("_", 1)

    proteina = parts[0]
    ligando = parts[1] if len(parts) > 1 else ""

    best_affinity = None
    estado = "fallido"

    for line in text.splitlines():
        m = re.match(r"\s*1\s+(-?\d+(?:\.\d+)?)", line)
        if m:
            best_affinity = float(m.group(1))
            estado = "ok"
            break

    if "PDBQT parsing error" in text:
        estado = "error_pdbqt"
    elif "ERROR" in text.upper() and best_affinity is None:
        estado = "error"

    rows.append({
        "proteina": proteina,
        "ligando": ligando,
        "afinidad_kcal_mol": best_affinity,
        "estado": estado,
        "log": str(log)
    })

rows.sort(
    key=lambda x: (
        x["proteina"],
        x["afinidad_kcal_mol"] if x["afinidad_kcal_mol"] is not None else 999
    )
)

with open("tabla_resultados_docking_final.csv", "w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=["proteina", "ligando", "afinidad_kcal_mol", "estado", "log"]
    )
    writer.writeheader()
    writer.writerows(rows)

print("Creado tabla_resultados_docking_final.csv")

ok = sum(1 for r in rows if r["estado"] == "ok")
fallidos = len(rows) - ok

print(f"Dockings OK: {ok}")
print(f"Dockings fallidos/error: {fallidos}")