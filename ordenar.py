from pathlib import Path
import shutil
import datetime

ROOT = Path.cwd()
STAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
BACKUP = ROOT / f"_backup_ordenar_{STAMP}"

print("Raíz detectada:", ROOT)
print("Backup:", BACKUP)

# =========================
# CARPETAS FINALES
# =========================

DIRS = [
    "datasets",
    "scripts/modelos",
    "scripts/membrana",
    "scripts/nanoparticulas",
    "scripts/diagnostico",
    "scripts/utilidades",
    "resultados/modelos",
    "resultados/diagnosticos",
    "resultados/nanoparticulas",
    "resultados/tablas",
    "resultados/xlsx_antiguos",
    "app",
    "_paquetes",
]

for d in DIRS:
    (ROOT / d).mkdir(parents=True, exist_ok=True)

BACKUP.mkdir(parents=True, exist_ok=True)


# =========================
# FUNCIONES
# =========================

def safe_copy_to_backup(path):
    if path.is_file():
        dst = BACKUP / path.name
        if not dst.exists():
            shutil.copy2(path, dst)


def safe_move(src, dst):
    src = Path(src)
    dst = Path(dst)

    if not src.exists():
        return

    dst.parent.mkdir(parents=True, exist_ok=True)

    if dst.exists():
        stem = dst.stem
        suffix = dst.suffix
        parent = dst.parent
        i = 1
        while True:
            candidate = parent / f"{stem}_old{i}{suffix}"
            if not candidate.exists():
                dst = candidate
                break
            i += 1

    print(f"MOVIENDO: {src.name} -> {dst.relative_to(ROOT)}")
    shutil.move(str(src), str(dst))


def disable_excel_in_py(path):
    """
    Comenta cualquier bloque .to_excel(...)
    para que los scripts solo guarden CSV.
    """
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        text = path.read_text(encoding="latin-1")

    lines = text.splitlines(keepends=True)
    new_lines = []

    skipping = False
    balance = 0
    changed = False

    for line in lines:
        if skipping:
            balance += line.count("(") - line.count(")")
            if line.lstrip().startswith("#"):
                new_lines.append(line)
            else:
                new_lines.append("# XLSX_DESACTIVADO " + line)
            changed = True
            if balance <= 0:
                skipping = False
            continue

        if ".to_excel(" in line:
            balance = line.count("(") - line.count(")")
            if line.lstrip().startswith("#"):
                new_lines.append(line)
            else:
                new_lines.append("# XLSX_DESACTIVADO " + line)
            changed = True
            if balance > 0:
                skipping = True
        else:
            new_lines.append(line)

    if changed:
        path.write_text("".join(new_lines), encoding="utf-8")
        print("XLSX desactivado en:", path.relative_to(ROOT))


def patch_paths_in_py(path):
    """
    Cambia rutas antiguas a la nueva estructura.
    Se asume que los scripts se ejecutan desde la raíz TFM_DOCKING.
    """
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        text = path.read_text(encoding="latin-1")

    replacements = {
        "dataset_descriptores_membrana.csv": "datasets/dataset_descriptores_membrana.csv",
        "dataset_descriptores_membrana_corregido.csv": "datasets/dataset_descriptores_membrana_corregido.csv",
        "dataset_descriptores.csv": "datasets/dataset_descriptores.csv",
        "membrana_features.csv": "datasets/membrana_features.csv",
        "membrana_features_corregido.csv": "datasets/membrana_features_corregido.csv",

        "correlaciones_features_targets.csv": "resultados/diagnosticos/correlaciones_features_targets.csv",
        "diagnostico_columnas.csv": "resultados/diagnosticos/diagnostico_columnas.csv",

        "resultados_docking_nanoparticulas_corregido.csv": "resultados/nanoparticulas/resultados_docking_nanoparticulas_corregido.csv",
        "resultados_docking_nanoparticulas.csv": "resultados/nanoparticulas/resultados_docking_nanoparticulas.csv",
        "ranking_contextual_ligando_nano_proteina.csv": "resultados/nanoparticulas/ranking_contextual_ligando_nano_proteina.csv",
        "dataset_ligando_nano_proteina_contexto_long.csv": "resultados/nanoparticulas/dataset_ligando_nano_proteina_contexto_long.csv",
        "tabla_nano_proteina_contexto_wide.csv": "resultados/nanoparticulas/tabla_nano_proteina_contexto_wide.csv",
        "tabla_nano_proteina_wide.csv": "resultados/nanoparticulas/tabla_nano_proteina_wide.csv",

        "resultados_comparacion_membrana": "resultados/modelos/resultados_comparacion_membrana",
        "resultados_modelos": "resultados/modelos/resultados_modelos",
        "modelo_final_rdkit": "resultados/modelos/modelo_final_rdkit",
    }

    original = text

    for old, new in replacements.items():
        # Evita duplicar si ya está cambiado
        text = text.replace(f'"{old}"', f'"{new}"')
        text = text.replace(f"'{old}'", f"'{new}'")

    if text != original:
        path.write_text(text, encoding="utf-8")
        print("Rutas actualizadas en:", path.relative_to(ROOT))


def destino_script(path):
    name = path.name.lower()

    if name == "ordenar_y_quitar_xlsx.py":
        return None

    if name.startswith("app_") or name == "app_tfm.py":
        return ROOT / "app" / path.name

    if "entrenar" in name or "modelo" in name:
        return ROOT / "scripts" / "modelos" / path.name

    if "membrana" in name or "logp" in name or "cosmo" in name:
        return ROOT / "scripts" / "membrana" / path.name

    if "nano" in name or "dlg" in name or "case" in name:
        return ROOT / "scripts" / "nanoparticulas" / path.name

    if "diagnost" in name or "correl" in name:
        return ROOT / "scripts" / "diagnostico" / path.name

    return ROOT / "scripts" / "utilidades" / path.name


def destino_csv(path):
    name = path.name.lower()

    if name.startswith("dataset_descriptores") or name.startswith("membrana_features"):
        return ROOT / "datasets" / path.name

    if "diagnostico" in name or "correlaciones" in name:
        return ROOT / "resultados" / "diagnosticos" / path.name

    if (
        "nano" in name
        or "nanoparticula" in name
        or "ranking_contextual" in name
        or "case_study" in name
        or "ligando_nano" in name
    ):
        return ROOT / "resultados" / "nanoparticulas" / path.name

    return ROOT / "resultados" / "tablas" / path.name


# =========================
# 1. BACKUP DE SCRIPTS Y CSV
# =========================

for p in ROOT.glob("*.py"):
    safe_copy_to_backup(p)

for p in ROOT.glob("*.csv"):
    safe_copy_to_backup(p)

print("\nBackup básico hecho.")


# =========================
# 2. QUITAR XLSX EN SCRIPTS DE RAÍZ
# =========================

for p in ROOT.glob("*.py"):
    disable_excel_in_py(p)


# =========================
# 3. MOVER SCRIPTS
# =========================

for p in list(ROOT.glob("*.py")):
    dst = destino_script(p)
    if dst is not None:
        safe_move(p, dst)


# =========================
# 4. MOVER CSV
# =========================

for p in list(ROOT.glob("*.csv")):
    safe_move(p, destino_csv(p))


# =========================
# 5. MOVER XLSX ANTIGUOS
# =========================

for p in list(ROOT.glob("*.xlsx")):
    safe_move(p, ROOT / "resultados" / "xlsx_antiguos" / p.name)


# =========================
# 6. MOVER ZIPS
# =========================

for p in list(ROOT.glob("*.zip")):
    safe_move(p, ROOT / "_paquetes" / p.name)


# =========================
# 7. MOVER CARPETAS DE RESULTADOS
# =========================

for d in list(ROOT.glob("resultados_comparacion*")):
    if d.is_dir():
        safe_move(d, ROOT / "resultados" / "modelos" / d.name)

for d in list(ROOT.glob("resultados_modelos*")):
    if d.is_dir():
        safe_move(d, ROOT / "resultados" / "modelos" / d.name)

if (ROOT / "modelo_final_rdkit").exists():
    safe_move(ROOT / "modelo_final_rdkit", ROOT / "resultados" / "modelos" / "modelo_final_rdkit")

if (ROOT / "TFM_cierre_modelos_app").exists():
    safe_move(ROOT / "TFM_cierre_modelos_app", ROOT / "_paquetes" / "TFM_cierre_modelos_app")


# =========================
# 8. ACTUALIZAR RUTAS EN SCRIPTS YA MOVIDOS
# =========================

for folder in ["scripts", "app"]:
    for p in (ROOT / folder).rglob("*.py"):
        disable_excel_in_py(p)
        patch_paths_in_py(p)


