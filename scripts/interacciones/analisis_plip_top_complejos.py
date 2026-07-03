# analisis_plip_top_complejos.py
# ============================================================
# Selecciona los mejores complejos ligando-proteína y prepara/runea PLIP.
#
# Ejecutar desde la raíz TFM_DOCKING:
#   python scripts\interacciones\analisis_plip_top_complejos.py
#
# Requisitos opcionales:
#   pip install plip
# o instalación por conda, si la tienes.
#
# Si PLIP no está instalado, el script igualmente genera:
#   resultados/plip/candidatos_plip.csv
#   resultados/plip/comandos_plip.bat
#
# IMPORTANTE:
# El script intenta encontrar automáticamente:
#   - receptores .pdbqt/.pdb
#   - poses de ligando .pdbqt/.pdb
# en carpetas como Docking/, datos_sergio/ y la raíz.
# Si no encuentra poses, te deja la tabla de candidatos para buscarlas manualmente.
# ============================================================

from pathlib import Path
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

ROOT = Path.cwd()

DATASET_CANDIDATES = [
    ROOT / "datasets" / "dataset_descriptores_membrana.csv",
    ROOT / "dataset_descriptores_membrana.csv",
]

OUTDIR = ROOT / "resultados" / "plip"
COMPLEX_DIR = OUTDIR / "complejos_pdb"
PLIP_DIR = OUTDIR / "plip_outputs"

OUTDIR.mkdir(parents=True, exist_ok=True)
COMPLEX_DIR.mkdir(parents=True, exist_ok=True)
PLIP_DIR.mkdir(parents=True, exist_ok=True)

TOP_PER_PROTEIN = 3

TARGETS = {
    "1ao6": "y_1ao6",
    "1hzh": "y_1hzh",
    "2hav": "y_2hav",
    "3ghg": "y_3ghg",
}

SEARCH_DIRS = [
    ROOT / "Docking",
    ROOT / "datos_sergio",
    ROOT / "resultados",
    ROOT,
]


def norm(s):
    s = str(s).lower()
    s = s.replace(".pdbqt", "")
    s = s.replace(".pdb", "")
    s = s.replace(".sdf", "")
    s = s.replace(".mol2", "")
    s = re.sub(r"[^a-z0-9]+", "", s)
    return s


def find_first_existing(candidates):
    for p in candidates:
        if p.exists():
            return p
    return None


def load_dataset():
    path = find_first_existing(DATASET_CANDIDATES)
    if path is None:
        raise FileNotFoundError("No encuentro dataset_descriptores_membrana.csv")
    print("Leyendo dataset:", path)
    return pd.read_csv(path)


def select_candidates(df):
    rows = []

    for prot, ycol in TARGETS.items():
        if ycol not in df.columns:
            print(f"AVISO: no existe {ycol}, salto {prot}")
            continue

        sub = df[["ligando", ycol]].dropna().copy()
        sub[ycol] = pd.to_numeric(sub[ycol], errors="coerce")
        sub = sub.dropna().sort_values(ycol, ascending=True).head(TOP_PER_PROTEIN)

        for _, r in sub.iterrows():
            rows.append({
                "ligando": r["ligando"],
                "proteina": prot,
                "energia_docking": r[ycol],
                "target_col": ycol,
            })

    cand = pd.DataFrame(rows)
    cand = cand.drop_duplicates(subset=["ligando", "proteina"])
    return cand


def list_structure_files():
    files = []
    for d in SEARCH_DIRS:
        if d.exists():
            files.extend(list(d.rglob("*.pdbqt")))
            files.extend(list(d.rglob("*.pdb")))
    return files


def find_receptor(protein, all_files):
    pn = norm(protein)

    # Prioridad: archivos que parezcan receptor/proteína y contengan el código
    scored = []
    for f in all_files:
        fn = norm(f.name)
        pathn = norm(str(f.parent))
        if pn not in fn and pn not in pathn:
            continue

        score = 0
        if "receptor" in fn or "protein" in fn or "proteina" in fn:
            score += 5
        if f.suffix.lower() == ".pdbqt":
            score += 2
        if fn == pn:
            score += 4
        if "lig" in fn or "pose" in fn or "out" in fn:
            score -= 4

        scored.append((score, f))

    if not scored:
        return None

    scored.sort(key=lambda x: x[0], reverse=True)
    return scored[0][1]


def find_pose(ligand, protein, all_files):
    ln = norm(ligand)
    pn = norm(protein)

    scored = []
    for f in all_files:
        fn = norm(f.name)
        pathn = norm(str(f.parent))

        if ln not in fn and ln not in pathn:
            continue
        if pn not in fn and pn not in pathn:
            continue

        score = 0
        if f.suffix.lower() == ".pdbqt":
            score += 5
        if "out" in fn or "pose" in fn or "vina" in fn or "ligand" in fn or "ligando" in fn:
            score += 4
        if "receptor" in fn or "protein" in fn or "proteina" in fn:
            score -= 10
        if pn in fn:
            score += 2
        if ln in fn:
            score += 2

        scored.append((score, f))

    if not scored:
        return None

    scored.sort(key=lambda x: x[0], reverse=True)
    return scored[0][1]


def guess_element(atom_name, line):
    # PDBQT suele tener elemento como último token; si no, usar atom_name.
    parts = line.split()
    if parts:
        last = parts[-1]
        if re.fullmatch(r"[A-Za-z]{1,2}", last):
            return last.capitalize()

    letters = re.sub(r"[^A-Za-z]", "", atom_name).strip()
    if not letters:
        return "C"
    if len(letters) >= 2 and letters[:2].capitalize() in {"Cl", "Br", "Zn", "Ag", "Au", "Na", "Ca", "Fe", "Mg"}:
        return letters[:2].capitalize()
    return letters[0].upper()


def pdbqt_atoms_to_pdb_lines(path, ligand=False):
    lines = []
    serial = 1

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            if not (raw.startswith("ATOM") or raw.startswith("HETATM")):
                continue

            atom_name = raw[12:16].strip() or f"C{serial}"
            x = raw[30:38].strip()
            y = raw[38:46].strip()
            z = raw[46:54].strip()

            try:
                x = float(x)
                y = float(y)
                z = float(z)
            except ValueError:
                continue

            element = guess_element(atom_name, raw)

            if ligand:
                record = "HETATM"
                resname = "LIG"
                chain = "Z"
                resseq = 1
            else:
                record = "ATOM  "
                resname = (raw[17:20].strip() or "RES")[:3]
                chain = (raw[21:22].strip() or "A")[:1]
                try:
                    resseq = int(raw[22:26].strip())
                except ValueError:
                    resseq = 1

            line = (
                f"{record}{serial:5d} {atom_name:<4s} {resname:>3s} {chain}"
                f"{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
                f"  1.00  0.00          {element:>2s}\n"
            )
            lines.append(line)
            serial += 1

    return lines


def pdb_atoms_to_pdb_lines(path, ligand=False):
    lines = []
    serial = 1

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            if not (raw.startswith("ATOM") or raw.startswith("HETATM")):
                continue

            if ligand:
                atom_name = raw[12:16].strip() or f"C{serial}"
                x = float(raw[30:38])
                y = float(raw[38:46])
                z = float(raw[46:54])
                element = raw[76:78].strip() or guess_element(atom_name, raw)
                line = (
                    f"HETATM{serial:5d} {atom_name:<4s} LIG Z"
                    f"{1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
                    f"  1.00  0.00          {element:>2s}\n"
                )
                lines.append(line)
                serial += 1
            else:
                lines.append(raw if raw.endswith("\n") else raw + "\n")

    return lines


def read_structure_as_pdb_lines(path, ligand=False):
    if path.suffix.lower() == ".pdbqt":
        return pdbqt_atoms_to_pdb_lines(path, ligand=ligand)
    return pdb_atoms_to_pdb_lines(path, ligand=ligand)


def make_complex(receptor_path, ligand_pose_path, out_path):
    rec_lines = read_structure_as_pdb_lines(receptor_path, ligand=False)
    lig_lines = read_structure_as_pdb_lines(ligand_pose_path, ligand=True)

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("REMARK Complex generated for PLIP\n")
        f.write(f"REMARK receptor={receptor_path}\n")
        f.write(f"REMARK ligand_pose={ligand_pose_path}\n")
        for line in rec_lines:
            f.write(line)
        f.write("TER\n")
        for line in lig_lines:
            f.write(line)
        f.write("END\n")

    return len(rec_lines), len(lig_lines)


def run_plip(complex_path, outdir):
    plip_exe = shutil.which("plip")

    if plip_exe is None:
        return False, "PLIP no está instalado o no está en PATH."

    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        plip_exe,
        "-f", str(complex_path),
        "-o", str(outdir),
        "-x",
        "-t",
    ]

    print("Ejecutando:", " ".join(cmd))

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        return False, result.stderr[-1000:]

    return True, result.stdout[-1000:]


def parse_plip_xml(outdir):
    xmls = list(outdir.rglob("*.xml"))
    if not xmls:
        return {}

    xml = xmls[0]
    counts = {}

    try:
        root = ET.parse(xml).getroot()
    except Exception:
        return {}

    # PLIP XML contiene bloques de interacciones con nombres tipo hydrophobic_interactions, hydrogen_bonds, etc.
    for interactions in root.iter():
        if interactions.tag.lower().endswith("interactions"):
            for group in list(interactions):
                tag = group.tag
                n = len(list(group))
                if n > 0:
                    counts[tag] = counts.get(tag, 0) + n

    return counts


def main():
    df = load_dataset()
    cand = select_candidates(df)

    if cand.empty:
        raise SystemExit("No hay candidatos.")

    all_files = list_structure_files()
    print("Archivos estructurales encontrados:", len(all_files))

    rows = []
    commands = []

    for _, row in cand.iterrows():
        ligand = row["ligando"]
        protein = row["proteina"]

        receptor = find_receptor(protein, all_files)
        pose = find_pose(ligand, protein, all_files)

        safe_lig = re.sub(r"[^A-Za-z0-9_]+", "_", str(ligand))[:80]
        complex_name = f"{protein}_{safe_lig}.pdb"
        complex_path = COMPLEX_DIR / complex_name

        status = "pendiente"
        msg = ""
        n_rec = np.nan
        n_lig = np.nan

        if receptor is None:
            status = "sin_receptor"
            msg = "No se encontró receptor."
        elif pose is None:
            status = "sin_pose"
            msg = "No se encontró pose de docking."
        else:
            try:
                n_rec, n_lig = make_complex(receptor, pose, complex_path)
                status = "complex_ok"

                plip_out = PLIP_DIR / f"{protein}_{safe_lig}"
                ok, msg = run_plip(complex_path, plip_out)

                if ok:
                    status = "plip_ok"
                else:
                    status = "complex_ok_sin_plip"

                commands.append(
                    f'plip -f "{complex_path}" -o "{plip_out}" -x -t'
                )

            except Exception as e:
                status = "error_complex"
                msg = str(e)

        rows.append({
            **row.to_dict(),
            "receptor_encontrado": str(receptor) if receptor else "",
            "pose_encontrada": str(pose) if pose else "",
            "complex_pdb": str(complex_path) if complex_path.exists() else "",
            "n_atomos_receptor": n_rec,
            "n_atomos_ligando": n_lig,
            "estado": status,
            "mensaje": msg,
        })

    out = pd.DataFrame(rows)
    out.to_csv(OUTDIR / "candidatos_plip.csv", index=False)

    with open(OUTDIR / "comandos_plip.bat", "w", encoding="utf-8") as f:
        f.write("@echo off\n")
        for cmd in commands:
            f.write(cmd + "\n")
        f.write("pause\n")

    # Resumen de XML si PLIP ha corrido.
    summary_rows = []
    for _, r in out.iterrows():
        if r["estado"] != "plip_ok":
            continue

        safe_lig = re.sub(r"[^A-Za-z0-9_]+", "_", str(r["ligando"]))[:80]
        plip_out = PLIP_DIR / f"{r['proteina']}_{safe_lig}"
        counts = parse_plip_xml(plip_out)

        base = {
            "ligando": r["ligando"],
            "proteina": r["proteina"],
            "energia_docking": r["energia_docking"],
        }
        base.update(counts)
        summary_rows.append(base)

    if summary_rows:
        resumen = pd.DataFrame(summary_rows).fillna(0)
        resumen.to_csv(OUTDIR / "resumen_plip_interacciones.csv", index=False)
        print("Creado:", OUTDIR / "resumen_plip_interacciones.csv")
    else:
        print("No hay resumen PLIP todavía. Si no tienes PLIP instalado, usa comandos_plip.bat tras instalarlo.")

    print("\nCreados:")
    print(OUTDIR / "candidatos_plip.csv")
    print(OUTDIR / "comandos_plip.bat")
    print("\nEstados:")
    print(out["estado"].value_counts().to_string())

    print("\nSi ves muchos 'sin_pose', hay que decirle al script dónde están los PDBQT de salida de Vina.")


if __name__ == "__main__":
    main()
