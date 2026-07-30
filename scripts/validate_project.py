#!/usr/bin/env python3
"""Validate the public ChIP-seq repository and headline results."""

from pathlib import Path
import gzip
import hashlib
import re

import nbformat
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks" / "01_eralpha_chipseq_analysis.ipynb"
notebook = nbformat.read(NOTEBOOK, as_version=4)
nbformat.validate(notebook)
code_cells = [cell for cell in notebook.cells if cell.cell_type == "code"]
errors = [
    output
    for cell in code_cells
    for output in cell.get("outputs", [])
    if output.get("output_type") == "error"
]
assert code_cells and all(cell.execution_count is not None for cell in code_cells)
assert not errors

expected_hashes = {
    "GSM365925_ER_minus_ligand_peaks.bed.gz": "a4c716754e9b313da8829e7a2dbc4ae1f90586e20cc1bcfdfae9cfab55a979b5",
    "GSM365926_ER_E2_peaks.bed.gz": "dffb9c5cf2cf8a6c7e0b010554764def9186b18d290a562f478a4bd33a000856",
    "Homo_sapiens.NCBI36.54.gtf.gz": "e8c8429741fe23f453e476f2d698d9eae3ae747a3d36d3cb372aa21f5ac07066",
    "MSigDB_Hallmark_2020.gmt": "4275592957a1587652092bb398cf77216fde5b8daa2aedaa0e016f7d10bbdb81",
}
for name, expected in expected_hashes.items():
    observed = hashlib.sha256((ROOT / "data" / name).read_bytes()).hexdigest()
    assert observed == expected, name

def bed_count(path: Path) -> int:
    with gzip.open(path, "rt") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith(("track", "browser", "#")))

assert bed_count(ROOT / "data" / "GSM365925_ER_minus_ligand_peaks.bed.gz") == 240
assert bed_count(ROOT / "data" / "GSM365926_ER_E2_peaks.bed.gz") == 10205

comparison = pd.read_csv(ROOT / "results" / "tables" / "condition_specific_peaks.csv")
counts = dict(zip(comparison["category"], comparison["count"], strict=True))
assert counts["E2-specific peaks"] == 9948

enrichment = pd.read_csv(ROOT / "results" / "tables" / "hallmark_enrichment.csv")
assert int((enrichment["fdr"] < 0.05).sum()) == 2
assert enrichment.iloc[0]["pathway"] == "Estrogen Response Early"

required = [
    ROOT / "reports" / "01_eralpha_chipseq_analysis.html",
    ROOT / "results" / "figures" / "01_alignment_qc.png",
    ROOT / "results" / "figures" / "02_peak_landscape.png",
    ROOT / "results" / "figures" / "03_motif_enrichment.png",
    ROOT / "results" / "figures" / "04_hallmark_enrichment.png",
]
for path in required:
    assert path.exists() and path.stat().st_size > 10_000, path

forbidden = re.compile(
    r"/Users/|/home/|r0[0-9]{6}|student/task|projects_biomed|anonymous|submission",
    flags=re.IGNORECASE,
)
for path in ROOT.rglob("*"):
    if path.resolve() == Path(__file__).resolve():
        continue
    if path.is_file() and path.suffix.lower() in {".md", ".py", ".sh", ".ipynb", ".txt", ".csv"}:
        assert not forbidden.search(path.read_text(errors="replace")), path

print("PASS: executed notebook, 0 cell errors, 9,948 E2-specific peaks and 2 FDR-significant Hallmark pathways.")
