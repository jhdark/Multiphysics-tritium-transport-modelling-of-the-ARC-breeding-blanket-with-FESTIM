"""Unpack the zipped input data on first use.

The large input files are stored as zip archives in data/ (each under GitHub's
100 MB file limit). Archive paths are relative to the repository root, so they
extract to the locations the scripts expect.
"""

from pathlib import Path
import zipfile

REPO_ROOT = Path(__file__).resolve().parent

OPENFOAM_CASE_ARCHIVES = [
    "data/NewOF_symmBCs_mesh.zip",  # arc.foam, system/, constant/
    "data/NewOF_symmBCs_fields.zip",  # 745/ (U, T, nut, tritProd, heatGen)
]
ORNL_ARCHIVES = ["data/OpenFOAM_results_ORNL.zip"]


def ensure_unzipped(archive):
    """Extract archive into the repository root unless all its files exist"""
    with zipfile.ZipFile(REPO_ROOT / archive) as zf:
        if all((REPO_ROOT / name).exists() for name in zf.namelist()):
            return
        print(f"Extracting {archive}")
        zf.extractall(REPO_ROOT)


def ensure_openfoam_case():
    """OpenFOAM case (data/NewOF_symmBCs) used as input to all FESTIM runs"""
    for archive in OPENFOAM_CASE_ARCHIVES:
        ensure_unzipped(archive)


def ensure_ornl_data():
    """ORNL OpenFOAM passive scalar results (results/OpenFOAM_results_ORNL)"""
    for archive in ORNL_ARCHIVES:
        ensure_unzipped(archive)
