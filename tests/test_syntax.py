"""Every shipped script must parse or compile in its own language."""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
C_SOURCES = sorted(SCRIPTS.glob("*.c"))
R_SCRIPTS = sorted(SCRIPTS.glob("*.R"))
PY_SCRIPTS = sorted(SCRIPTS.glob("*.py"))


@pytest.mark.parametrize("path", PY_SCRIPTS, ids=lambda p: p.name)
def test_python_parses(path):
    subprocess.run([sys.executable, "-m", "py_compile", str(path)], check=True)


@pytest.mark.parametrize("path", R_SCRIPTS, ids=lambda p: p.name)
def test_r_parses(path):
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript not installed")
    subprocess.run([rscript, "-e", f"invisible(parse('{path}'))"], check=True)


@pytest.mark.parametrize("path", C_SOURCES, ids=lambda p: p.name)
def test_c_compiles(path):
    cc = shutil.which("cc")
    if cc is None:
        pytest.skip("no C compiler")
    subprocess.run([cc, "-O3", "-std=c99", "-fsyntax-only", str(path)], check=True)


@pytest.mark.parametrize("path", sorted((ROOT / "notebooks").glob("*.ipynb")), ids=lambda p: p.name)
def test_notebook_is_valid_json(path):
    nb = json.loads(path.read_text())
    assert nb["nbformat"] >= 4
    assert nb["cells"]
