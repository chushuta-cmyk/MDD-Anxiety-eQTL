"""Repository hygiene: nothing personal, nothing oversized, nothing unbuildable."""
from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
TEXT_SUFFIXES = {".py", ".R", ".c", ".sh", ".md", ".yml", ".yaml", ".txt", ".cff", ".gitignore"}
CJK = re.compile(r"[\u4e00-\u9fff]")
MAX_BYTES = 2_000_000


def text_files():
    for p in ROOT.rglob("*"):
        if p.is_file() and (p.suffix in TEXT_SUFFIXES or p.name == ".gitignore"):
            if "data" in p.relative_to(ROOT).parts[:1]:
                continue
            yield p


def test_no_cjk_characters():
    offenders = [str(p.relative_to(ROOT)) for p in text_files() if CJK.search(p.read_text(errors="ignore"))]
    assert not offenders, f"CJK characters found in: {offenders}"


def test_no_personal_absolute_paths():
    offenders = []
    for p in list((ROOT / "scripts").rglob("*")) + list((ROOT / "notebooks").rglob("*")):
        if p.is_file() and p.suffix in TEXT_SUFFIXES:
            if "/Users/" in p.read_text(errors="ignore"):
                offenders.append(str(p.relative_to(ROOT)))
    assert not offenders, f"personal absolute paths in: {offenders}"


def test_no_large_files_tracked():
    big = [
        str(p.relative_to(ROOT))
        for p in ROOT.rglob("*")
        if p.is_file() and p.stat().st_size > MAX_BYTES
    ]
    assert not big, f"files over {MAX_BYTES} bytes: {big}"


@pytest.mark.parametrize("pattern", ["data/", "*.bin", "*.mtx", "outputs/"])
def test_gitignore_covers_data(pattern):
    ignored = (ROOT / ".gitignore").read_text()
    assert pattern in ignored


def test_required_files_exist():
    for name in ["README.md", "LICENSE", "Makefile", "requirements.txt", "environment.yml"]:
        assert (ROOT / name).is_file(), f"missing {name}"


def test_no_manuscript_or_submission_builder():
    """The manuscript and the docx builder that embeds it must never be committed."""
    banned = [
        p for p in ROOT.rglob("*")
        if p.is_file() and ("MANUSCRIPT" in p.name.upper() or "submission_docx" in p.name)
    ]
    assert not banned, f"manuscript artefacts inside repo: {[str(p.relative_to(ROOT)) for p in banned]}"
