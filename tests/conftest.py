"""Shared pytest fixtures for V2F Reporter tests."""

import os
import shutil
import sys
import tempfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))


@pytest.fixture
def repo_root():
    return REPO_ROOT


@pytest.fixture
def isolated_workdir(tmp_path, monkeypatch):
    """Run a test with PROJECT_DIR pointed at a temp dir.

    Used by tests that exercise the demo load/reset endpoints — they write
    to data/, reports/, etc. We don't want them clobbering the developer's
    real work tree.
    """
    # Mirror demo/ into the temp dir so dashboard/app.py can find it
    src = REPO_ROOT / "demo"
    if src.exists():
        shutil.copytree(src, tmp_path / "demo")
    monkeypatch.chdir(tmp_path)
    return tmp_path
