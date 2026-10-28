import os, sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


@pytest.fixture(autouse=True)
def cache(tmp_path, monkeypatch):
    """Point the measurement cache at tmp_path (env var reaches spawned workers too)."""
    import subselect
    d = str(tmp_path / "cache")
    monkeypatch.setenv("SUBSELECT_CACHE", d)
    monkeypatch.setattr(subselect, "CACHE", d)
    return d
