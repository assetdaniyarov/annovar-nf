"""Make the pipeline's bin/ scripts importable as modules in the tests."""

import importlib.util
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[2] / "bin"
sys.path.insert(0, str(BIN))        # scripts import their sibling modules (e.g. inheritance)


def load_script(name: str):
    """Import bin/<name>.py as a module (the scripts are not a package)."""
    path = BIN / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def detect():
    return load_script("detect_genome_build")


@pytest.fixture(scope="session")
def plan():
    return load_script("plan_chunks")


@pytest.fixture(scope="session")
def fmt():
    return load_script("format_annovar_table")


@pytest.fixture(scope="session")
def summarize():
    return load_script("summarize_annotations")


@pytest.fixture(scope="session")
def sqc():
    return load_script("sample_qc")


@pytest.fixture(scope="session")
def inherit():
    return load_script("inheritance")
