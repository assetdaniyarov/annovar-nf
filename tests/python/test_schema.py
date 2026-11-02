"""nextflow_schema.json must describe exactly the parameters of nextflow.config."""

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def config_params() -> dict[str, str]:
    """Names and raw default values of the top-level `params { }` block."""
    text = (ROOT / "nextflow.config").read_text(encoding="utf-8")
    block = re.search(r"^params \{\n(.*?)^\}", text, re.S | re.M).group(1)
    return dict(re.findall(r"^    (\w+)\s*=\s*(.+?)\s*(?://.*)?$", block, re.M))


def schema_params() -> dict[str, dict]:
    schema = json.loads((ROOT / "nextflow_schema.json").read_text(encoding="utf-8"))
    out = {}
    for group in schema["$defs"].values():
        out.update(group["properties"])
    assert [ref["$ref"].split("/")[-1] for ref in schema["allOf"]] == list(schema["$defs"])
    return out


def test_schema_covers_config():
    # annovar_protocols is set in conf/annovar_protocols.config
    assert set(schema_params()) == set(config_params()) | {"annovar_protocols"}


def test_schema_defaults_match_config():
    config = config_params()
    for name, prop in schema_params().items():
        if "default" not in prop or name not in config:
            continue
        raw = config[name].strip("'\"")
        assert str(prop["default"]).lower() == raw.lower(), name


def test_one_version_everywhere():
    """The release version is the same in the manifest, the container tag, CITATION.cff and CHANGELOG.md."""
    config = (ROOT / "nextflow.config").read_text(encoding="utf-8")
    version = re.search(r"^\s+version\s*=\s*'([^']+)'", config, re.M).group(1)
    assert config_params()["container"].strip("'") == f"ghcr.io/assetdaniyarov/annovar-nf:{version}"
    assert f"version: {version}\n" in (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert re.search(r"^## \[([^\]]+)\]", changelog, re.M).group(1) == version
