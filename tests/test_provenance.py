import json
from pathlib import Path

from bliss_alma.provenance import (
    describe_path,
    sha256_file,
    write_provenance,
    write_script_provenance,
)


def test_describe_and_hash_file(tmp_path):
    source = tmp_path / "input.csv"
    source.write_text("x,y\n1,2\n")
    record = describe_path(source, root=tmp_path, checksum=True)
    assert record["path"] == "input.csv"
    assert record["size_bytes"] == source.stat().st_size
    assert record["sha256"] == sha256_file(source)


def test_write_provenance_uses_relative_paths(tmp_path):
    script = tmp_path / "scripts" / "plot.py"
    product = tmp_path / "data" / "project" / "figures" / "figure.pdf"
    script.parent.mkdir()
    product.parent.mkdir(parents=True)
    script.write_text("pass\n")
    product.write_bytes(b"figure")
    sidecar = product.with_suffix(".provenance.json")

    write_provenance(
        sidecar,
        root=tmp_path,
        script=script,
        command=["python", str(script)],
        outputs=[product],
        random_seeds=[10, 11],
    )
    manifest = json.loads(sidecar.read_text())
    assert manifest["script"] == "scripts/plot.py"
    assert manifest["outputs"][0]["path"] == "data/project/figures/figure.pdf"
    assert manifest["random_seeds"] == [10, 11]


def test_write_script_provenance_records_current_command(tmp_path):
    source = tmp_path / "input.csv"
    source.write_text("x\n1\n")
    output = tmp_path / "figure.png"
    output.write_bytes(b"png")
    destination = tmp_path / "figure.provenance.json"
    payload = write_script_provenance(
        destination,
        root=tmp_path,
        script=tmp_path / "plot.py",
        inputs=[source],
        outputs=[output],
    )
    assert payload["command"]
    assert payload["inputs"][0]["path"] == "input.csv"
    assert payload["outputs"][0]["path"] == "figure.png"
