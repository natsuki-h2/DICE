"""Record reproducible provenance for DICE products.

The manifest uses relative paths and compact metadata for
large CASA tables.  Full SHA-256 hashes are practical for catalogues and
figures, but are opt-in for measurement sets and image directories.
"""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
from typing import Iterable, Mapping, Sequence


def _relative(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path.resolve())


def _directory_size(path: Path) -> int:
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    """Return the SHA-256 checksum of one regular file."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def describe_path(
    path: str | Path,
    *,
    root: str | Path,
    checksum: bool = False,
) -> dict:
    """Describe a file or directory without hashing large CASA tables by default."""

    item = Path(path)
    project_root = Path(root)
    record = {
        "path": _relative(item, project_root),
        "exists": item.exists(),
    }
    if not item.exists():
        return record
    stat = item.stat()
    record.update(
        {
            "kind": "directory" if item.is_dir() else "file",
            "size_bytes": _directory_size(item) if item.is_dir() else stat.st_size,
            "modified_utc": datetime.fromtimestamp(
                stat.st_mtime, timezone.utc
            ).isoformat(),
        }
    )
    if checksum and item.is_file():
        record["sha256"] = sha256_file(item)
    return record


def git_commit(root: str | Path) -> str | None:
    """Return the current commit, or ``None`` outside a Git checkout."""

    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=Path(root),
        capture_output=True,
        text=True,
        check=False,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def write_provenance(
    output: str | Path,
    *,
    root: str | Path,
    script: str | Path,
    command: Sequence[str],
    inputs: Iterable[str | Path] = (),
    outputs: Iterable[str | Path] = (),
    config: str | Path | None = None,
    parameters: Mapping | None = None,
    random_seeds: Sequence[int] = (),
    casa_version: str | None = None,
    checksum_inputs: bool = False,
) -> dict:
    """Write a JSON sidecar describing how a DICE product was generated."""

    project_root = Path(root)
    manifest = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit(project_root),
        "script": _relative(Path(script), project_root),
        "command": list(command),
        "config": _relative(Path(config), project_root) if config else None,
        "inputs": [
            describe_path(item, root=project_root, checksum=checksum_inputs)
            for item in inputs
        ],
        "outputs": [
            describe_path(item, root=project_root, checksum=True)
            for item in outputs
        ],
        "parameters": dict(parameters or {}),
        "random_seeds": [int(seed) for seed in random_seeds],
        "software": {
            "python": platform.python_version(),
            "casa": casa_version,
        },
    }
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def write_script_provenance(
    output: str | Path,
    *,
    root: str | Path,
    script: str | Path,
    inputs: Iterable[str | Path] = (),
    outputs: Iterable[str | Path] = (),
    config: str | Path | None = None,
    parameters: Mapping | None = None,
    random_seeds: Sequence[int] = (),
    casa_version: str | None = None,
) -> dict:
    """Write a sidecar for a directly executed workflow script.

    This compact wrapper standardizes ``sys.argv`` recording in plotting
    entry points while retaining the full :func:`write_provenance` schema.
    """

    return write_provenance(
        output,
        root=root,
        script=script,
        command=[sys.executable, *sys.argv],
        inputs=inputs,
        outputs=outputs,
        config=config,
        parameters=parameters,
        random_seeds=random_seeds,
        casa_version=casa_version,
    )
