"""Small shared I/O and process helpers."""
import csv
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
BIN = ROOT / "bin"


def executable(name):
    sibling = Path(sys.executable).parent / name
    return str(sibling) if sibling.is_file() else (shutil.which(name) or str(name))


def run_checked(command, **kwargs):
    env = os.environ.copy()
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env.get("PATH", "")
    env["PYTHONHASHSEED"] = "0"
    env.update(kwargs.pop("env", {}))
    return subprocess.run([str(value) for value in command], env=env, check=True, **kwargs)


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row[field] for field in fields} for row in rows)


def prefix_value(prefix):
    return "" if prefix == "." else prefix


def unique_file(directory, pattern):
    paths = sorted(Path(directory).glob(pattern))
    if len(paths) != 1:
        raise ValueError("Expected one {} in {}, found {}".format(pattern, directory, len(paths)))
    return paths[0]
