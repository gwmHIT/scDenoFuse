"""Export final, row-preserving public tables without changing native results.

The native result's direction fields describe retained genomic sides.  Public
directions describe transcript sides: equal to gene strand is downstream (d),
opposite is upstream (u).  An explicit d/u row is presented as u/d by swapping
both partners.  Candidate identifiers remain in the export manifest only.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile


ABRIDGED_COLUMNS = (
    "FusionName", "Position1", "Position2", "strands", "directions",
    "CellSupport", "JunctionReadCount", "SpanningFragCount", "Pvalue",
)
FULL_COLUMNS = ABRIDGED_COLUMNS + (
    "SupportingCells", "SplitReadCountsByCell", "SpanningCells", "SpanningReadCountsByCell",
)
NATIVE_COLUMNS = (
    "FusionName", "CellSupport", "JunctionReadCount", "SpanningFragCount",
    "Position1", "Position2", "Pv23", "CandidateID",
    "direction1", "direction2", "strand1", "strand2",
)
COUNT_COLUMNS = ("CandidateID", "Cell", "SplitCount", "SpanningCount")
MANIFEST_NAME = "result_export_manifest.json"
FORMAT_VERSION = 1
EXPORT_VERSION = 1
_RULES = {
    "abridged_columns": ABRIDGED_COLUMNS,
    "full_columns": FULL_COLUMNS,
    "direction": "retained side equals gene strand: d; opposite: u; unknown: NA",
    "swap": "explicit d/u only; swap fusion partners, positions, strands, directions",
    "rows": "native order and count, no aggregation",
    "score": "Pvalue is the unchanged native Pv23 string",
    "cells": "positive counts only, numeric cell order, all spanning cells included",
}
COMPATIBILITY_SHA256 = hashlib.sha256(
    json.dumps(_RULES, sort_keys=True, separators=(",", ":")).encode("utf-8")
).hexdigest()
_INTEGER = re.compile(r"[0-9]+\Z")


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _integer(value, context):
    if not isinstance(value, str) or not _INTEGER.fullmatch(value):
        raise ValueError(f"{context}: expected a nonnegative integer, got {value!r}")
    return int(value)


def _reader(handle, columns, path):
    reader = csv.DictReader(handle, delimiter="\t")
    names = reader.fieldnames or []
    if len(names) != len(set(names)):
        raise ValueError(f"{path}: duplicate column names")
    missing = set(columns) - set(names)
    if missing:
        raise ValueError(f"{path}: missing named columns: {', '.join(sorted(missing))}")
    return reader


def _check_row(row, path, number):
    if None in row or any(value is None for value in row.values()):
        raise ValueError(f"{path}: malformed TSV row {number}")


def _snapshot(path):
    value = path.stat()
    return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns


def _load_native(work):
    """Read selected counts only, checking source stability and streaming hashes."""
    result_path = work / "Result.full.txt"
    counts_path = work / "candidate_cell_counts.tsv"
    before = {path: _snapshot(path) for path in (result_path, counts_path)}
    rows = []
    by_candidate = {}
    with result_path.open(encoding="utf-8", newline="") as handle:
        for number, row in enumerate(_reader(handle, NATIVE_COLUMNS, result_path), 2):
            _check_row(row, result_path, number)
            cid = row["CandidateID"]
            if not cid or cid in by_candidate:
                raise ValueError(f"{result_path}: empty or duplicate CandidateID: {cid!r}")
            if len(row["FusionName"].split("--")) != 2 or not all(row["FusionName"].split("--")):
                raise ValueError(f"{cid}: expected exactly two fusion partners in FusionName")
            for column in ("CellSupport", "JunctionReadCount", "SpanningFragCount"):
                _integer(row[column], f"{cid}: {column}")
            rows.append(row)
            by_candidate[cid] = {}
    with counts_path.open(encoding="utf-8", newline="") as handle:
        for number, row in enumerate(_reader(handle, COUNT_COLUMNS, counts_path), 2):
            # Unselected counts can be very large; retain and parse selected IDs only.
            if row.get("CandidateID") not in by_candidate:
                continue
            _check_row(row, counts_path, number)
            cid = row["CandidateID"]
            cell = _integer(row["Cell"], f"{cid}: Cell")
            if cell in by_candidate[cid]:
                raise ValueError(f"Repeated candidate/cell count: {cid}, {cell}")
            split = _integer(row["SplitCount"], f"{cid}, cell {cell}: SplitCount")
            span = _integer(row["SpanningCount"], f"{cid}, cell {cell}: SpanningCount")
            by_candidate[cid][cell] = split, span
    for row in rows:
        cid = row["CandidateID"]
        values = by_candidate[cid].values()
        positive_splits = [split for split, _ in values if split > 0]
        if (len(positive_splits) != int(row["CellSupport"])
                or sum(positive_splits) != int(row["JunctionReadCount"])):
            raise ValueError(f"{cid}: native split totals differ from per-cell split counts")
        if sum(span for _, span in values) != int(row["SpanningFragCount"]):
            raise ValueError(f"{cid}: native SpanningFragCount differs from all per-cell spanning counts")
    hashes = {path.name: _sha256(path) for path in (result_path, counts_path)}
    for path, old in before.items():
        if _snapshot(path) != old:
            raise ValueError(f"Native export input changed while being read: {path}")
    return rows, by_candidate, hashes, before


def _public_rows(rows, counts):
    public, mapping = [], []
    for number, row in enumerate(rows, 1):
        partners = row["FusionName"].split("--")
        positions = [row["Position1"], row["Position2"]]
        strands = [row[f"strand{index}"] for index in (1, 2)]
        directions = []
        for index, strand in enumerate(strands, 1):
            side = row[f"direction{index}"]
            directions.append("NA" if side not in ("+", "-") or strand not in ("+", "-")
                              else "d" if side == strand else "u")
        strands = [strand if strand in ("+", "-") else "NA" for strand in strands]
        swapped = directions == ["d", "u"]
        if swapped:
            partners.reverse()
            positions.reverse()
            strands.reverse()
            directions.reverse()
        item = dict(FusionName="--".join(partners), Position1=positions[0], Position2=positions[1],
                    strands="/".join(strands), directions="/".join(directions), Pvalue=row["Pv23"])
        for column in ("CellSupport", "JunctionReadCount", "SpanningFragCount"):
            item[column] = row[column]
        cell_counts = sorted(counts[row["CandidateID"]].items())
        split_cells = [(cell, split) for cell, (split, _) in cell_counts if split > 0]
        span_cells = [(cell, span) for cell, (_, span) in cell_counts if span > 0]
        for cells, cell_column, count_column in (
                (split_cells, "SupportingCells", "SplitReadCountsByCell"),
                (span_cells, "SpanningCells", "SpanningReadCountsByCell")):
            item[cell_column] = json.dumps([cell for cell, _ in cells])
            item[count_column] = json.dumps([count for _, count in cells])
        public.append(item)
        mapping.append({"native_row": number, "public_row": number,
                        "CandidateID": row["CandidateID"], "swapped": swapped})
    return public, mapping


class _HashWriter:
    def __init__(self):
        self.digest = hashlib.sha256()

    def write(self, value):
        self.digest.update(value.encode("utf-8"))
        return len(value)


def _write_table(handle, columns, rows):
    writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n",
                            extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)


def _table_hash(columns, rows):
    handle = _HashWriter()
    _write_table(handle, columns, rows)
    return handle.digest.hexdigest()


def _manifest(rows, mapping, native_hashes):
    return {
        "format_version": FORMAT_VERSION,
        "export_version": EXPORT_VERSION,
        "compatibility_sha256": COMPATIBILITY_SHA256,
        "rows": len(rows),
        "row_mapping": mapping,
        "hashes": {
            "native": native_hashes,
            "public": {"Result.abridged.txt": _table_hash(ABRIDGED_COLUMNS, rows),
                       "Result.full.txt": _table_hash(FULL_COLUMNS, rows)},
        },
    }


def _destinations(work, output, create=False):
    """Reject symlink destinations and accidental replacement of native tables."""
    if work.resolve() == output.resolve():
        raise ValueError("Public output directory must differ from the native work directory")
    directories = (output, output / "logs")
    targets = (output / "Result.abridged.txt", output / "Result.full.txt",
               output / "logs" / MANIFEST_NAME)
    for directory in directories:
        if directory.is_symlink():
            raise ValueError(f"Refusing public symlink directory: {directory}")
        if directory.exists() and not directory.is_dir():
            raise ValueError(f"Public output directory is not a directory: {directory}")
    native_paths = {(work / name).resolve() for name in ("Result.full.txt", "candidate_cell_counts.tsv")}
    for target in targets:
        if target.is_symlink():
            raise ValueError(f"Refusing public symlink target: {target}")
        if target.exists() and not target.is_file():
            raise ValueError(f"Public target is not a regular file: {target}")
        if target.resolve() in native_paths:
            raise ValueError(f"Public target would overwrite native input: {target}")
    if create:
        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)
    return targets


def export_results(work: Path, output: Path) -> dict:
    """Atomically publish two TSVs and commit their deterministic manifest last.

    Files are replaced individually.  If interrupted between replacements, the
    manifest cannot validate the mixed export, so cleanup must call validate_export.
    No native files are created, modified, or removed.
    """
    work, output = Path(work), Path(output)
    targets = _destinations(work, output)
    native, counts, hashes, snapshots = _load_native(work)
    rows, mapping = _public_rows(native, counts)
    manifest = _manifest(rows, mapping, hashes)
    _destinations(work, output, create=True)
    staged = []
    try:
        for target, columns in zip(targets[:2], (ABRIDGED_COLUMNS, FULL_COLUMNS)):
            with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                             dir=target.parent, prefix=f".{target.name}.",
                                             suffix=".tmp", delete=False) as handle:
                staged.append(Path(handle.name))
                _write_table(handle, columns, rows)
                handle.flush()
                os.fsync(handle.fileno())
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                         dir=targets[2].parent, prefix=".result_export_manifest.",
                                         suffix=".tmp", delete=False) as handle:
            staged.append(Path(handle.name))
            json.dump(manifest, handle, sort_keys=True, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        for path, old in snapshots.items():
            if _snapshot(path) != old:
                raise ValueError(f"Native export input changed before publication: {path}")
        _destinations(work, output)
        for temporary, target in zip(staged, targets):
            os.replace(temporary, target)
    finally:
        for temporary in staged:
            temporary.unlink(missing_ok=True)
    return manifest


def validate_export(work: Path, output: Path) -> dict:
    """Return provenance only if native inputs, public bytes and mapping agree.

    Native counts are revalidated and the expected public table hashes are
    reconstructed.  Merely updating a recorded hash cannot bless a changed table.
    """
    work, output = Path(work), Path(output)
    targets = _destinations(work, output)
    try:
        with targets[2].open(encoding="utf-8") as handle:
            manifest = json.load(handle)
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"Missing or invalid export manifest: {targets[2]}") from error
    if not isinstance(manifest, dict) or manifest.get("compatibility_sha256") != COMPATIBILITY_SHA256:
        raise ValueError("Export manifest has incompatible presentation rules")
    try:
        native, counts, hashes, snapshots = _load_native(work)
        rows, mapping = _public_rows(native, counts)
        expected = _manifest(rows, mapping, hashes)
        if manifest != expected:
            raise ValueError("Export manifest is stale or differs from native results/counts")
        for target in targets[:2]:
            if _sha256(target) != expected["hashes"]["public"][target.name]:
                raise ValueError(f"Public export is stale or has been modified: {target}")
        for path, old in snapshots.items():
            if _snapshot(path) != old:
                raise ValueError(f"Native export input changed during validation: {path}")
    except OSError as error:
        raise ValueError(f"Export validation could not read a required file: {error}") from error
    return manifest
