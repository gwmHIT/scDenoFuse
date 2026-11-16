"""Exact, label-blind matching against a frozen recurrent-breakpoint resource.

Only the Arriba ``## recurrent breakpoints`` section and unsigned, single-base
coordinate pairs are used. This is not a reimplementation of Arriba's filters.
"""
import gzip
import hashlib
from pathlib import Path

from .common import ROOT


DEFAULT_BLACKLIST = ROOT / "data" / "blacklist_hg38_GRCh38_v2.5.1.tsv.gz"
SECTION = "## recurrent breakpoints"
REASON = "recurrent_background_breakpoint"
STATUS_COLUMNS = ["BackgroundStatus", "BackgroundMatchedPair", "BackgroundSourceLine"]


def coordinate(value):
    """Normalize chr prefixes without converting assemblies or shifting bases."""
    chrom, pos = value.rsplit(":", 1)
    if chrom.startswith(("+", "-")) or not pos.isdecimal() or int(pos) < 1:
        raise ValueError(f"Expected an unsigned 1-based breakpoint, got {value!r}")
    return f"{chrom.removeprefix('chr')}:{int(pos)}"


def validate_resource(path, genome_dir):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Background blacklist not found: {path}; provide --background-blacklist or use --background-filter off")
    if path.resolve() == DEFAULT_BLACKLIST.resolve():
        # The bundled coordinate resource is GRCh38-specific. Do not silently
        # apply it to hg19/mouse or write an index in a shared reference folder.
        index = Path(genome_dir) / "ref.fa.fai"
        with index.open() as handle:
            lengths = {fields[0].removeprefix('chr'): int(fields[1])
                       for fields in (line.split('\t') for line in handle)}
        expected = {"1": 248956422, "2": 242193529, "X": 156040895}
        if any(lengths.get(chrom) != length for chrom, length in expected.items()):
            raise ValueError("Bundled background blacklist requires GRCh38 (chr1/chr2/chrX lengths); provide a matching --background-blacklist or use --background-filter off")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def recurrent_matches(rows, path):
    """Stream the resource once; memory scales with candidates, not database size.

    Return CandidateID -> (resource pair, 1-based source line). The first matching
    record is retained for audit. Names, labels, directions and read counts do
    not enter this exact unordered coordinate-pair match.
    """
    wanted = {}
    for row in rows:
        a, b = coordinate(row["Position1"]), coordinate(row["Position2"])
        for key in {(a, b), (b, a)}:
            wanted.setdefault(key, []).append(row["CandidateID"])
    if not wanted:
        return {}
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    found, active, section_seen = {}, False, False
    with opener(path, "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("##"):
                active = line.strip().lower() == SECTION
                section_seen |= active
                continue
            if not active or line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"{path}:{line_number}: expected two blacklist columns")
            # Ranges, gene names, strand-qualified entries and Arriba keywords
            # are deliberately NOT interpreted by this narrow mechanism.
            key = tuple(field.removeprefix("chr") for field in fields)
            for candidate_id in wanted.get(key, ()):
                found.setdefault(candidate_id, ("|".join(fields), line_number))
    if not section_seen:
        raise ValueError(f"{path}: missing {SECTION!r} section")
    return found


def apply(statuses, mode, path):
    """Last filtering stage; never influence bad groups or the Pv23 cutoff."""
    eligible = [row for row in statuses if row["SelectionStatus"] == "SELECTED"]
    matches = recurrent_matches(eligible, path) if mode == "on" else {}
    for row in statuses:
        status, pair, line = ("DISABLED" if mode == "off" else "NOT_APPLIED"), "", ""
        if mode == "on" and row["SelectionStatus"] == "SELECTED":
            status = "PASS"
            if row["CandidateID"] in matches:
                status = "FILTERED"
                pair, line = matches[row["CandidateID"]]
                row["SelectionStatus"] = "FILTERED"
                row["SelectionReason"] = REASON
        row.update(BackgroundStatus=status, BackgroundMatchedPair=pair, BackgroundSourceLine=line)
    return {"BackgroundFilterMode": mode,
            "BackgroundBlacklist": str(Path(path).resolve()) if mode == "on" else "NA",
            "BackgroundBlacklistSHA256": sha256(path) if mode == "on" else "NA",
            "BackgroundMatching": "exact_unordered_pair" if mode == "on" else "NA",
            "BeforeBackgroundSelected": len(eligible),
            "BackgroundChecked": len(eligible) if mode == "on" else 0,
            "BackgroundRemoved": sum(row["BackgroundStatus"] == "FILTERED" for row in statuses)}
