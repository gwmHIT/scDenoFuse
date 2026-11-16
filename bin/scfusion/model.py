"""Named-table interface to the original, sequential nine-parameter ZINB model."""
import ast
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil

from .common import executable, run_checked


RESULT_COLUMNS = [
    "FusionName", "CellSupport", "JunctionReadCount", "SpanningFragCount",
    "Position1", "Position2", "Pv23", "CandidateID",
    "direction1", "direction2", "strand1", "strand2",
]
INPUT_COLUMNS = [
    "Gene1", "Gene2", "CellSupport", "SplitReadCountsByCell",
    "SpanningSupportingCellCount", "SpanningReadCountsByCell", "Chrom1", "Chrom2",
    "Breakpoint1", "Breakpoint2", "ClusterIndex", "GC1", "GC2",
    "Gene1ExpressionBySplitCells", "Gene2ExpressionBySplitCells", "CandidateID",
    "SupportingCells", "PairSpanningCells", "PairSpanningCounts",
]


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_log_metadata(output, log_output, metadata):
    """Bind optional log scores to the exact legacy result bytes."""
    data = {"Version": 1, "AllresultSHA256": file_sha256(output),
            "LogScoresSHA256": file_sha256(log_output)}
    Path(metadata).write_text(json.dumps(data, indent=2) + "\n")


def iter_candidates(path):
    """Validate one row at a time, without retaining a second model input table."""
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(INPUT_COLUMNS) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing named columns: {', '.join(sorted(missing))}")
        seen = set()
        for row in reader:
            cid = row["CandidateID"]
            if cid in seen:
                raise ValueError(f"Duplicate CandidateID: {cid}")
            seen.add(cid)
            cells = ast.literal_eval(row["SupportingCells"])
            counts = ast.literal_eval(row["SplitReadCountsByCell"])
            expr1 = ast.literal_eval(row["Gene1ExpressionBySplitCells"])
            expr2 = ast.literal_eval(row["Gene2ExpressionBySplitCells"])
            if not (len(cells) == len(counts) == len(expr1) == len(expr2) == int(row["CellSupport"])):
                raise ValueError(f"{cid}: split cells, counts and expression vectors differ in length")
            if len(ast.literal_eval(row["PairSpanningCells"])) != len(ast.literal_eval(row["PairSpanningCounts"])):
                raise ValueError(f"{cid}: spanning cells and counts differ in length")
            yield row


def read_candidates(path):
    return list(iter_candidates(path))


def run(args):
    """Produce Allresult.txt; the caller subsequently runs the annotation filter."""
    outdir = Path(args.outdir)
    prefix = "" if args.prefix == "." else args.prefix
    candidates = outdir / "Candidates" / f"{prefix}ChiDist.tsv"
    background = outdir / "Candidates" / f"{prefix}ChiDist_background.tsv"
    output = outdir / f"{prefix}Allresult.txt"
    parameters = outdir / "Model" / f"{prefix}model_parameters.RData"
    log_output = outdir / "Model" / f"{prefix}pvalue_log.tsv"
    log_metadata = outdir / "Model" / f"{prefix}pvalue_log.metadata.json"
    parameters.parent.mkdir(parents=True, exist_ok=True)
    candidate_count = sum(1 for _ in iter_candidates(candidates))
    sum(1 for _ in iter_candidates(background))
    if not candidate_count:
        with output.open("w", newline="") as handle:
            csv.writer(handle, delimiter="\t", lineterminator="\n").writerow(RESULT_COLUMNS)
        # A previous nonempty run must not leave a sidecar for these empty bytes.
        log_output.unlink(missing_ok=True)
        log_metadata.unlink(missing_ok=True)
        return output
    rscript = os.environ.get("SCFUSION_RSCRIPT") or getattr(args, "rscript", None) or "Rscript"
    if str(rscript) == "Rscript":
        rscript = executable("Rscript")
    resolved = shutil.which(str(rscript))
    if resolved is None:
        raise FileNotFoundError("Rscript not found; activate the scFusion environment or set --rscript")
    script = Path(__file__).resolve().parents[2] / "bin" / "zinb_model.R"
    command = [resolved, str(script), "--candidates", str(candidates),
               "--background", str(background), "--cells", str(args.end - args.begin + 1),
                "--output", str(output), "--parameters", str(parameters),
                "--log-output", str(log_output)]
    run_checked(command)
    write_log_metadata(output, log_output, log_metadata)
    return output
