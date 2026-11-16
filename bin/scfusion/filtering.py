"""Two-of-three threshold, annotation/R2/R3, orientation, then background."""
from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from .model import RESULT_COLUMNS, file_sha256
from .common import read_tsv as read_table, write_tsv as write_table
from .baseline3 import Gene, Screens, STATUS_COLUMNS, empty_status
from . import background, orientation


def read_cell_counts(path, candidate_ids=None):
    counts = defaultdict(dict)
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            cid = row["CandidateID"]
            if candidate_ids is not None and cid not in candidate_ids:
                continue
            cell = int(row["Cell"])
            if cell in counts[cid]:
                raise ValueError(f"Repeated candidate/cell count: {cid}, {cell}")
            split, span = int(row["SplitCount"]), int(row["SpanningCount"])
            if split < 0 or span < 0:
                raise ValueError(f"Negative candidate/cell count: {cid}, {cell}")
            if split > 0:
                counts[cid][cell] = (split, span)
    return counts


def bad_statistics(row, cell_counts, total_cells):
    # cell_counts maps Cell -> (SplitCount, SpanningCount).
    positive = [value for value in cell_counts.values() if value[0] > 0]
    c = len(positive)
    junctions = sum(split for split, _ in positive)
    if c != int(row["CellSupport"]) or junctions != int(row["JunctionReadCount"]):
        raise ValueError(f"{row['CandidateID']}: model summary differs from per-cell split counts")
    a = junctions / c if c else 0.0
    r = sum(span > 0 for _, span in positive) / c if c else 0.0
    k = max(2, math.floor(0.05 * total_cells))
    flags = (c < k, a < 1.25, r < 0.50)
    return {"C": c, "A": a, "R": r, "K": k, "LowCellSupport": flags[0],
            "LowAverageSplit": flags[1], "LowSpanningCellFraction": flags[2],
            "Bad": sum(flags) >= 2}


@dataclass(frozen=True)
class Threshold:
    cutoff: float | None
    exponent: float | None
    estimated_fdr: float | None
    good_count: int
    bad_count: int
    log_cutoff: float | None = None

    @property
    def valid(self):
        return self.cutoff is not None or self.log_cutoff is not None

    def accepts(self, pvalue, log_pvalue=None):
        if self.log_cutoff is not None:
            return log_pvalue is not None and log_pvalue <= self.log_cutoff
        return self.cutoff is not None and pvalue <= self.cutoff


def dynamic_threshold(pvalues, bad_flags, target=0.05, log_pvalues=None):
    if len(pvalues) != len(bad_flags):
        raise ValueError("P-values and bad flags differ in length")
    if any(not math.isfinite(p) or p < 0 or p > 1 for p in pvalues):
        raise ValueError("Pv23 must be a finite probability in [0, 1]")
    good = sorted(p for p, bad in zip(pvalues, bad_flags) if not bad)
    bad = sorted(p for p, bad in zip(pvalues, bad_flags) if bad)
    from bisect import bisect_left
    for level in range(1000):
        cutoff = 10.0 ** (-level - 0.5)
        if cutoff == 0.0:
            break  # Underflow is not a successful threshold.
        n_good, n_bad = bisect_left(good, cutoff), bisect_left(bad, cutoff)
        selected = n_good + n_bad
        if selected == 0:
            continue
        estimate = n_bad / max(1, len(bad)) * len(pvalues) / selected
        if estimate < target:
            return Threshold(cutoff, level + 0.5, estimate, len(good), len(bad))
    if log_pvalues is None:
        return Threshold(None, None, None, len(good), len(bad))
    if len(log_pvalues) != len(pvalues):
        raise ValueError("Log P-values and probabilities differ in length")
    if any(math.isnan(p) or p > 0 for p in log_pvalues):
        raise ValueError("LogPv23 must be nonpositive and not NaN")
    # Keep all legacy-representable probabilities exact. Only zero probabilities
    # require the sidecar's lost tail information.
    logs = [math.log(p) if p > 0 else lp for p, lp in zip(pvalues, log_pvalues)]
    log_good = sorted(p for p, flag in zip(logs, bad_flags) if not flag)
    log_bad = sorted(p for p, flag in zip(logs, bad_flags) if flag)
    # The half-integer grid changes membership only when it crosses a score.
    # Visit these changes directly instead of allocating millions of grid steps.
    levels = {0}
    for lp in logs:
        if math.isfinite(lp):
            crossing = max(0, math.ceil(-lp / math.log(10) - .5))
            # Check adjacent grid points too, including floating-point ties.
            levels.update((max(0, crossing - 1), crossing, crossing + 1))
    for level in sorted(levels):
        exponent = level + .5
        log_cutoff = -exponent * math.log(10)
        n_good = bisect_left(log_good, log_cutoff)
        n_bad = bisect_left(log_bad, log_cutoff)
        selected = n_good + n_bad
        if selected:
            estimate = n_bad / max(1, len(bad)) * len(pvalues) / selected
            if estimate < target:
                return Threshold(None, exponent, estimate, len(good), len(bad), log_cutoff)
    return Threshold(None, None, None, len(good), len(bad))


def read_log_pvalues(outdir, prefix, rows):
    path = outdir / "Model" / f"{prefix}pvalue_log.tsv"
    if not path.is_file():
        return None
    metadata_path = path.with_name(f"{prefix}pvalue_log.metadata.json")
    metadata = json.loads(metadata_path.read_text())
    if metadata.get("Version") != 1 or metadata.get("AllresultSHA256") != file_sha256(outdir / f"{prefix}Allresult.txt") or metadata.get("LogScoresSHA256") != file_sha256(path):
        raise ValueError("Log score sidecar does not match the saved result bytes")
    values = {}
    for record in read_table(path):
        cid = record["CandidateID"]
        if cid in values:
            raise ValueError(f"Duplicate log-score CandidateID: {cid}")
        values[cid] = (float(record["Pv23"]), float(record["LogPv23"]))
    if set(values) != {r["CandidateID"] for r in rows}:
        raise ValueError("Log-score CandidateIDs differ from saved results")
    if any(values[r["CandidateID"]][0] != float(r["Pv23"]) for r in rows):
        raise ValueError("Log-score probabilities differ from saved results")
    return [values[r["CandidateID"]][1] for r in rows]


def attributes(text):
    values = {}
    for item in text.split(";"):
        item = item.strip()
        if " " in item:
            key, value = item.split(" ", 1)
            values[key] = value.strip('"')
    return values


def read_annotation(path, gene_records=None):
    """Keep individual gene features: do not union spans over separate loci."""
    types, exons, genes = defaultdict(set), defaultdict(list), defaultdict(list)
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            chrom, _, feature, start, end, _, strand, _, text = fields
            if feature not in ("gene", "exon"):
                continue
            info = attributes(text)
            names = {value for value in (info.get("gene_name"), info.get("gene_id")) if value}
            kind = info.get("gene_type", info.get("gene_biotype", ""))
            interval = (chrom, int(start), int(end), strand)
            if feature == "gene" and gene_records is not None:
                gene_records.append(Gene.from_fields(fields, kind))
            for gene in names:
                if kind:
                    types[gene].add(kind)
                (genes if feature == "gene" else exons)[gene].append(interval)
    return types, exons, genes


def position(text):
    chrom, coordinate = text.rsplit(":", 1)
    return chrom, int(coordinate)


def overlaps_exon(endpoint, gene, exons, window=100):
    chrom, coordinate = position(endpoint)
    return any(chrom == c and start - window <= coordinate <= end + window
               for c, start, end, _ in exons.get(gene, ()))


def annotation_reasons(row, annotation, exon_window=100):
    types, exons, genes = annotation
    gene1, gene2 = row["FusionName"].split("--", 1)
    t1, t2 = types[gene1], types[gene2]
    reasons = []
    if t1 and t2 and "protein_coding" not in t1 and "protein_coding" not in t2:
        reasons.append("both_ends_noncoding")
    if any("pseudogene" in kind for kind in t1 | t2):
        reasons.append("pseudogene")
    if not overlaps_exon(row["Position1"], gene1, exons, exon_window) or not overlaps_exon(row["Position2"], gene2, exons, exon_window):
        reasons.append("outside_own_exon_window")
    chrom1, bp1 = position(row["Position1"])
    chrom2, bp2 = position(row["Position2"])
    if chrom1 == chrom2 and any(chrom == chrom1 and start <= min(bp1, bp2) and max(bp1, bp2) <= end
                               for gene in (gene1, gene2)
                               for chrom, start, end, _ in genes.get(gene, ())):
        reasons.append("partner_gene_covers_both_breakpoints")
    return reasons


def run(args):
    background_mode = getattr(args, "background_filter", "on")
    blacklist = getattr(args, "background_blacklist", background.DEFAULT_BLACKLIST)
    if background_mode == "on":
        background.validate_resource(blacklist, args.genome_dir)
    outdir = Path(args.outdir)
    prefix = "" if args.prefix == "." else args.prefix
    rows = read_table(outdir / f"{prefix}Allresult.txt")
    orientation_columns = ["direction1", "direction2", "strand1", "strand2"]
    candidate_path = outdir / "Candidates" / f"{prefix}ChiDist.tsv"
    direction_columns = ["DirectionCountPP", "DirectionCountPM", "DirectionCountMP",
                         "DirectionCountMM", "DirectionAmbiguousCount", "DirectionStatus"]
    direction_metadata = {}
    if candidate_path.is_file():
        candidate_ids = {row["CandidateID"] for row in rows}
        with candidate_path.open(newline="") as handle:
            for candidate in csv.DictReader(handle, delimiter="\t"):
                if candidate["CandidateID"] in candidate_ids and "DirectionStatus" in candidate:
                    direction_metadata[candidate["CandidateID"]] = {
                        key: candidate.get(key, "?") for key in direction_columns + orientation_columns}
    for row in rows:
        metadata = direction_metadata.get(row["CandidateID"], {})
        for field in orientation_columns:
            row.setdefault(field, metadata.get(field, "?"))
    counts = read_cell_counts(outdir / f"{prefix}candidate_cell_counts.tsv", {row["CandidateID"] for row in rows})
    n = args.end - args.begin + 1
    stats = [bad_statistics(row, counts[row["CandidateID"]], n) for row in rows]
    target = getattr(args, "fdr_cutoff", 0.05)
    threshold = dynamic_threshold([float(row["Pv23"]) for row in rows], [s["Bad"] for s in stats], target)
    log_pvalues = None
    if not threshold.valid:
        log_pvalues = read_log_pvalues(outdir, prefix, rows)
        if log_pvalues is not None:
            threshold = dynamic_threshold([float(row["Pv23"]) for row in rows],
                                          [s["Bad"] for s in stats], target, log_pvalues)
    gtf = getattr(args, "gtf", None) or Path(args.genome_dir) / "ref_annot.gtf.added"
    mode = getattr(args, "baseline3", "on")
    gene_records = []
    annotation = read_annotation(gtf, gene_records if mode == "on" else None)
    window = getattr(args, "exon_window", 100)
    genome_dir = getattr(args, "genome_dir", Path(gtf).parent)
    statuses, selected = [], []
    with Screens(Path(genome_dir) / "ref.fa", gene_records) as screens:
        for index, (row, stat) in enumerate(zip(rows, stats)):
            significant = threshold.accepts(float(row["Pv23"]),
                                            log_pvalues[index] if log_pvalues is not None else None)
            reasons = ([] if significant else ["model_not_significant"]) + annotation_reasons(row, annotation, window)
            baseline1_selected = not reasons
            details = dict(empty_status("DISABLED" if mode == "off" else "NOT_APPLIED"),
                           Baseline1Selected="YES" if baseline1_selected else "NO", Baseline3Rules="none")
            if mode == "on" and baseline1_selected:
                extra_details, extra_reasons = screens.evaluate(row)
                details.update(extra_details)
                reasons.extend(extra_reasons)
            status = dict(row, **stat, **details, ModelSignificant=significant,
                          SelectionStatus="FILTERED" if reasons else "SELECTED",
                          SelectionReason=";".join(reasons) if reasons else "pass")
            metadata = direction_metadata.get(row["CandidateID"], {})
            status.update({key: metadata.get(key, "UNKNOWN" if key == "DirectionStatus" else 0)
                           for key in direction_columns})
            for flag in ("LowCellSupport", "LowAverageSplit", "LowSpanningCellFraction", "Bad", "ModelSignificant"):
                status[flag] = "YES" if status[flag] else "NO"
            statuses.append(status)
            if not reasons:
                selected.append(row)
    baseline3_removed = sum(row["Baseline1Selected"] == "YES" and row["SelectionStatus"] == "FILTERED"
                            for row in statuses)
    orientation_summary = orientation.apply(statuses, getattr(args, "orientation_filter", "on"))
    background_summary = background.apply(statuses, background_mode, blacklist)
    selected = [row for row in statuses if row["SelectionStatus"] == "SELECTED"]
    write_table(outdir / f"{prefix}Result.full.txt", RESULT_COLUMNS, selected)
    extra = ["C", "A", "R", "K", "LowCellSupport", "LowAverageSplit", "LowSpanningCellFraction",
             "Bad", "ModelSignificant", "SelectionStatus", "SelectionReason"] + direction_columns + STATUS_COLUMNS + orientation.STATUS_COLUMNS + background.STATUS_COLUMNS
    if getattr(args, "diagnostics", True):
        write_table(outdir / f"{prefix}Allresult.selection_status.tsv", RESULT_COLUMNS + extra, statuses)
    summary = {"InputCandidates": len(rows), "TotalCells": n, "MinimumCellsForGood": max(2, math.floor(.05*n)),
               "GoodCandidates": threshold.good_count, "BadCandidates": threshold.bad_count,
               "TargetFDR": target, "ThresholdStatus": "VALID" if threshold.valid else "NO_VALID_THRESHOLD",
               "Pv23Cutoff": threshold.cutoff, "MinusLog10Pv23Cutoff": threshold.exponent,
               "EstimatedFDR": threshold.estimated_fdr,
               "ModelSignificant": sum(row["ModelSignificant"] == "YES" for row in statuses),
               "Baseline1Selected": sum(row["Baseline1Selected"] == "YES" for row in statuses),
               "Baseline3Mode": mode,
               "Baseline3R2Hits": sum(row["R2Status"] == "FILTERED" for row in statuses),
               "Baseline3R3Hits": sum(row["R3Status"] == "FILTERED" for row in statuses),
               "Baseline3Removed": baseline3_removed,
               **orientation_summary,
               **background_summary,
               "Selected": len(selected)}
    if log_pvalues is not None:
        summary["ThresholdProbabilitySpace"] = "LOG_RESCUE" if threshold.valid else "LOG_NO_VALID_THRESHOLD"
        summary["NaturalLogPv23Cutoff"] = threshold.log_cutoff
    write_table(outdir / f"{prefix}threshold_summary.tsv", ["Metric", "Value"],
                [{"Metric": key, "Value": value if value is not None else "NA"} for key, value in summary.items()])
    return outdir / f"{prefix}Result.full.txt"
