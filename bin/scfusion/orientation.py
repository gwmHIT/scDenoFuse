"""Gene-strand annotation and lightweight retained-side compatibility filtering."""
import csv
from collections import defaultdict
from pathlib import Path


REASON = "incompatible_gene_orientation"
STATUS_COLUMNS = ["OrientationStatus", "TranscriptPart1", "TranscriptPart2", "JunctionTopology"]


def transcript_part(retained_side, gene_strand):
    """Return the retained transcript portion, not the read's SAM strand.

    '-' keeps the genomic left side; '+' keeps the genomic right side.
    For a plus-strand gene left/right are 5prime/3prime; for minus, reversed.
    """
    if retained_side not in ("+", "-") or gene_strand not in ("+", "-"):
        return "unknown"
    return "3prime" if retained_side == gene_strand else "5prime"


def apply(statuses, mode):
    """After baseline3, remove known 5prime--5prime or 3prime--3prime pairs.

    Unknown directions (including cluster ties) are not grounds for removal.
    No read scan, evidence recount, annotation lookup or external label is used.
    """
    before = sum(row["SelectionStatus"] == "SELECTED" for row in statuses)
    for row in statuses:
        status = "DISABLED" if mode == "off" else "NOT_APPLIED"
        part1 = part2 = topology = ""
        if mode == "on" and row["SelectionStatus"] == "SELECTED":
            part1 = transcript_part(row.get("direction1"), row.get("strand1"))
            part2 = transcript_part(row.get("direction2"), row.get("strand2"))
            if "unknown" in (part1, part2):
                status, topology = "UNKNOWN", "unknown"
            else:
                topology = part1 + "--" + part2
                status = "PASS" if part1 != part2 else "FILTERED"
                if status == "FILTERED":
                    row["SelectionStatus"] = "FILTERED"
                    row["SelectionReason"] = REASON
        row.update(OrientationStatus=status, TranscriptPart1=part1,
                   TranscriptPart2=part2, JunctionTopology=topology)
    return {"OrientationFilterMode": mode,
            "BeforeOrientationSelected": before,
            "OrientationChecked": before if mode == "on" else 0,
            "OrientationUnknown": sum(row["OrientationStatus"] == "UNKNOWN" for row in statuses),
            "OrientationRemoved": sum(row["OrientationStatus"] == "FILTERED" for row in statuses)}


def gene_aliases(attributes):
    values = {}
    for item in attributes.split(";"):
        item = item.strip()
        if " " in item:
            key, value = item.split(None, 1)
            if key in ("gene_name", "gene_id"):
                values[key] = value.strip('"')
    names = set(values.values())
    if "gene_id" in values:
        names.add(values["gene_id"].split(".")[0])
    return names


class GeneStrands:
    def __init__(self):
        self.spans = defaultdict(set)
        self.exon_strands = defaultdict(set)

    def observe(self, fields):
        feature, strand = fields[2], fields[6]
        if feature not in ("gene", "exon") or strand not in ("+", "-"):
            return
        for name in gene_aliases(fields[8]):
            key = (name, fields[0])
            if feature == "gene":
                self.spans[key].add((int(fields[3]), int(fields[4]), strand))
            else:
                self.exon_strands[key].add(strand)

    def write(self, path):
        with Path(path).open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["Gene", "Chrom", "Start", "End", "Strand"])
            for (gene, chrom), spans in sorted(self.spans.items()):
                for start, end, strand in sorted(spans):
                    writer.writerow([gene, chrom, start, end, strand])
            for (gene, chrom), strands in sorted(self.exon_strands.items()):
                if (gene, chrom) not in self.spans:
                    for strand in sorted(strands):
                        writer.writerow([gene, chrom, 0, 0, strand])


def read_gene_strands(path):
    spans = defaultdict(list)
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            spans[(row["Gene"], row["Chrom"])].append(
                (int(row["Start"]), int(row["End"]), row["Strand"]))
    return spans


def resolve_strand(spans, gene, chrom, breakpoint):
    records = spans.get((gene, chrom), ())
    matched = {strand for start, end, strand in records if start <= breakpoint <= end}
    strands = matched or {strand for _, _, strand in records}
    return next(iter(strands)) if len(strands) == 1 else "?"


def annotate(rows, spans):
    for row in rows:
        for side in (1, 2):
            row["strand{}".format(side)] = resolve_strand(
                spans, row["Gene{}".format(side)], row["Chrom{}".format(side)],
                int(row["Breakpoint{}".format(side)]))
    return rows
