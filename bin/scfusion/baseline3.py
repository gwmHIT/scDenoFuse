"""R2 reference-side competition and R3 lncRNA readthrough-like topology.

Uses candidate summaries, reference sequence and gene annotation only.
No read names, tracing, external support labels or model refitting.
"""
from bisect import bisect_left, bisect_right
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from .orientation import gene_aliases


WINDOW = 50
MIN_MATCHES = 40
MIN_IDENTITY = 0.90
MIN_SPLIT_SPANNING_RATIO = 5.0
R2_REASON = "baseline3_R2_one_sided_homology"
R3_REASON = "baseline3_R3_lncRNA_readthrough"
STATUS_COLUMNS = [
    "Baseline1Selected", "Baseline3Rules", "R2Status", "R2SplitSpanningRatio",
    "R2LeftMatches", "R2LeftLength", "R2LeftIdentity",
    "R2RightMatches", "R2RightLength", "R2RightIdentity",
    "R3Status", "R3InterveningProteinCodingGenes",
]


def anchored_prefix(a, b):
    """Global at the breakpoint, free distal end; same scoring as legacy R2.

Match +2, mismatch -2, gap -3. N and gaps count against identity.
Choose the highest-scoring qualifying distal endpoint if one exists;
otherwise report the highest-scoring prefix. Ties retain legacy ordering.
"""
    n, m = len(a), len(b)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    matches = [[0] * (m + 1) for _ in range(n + 1)]
    lengths = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0], lengths[i][0] = -3 * i, i
    for j in range(1, m + 1):
        score[0][j], lengths[0][j] = -3 * j, j
    best, passing = None, None
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            same = a[i - 1] == b[j - 1] and a[i - 1] in "ACGT"
            choices = (score[i - 1][j - 1] + (2 if same else -2),
                       score[i - 1][j] - 3, score[i][j - 1] - 3)
            step = max(range(3), key=choices.__getitem__)
            x, y = i - (step != 2), j - (step != 1)
            score[i][j] = choices[step]
            matches[i][j] = matches[x][y] + int(step == 0 and same)
            lengths[i][j] = lengths[x][y] + 1
            record = (score[i][j], matches[i][j], -lengths[i][j], i, j)
            if best is None or record > best:
                best = record
            if matches[i][j] >= MIN_MATCHES and matches[i][j] / lengths[i][j] >= MIN_IDENTITY:
                if passing is None or record > passing:
                    passing = record
    if best is None:
        return dict(pass_sequence=False, score=0, matches=0, length=0, identity=0.0)
    value, matched, negative_length, _, _ = passing or best
    length = -negative_length
    return dict(pass_sequence=passing is not None, score=value, matches=matched,
                length=length, identity=matched / length)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def position(text):
    chrom, coordinate = text.rsplit(":", 1)
    return chrom, int(coordinate)


def oriented_context(reference, endpoint, direction, first):
    """A-retained|A-unused or B-unused|B-retained, 50 bases per side.

The breakpoint is the last retained base for '-' and the first retained
base for '+'. Fetch is zero-based and half-open. Pad contig edges with N,
which can never contribute a matching base to R2.
"""
    chrom, bp = position(endpoint)
    if direction not in ("+", "-") or chrom not in reference.references:
        return None
    length = reference.get_reference_length(chrom)
    if not 1 <= bp <= length:
        return None
    cut = bp if direction == "-" else bp - 1
    left = reference.fetch(chrom, max(0, cut - WINDOW), cut).upper().rjust(WINDOW, "N")
    right = reference.fetch(chrom, cut, min(length, cut + WINDOW)).upper().ljust(WINDOW, "N")
    context = left + right
    return reverse_complement(context) if (first and direction == "+") or (not first and direction == "-") else context


def sequence_metrics(a, b):
    # Reverse both left halves so the first base is next to the breakpoint.
    # Do not complement them: both are already in the same fusion orientation.
    return {"left": anchored_prefix(a[:WINDOW][::-1], b[:WINDOW][::-1]),
            "right": anchored_prefix(a[WINDOW:], b[WINDOW:])}


@dataclass(frozen=True)
class Gene:
    aliases: frozenset
    chrom: str
    start: int
    end: int
    strand: str
    biotype: str

    @classmethod
    def from_fields(cls, fields, biotype):
        return cls(frozenset(gene_aliases(fields[8])), fields[0], int(fields[3]),
                   int(fields[4]), fields[6], biotype)


class GeneIndex:
    """Gene features only: no union across loci and no inferred gene spans."""
    def __init__(self, genes):
        self.partners = defaultdict(list)
        coding = defaultdict(list)
        for gene in genes:
            for name in gene.aliases:
                self.partners[name].append(gene)
            if gene.biotype == "protein_coding" and gene.strand in ("+", "-"):
                coding[(gene.chrom, gene.strand)].append(gene)
        self.coding = {}
        for key, intervals in coding.items():
            intervals.sort(key=lambda gene: (gene.start, gene.end, sorted(gene.aliases)))
            maximum, maxima = 0, []
            for gene in intervals:
                maximum = max(maximum, gene.end)
                maxima.append(maximum)
            self.coding[key] = (intervals, [gene.start for gene in intervals], maxima)

    def intervening(self, left, right):
        start, end = left.end + 1, right.start - 1
        intervals, starts, maxima = self.coding.get((left.chrom, left.strand), ([], [], []))
        excluded = left.aliases | right.aliases
        hits = intervals[bisect_left(maxima, start):bisect_right(starts, end)] if start <= end else []
        return [gene for gene in hits if gene.end >= start and not gene.aliases & excluded]

    def readthrough(self, row):
        if any(row.get("direction" + str(side)) not in ("+", "-") for side in (1, 2)):
            return "UNKNOWN_DIRECTION", "NA"
        endpoints = [position(row["Position" + str(side)]) for side in (1, 2)]
        (c1, p1), (c2, p2) = endpoints
        low_side, high_side = ((row["direction1"], row["direction2"]) if p1 <= p2
                               else (row["direction2"], row["direction1"]))
        if c1 != c2 or (low_side, high_side) != ("-", "+"):
            return "PASS", ""
        names = row["FusionName"].split("--", 1)
        partners = [[gene for gene in self.partners[name]
                     if gene.chrom == chrom and gene.start <= bp <= gene.end]
                    for name, (chrom, bp) in zip(names, endpoints)]
        if not all(partners):
            return "UNKNOWN_ANNOTATION", "NA"
        # Missing or ambiguous table strands are resolved only when the
        # breakpoint's annotated gene features give a unique known strand.
        for side, records in enumerate(partners, 1):
            strand = row.get("strand" + str(side), "?")
            if strand not in ("+", "-"):
                strands = {gene.strand for gene in records}
                strand = next(iter(strands)) if len(strands) == 1 else "?"
            if strand not in ("+", "-"):
                return "UNKNOWN_ANNOTATION", "NA"
            partners[side - 1] = [gene for gene in records if gene.strand == strand]
        if not all(partners):
            return "UNKNOWN_ANNOTATION", "NA"
        blockers = set()
        for a in partners[0]:
            for b in partners[1]:
                left, right = sorted((a, b), key=lambda gene: gene.start)
                if left.strand != right.strand or left.end >= right.start:
                    continue
                if "lncRNA" not in (left.biotype, right.biotype):
                    continue
                intervening = self.intervening(left, right)
                if not intervening:
                    return "FILTERED", ""
                blockers.update(name for gene in intervening for name in gene.aliases)
        return "PASS", ";".join(sorted(blockers))


def empty_status(value="NOT_APPLIED"):
    return {field: value if field in ("R2Status", "R3Status") else "NA"
            for field in STATUS_COLUMNS if field not in ("Baseline1Selected", "Baseline3Rules")}


class Screens:
    def __init__(self, fasta_path, genes):
        self.fasta_path = Path(fasta_path) if fasta_path is not None else None
        self.reference = None
        self.genes = GeneIndex(genes)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        if self.reference is not None:
            self.reference.close()

    def r2(self, row):
        details = empty_status()
        details["R2SplitSpanningRatio"] = int(row["JunctionReadCount"]) / (int(row["SpanningFragCount"]) + 1)
        if any(row.get("direction" + str(side)) not in ("+", "-") for side in (1, 2)):
            details["R2Status"] = "UNKNOWN_DIRECTION"
            return details
        if details["R2SplitSpanningRatio"] < MIN_SPLIT_SPANNING_RATIO:
            details["R2Status"] = "RATIO_BELOW_CUTOFF"
            return details
        if self.reference is None:
            import pysam
            if self.fasta_path is None or not Path(str(self.fasta_path) + ".fai").is_file():
                raise FileNotFoundError("R2 requires GenomeDir/ref.fa with an existing .fai; prepare it with scFusion.py Index")
            self.reference = pysam.FastaFile(str(self.fasta_path))
        a = oriented_context(self.reference, row["Position1"], row["direction1"], True)
        b = oriented_context(self.reference, row["Position2"], row["direction2"], False)
        if a is None or b is None:
            details["R2Status"] = "UNKNOWN_REFERENCE"
            return details
        metrics = sequence_metrics(a, b)
        for side in ("left", "right"):
            for name in ("matches", "length", "identity"):
                details["R2" + side.title() + name.title()] = metrics[side][name]
        details["R2Status"] = "FILTERED" if any(value["pass_sequence"] for value in metrics.values()) else "PASS"
        return details

    def evaluate(self, row):
        details = self.r2(row)
        details["R3Status"], details["R3InterveningProteinCodingGenes"] = self.genes.readthrough(row)
        rules = [rule for rule in ("R2", "R3") if details[rule + "Status"] == "FILTERED"]
        details["Baseline3Rules"] = ";".join(rules) if rules else "none"
        reasons = [reason for rule, reason in (("R2", R2_REASON), ("R3", R3_REASON)) if rule in rules]
        return details, reasons
