#!/usr/bin/env python3
"""Dependency-free SAM/CIGAR helpers for HotSpotCluster."""

from __future__ import annotations

import re
from dataclasses import dataclass
from itertools import combinations
from typing import Iterable, List, Optional, Sequence, Tuple


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
REFERENCE_OPS = frozenset("MDN=X")
QUERY_OPS = frozenset("MIS=X")


def parse_cigar(cigar: str) -> List[Tuple[int, str]]:
    if cigar == "*":
        return []
    result = [(int(length), op) for length, op in CIGAR_RE.findall(cigar)]
    if not result or "".join(f"{length}{op}" for length, op in result) != cigar:
        raise ValueError(f"invalid CIGAR: {cigar}")
    return result


def reference_interval(pos: int, cigar: str) -> Tuple[int, int]:
    """Return a 1-based closed interval; M/D/N/=/X consume reference."""
    length = sum(n for n, op in parse_cigar(cigar) if op in REFERENCE_OPS)
    if pos < 1 or length < 1:
        raise ValueError(f"empty alignment: POS={pos}, CIGAR={cigar}")
    return pos, pos + length - 1


def terminal_clip_side(cigar: str) -> Optional[str]:
    """Return genomic CIGAR side: '+' is right and '-' is left."""
    ops = parse_cigar(cigar)
    if not ops:
        return None
    left = ops[0][0] if ops[0][1] in "SH" else 0
    right = ops[-1][0] if ops[-1][1] in "SH" else 0
    if left == right:
        return None
    return "+" if right > left else "-"


def query_bounds(cigar: str) -> Tuple[int, int, int]:
    ops = parse_cigar(cigar)
    leading = 0
    for n, op in ops:
        if op not in "SH":
            break
        leading += n
    trailing = 0
    for n, op in reversed(ops):
        if op not in "SH":
            break
        trailing += n
    total = sum(n for n, op in ops if op in QUERY_OPS or op == "H")
    return leading, total - trailing, total


def mate_number(flag: int) -> int:
    if flag & 0x40:
        return 1
    if flag & 0x80:
        return 2
    return 0


def strand(flag: int) -> str:
    return "-" if flag & 0x10 else "+"


@dataclass(frozen=True)
class SamRecord:
    raw: str
    fields: Tuple[str, ...]

    @property
    def qname(self):
        return self.fields[0]

    @property
    def flag(self):
        return int(self.fields[1])

    @property
    def chrom(self):
        return self.fields[2]

    @property
    def pos(self):
        return int(self.fields[3])

    @property
    def cigar(self):
        return self.fields[5]

    @property
    def sequence(self):
        return self.fields[9]


def aligned_query_length(record: SamRecord) -> int:
    """Return query bases aligned by M, = or X; clips and insertions are excluded."""
    return sum(n for n, op in parse_cigar(record.cigar) if op in "M=X")


def aligned_edge_sequence(record: SamRecord, edge: str, length: int) -> str:
    """Return aligned query bases nearest the genomic left or right interval edge."""
    if record.sequence == "*" or edge not in {"LEFT", "RIGHT"}:
        return ""
    query_position = 0
    aligned_bases = []
    for count, operation in parse_cigar(record.cigar):
        if operation in "M=X":
            aligned_bases.append(record.sequence[query_position:query_position + count])
            query_position += count
        elif operation in "SI":
            query_position += count
        # D/N/H/P consume no query sequence.
    sequence = "".join(aligned_bases).upper()
    if len(sequence) < length:
        return ""
    return sequence[:length] if edge == "LEFT" else sequence[-length:]


def grouped_sam(handle: Iterable[str]):
    name = None
    records: List[SamRecord] = []
    for line in handle:
        if not line or line.startswith("@"):
            continue
        fields = tuple(line.rstrip("\r\n").split("\t"))
        if len(fields) < 11:
            continue
        record = SamRecord(line, fields)
        if name is not None and record.qname != name:
            yield name, records
            records = []
        name = record.qname
        records.append(record)
    if name is not None:
        yield name, records


def choose_same_mate_split(records: Sequence[SamRecord]):
    """Choose the complementary clipped alignments from a three-record group."""
    if len(records) != 3:
        return None
    candidates = []
    for first, second in combinations(records, 2):
        mate = mate_number(first.flag)
        if mate == 0 or mate != mate_number(second.flag):
            continue
        try:
            side1 = terminal_clip_side(first.cigar)
            side2 = terminal_clip_side(second.cigar)
            if side1 is None or side2 is None:
                continue
            a1, b1, length1 = query_bounds(first.cigar)
            a2, b2, length2 = query_bounds(second.cigar)
            if length1 != length2:
                continue
            overlap = max(0, min(b1, b2) - max(a1, a2))
            gap = max(0, max(a1, a2) - min(b1, b2))
            same_locus = first.chrom == second.chrom and abs(first.pos - second.pos) < 1000
            candidates.append((overlap + gap + (1000000 if same_locus else 0), first, second))
        except ValueError:
            continue
    if not candidates:
        return None
    candidates.sort(key=lambda item: item[0])
    return candidates[0][1], candidates[0][2]


def alignment_values(record: SamRecord):
    start, end = reference_interval(record.pos, record.cigar)
    genomic_side = terminal_clip_side(record.cigar)
    if genomic_side is None:
        raise ValueError("alignment has no unambiguous terminal clip")
    read_strand = strand(record.flag)
    reported_side = genomic_side
    if read_strand == "-":
        reported_side = "-" if genomic_side == "+" else "+"
    junction = end if genomic_side == "+" else start
    return start, end, read_strand, reported_side, junction

