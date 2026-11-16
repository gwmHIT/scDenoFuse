"""Global hotspot filtering ported from scFusion v2 without persistent read tracing."""

from __future__ import annotations

import bisect
import csv
import os
from concurrent.futures import ThreadPoolExecutor
from tempfile import TemporaryDirectory
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

from .sam_interval import (
    aligned_edge_sequence, aligned_query_length, alignment_values,
    choose_same_mate_split, grouped_sam,
)

IntervalKey = Tuple[str, int, int]
EdgeSignature = Tuple[str, str, str, str]

@dataclass
class EdgeFamily:
    """一个固定区间边界下，由不同长度的嵌套短匹配组成的候选家族。"""

    provisional_id: int
    chrom: str
    strand: str
    clip: str
    fixed_edge: str
    fixed_coordinate: int
    member_coordinates: Set[int]
    observation_count: int
    local_observation_count: int
    boundary_fraction: float
    dominant_sequence: str
    sequence_count: int
    sequence_fraction: float
    reads: Set[Tuple[str, str]] = field(default_factory=set)
    cells: Set[str] = field(default_factory=set)
    partners: Set[Tuple[str, int]] = field(default_factory=set)
    anchor_lengths: Counter = field(default_factory=Counter)
    decision: str = "PENDING"
    reject_reason: str = ""

@dataclass
class Mode:
    provisional_id: int
    chrom: str
    modal_start: int
    modal_end: int
    member_keys: Set[IntervalKey]
    interval_start: int
    interval_end: int
    observation_count: int
    dominant_strand: str = ""
    dominant_clip: str = ""
    modal_breakpoint: int = -1
    modal_observation_count: int = 0
    core_observation_count: int = 0
    local_observation_count: int = 0
    modal_fraction: float = 0.0
    core_fraction: float = 0.0
    decision: str = "PENDING"
    reject_reason: str = ""
    left_gap: str = "NA"
    right_gap: str = "NA"

@dataclass
class ModeSummary:
    """由区间成员记录派生的聚类级汇总，不保存成员来源本身。"""

    reads: Set[Tuple[str, str]] = field(default_factory=set)
    cells: Set[str] = field(default_factory=set)
    core_reads: Set[Tuple[str, str]] = field(default_factory=set)
    core_cells: Set[str] = field(default_factory=set)
    expanded_reads: Set[Tuple[str, str]] = field(default_factory=set)
    partners: Set[Tuple[str, int]] = field(default_factory=set)
    strands: Set[str] = field(default_factory=set)
    clips: Set[str] = field(default_factory=set)
    breakpoint_min: int = 2**63 - 1
    breakpoint_max: int = -1

def resolve_input(outdir: Path, pattern: str, cell: int) -> Path:
    rendered = Path(pattern.format(cell=cell))
    return rendered if rendered.is_absolute() else outdir / rendered

def extract_evidence(
    outdir: Path,
    pattern: str,
    cells: Sequence[int],
    evidence_path: Path,
    edge_sequence_length: int,
):
    interval_counts: Counter = Counter()
    profile_counts: Counter = Counter()
    breakpoint_counts: Counter = Counter()
    source_paths = {}
    group_stats = Counter()
    with evidence_path.open("w", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "ReadName", "Cell", "Chrom1", "Start1", "End1", "Chrom2", "Start2", "End2",
                "Strand1", "Strand2", "SoftclipSide1", "SoftclipSide2", "Breakpoint1", "Breakpoint2",
                "Cigar1", "Cigar2", "AnchorLength1", "AnchorLength2",
                "LeftEdgeSequence1", "RightEdgeSequence1",
                "LeftEdgeSequence2", "RightEdgeSequence2", "Evidence",
            ]
        )
        for cell in cells:
            source = resolve_input(outdir, pattern, cell)
            if not source.exists():
                raise FileNotFoundError(f"cell {cell} input SAM not found: {source}")
            source_paths[str(cell)] = source
            with source.open(errors="replace") as handle:
                for readname, records in grouped_sam(handle):
                    group_stats[f"groups_{len(records)}"] += 1
                    pair = choose_same_mate_split(records)
                    if pair is None:
                        continue
                    try:
                        first, second = pair
                        start1, end1, strand1, clip1, point1 = alignment_values(first)
                        start2, end2, strand2, clip2, point2 = alignment_values(second)
                    except (ValueError, IndexError):
                        group_stats["invalid_split"] += 1
                        continue
                    interval_counts[(first.chrom, start1, end1)] += 1
                    interval_counts[(second.chrom, start2, end2)] += 1
                    profile_counts[(first.chrom, start1, end1, strand1, clip1, point1)] += 1
                    profile_counts[(second.chrom, start2, end2, strand2, clip2, point2)] += 1
                    breakpoint_counts[(first.chrom, strand1, clip1, point1)] += 1
                    breakpoint_counts[(second.chrom, strand2, clip2, point2)] += 1
                    evidence = ";".join(
                        [
                            readname,
                            str(cell),
                            f"{first.chrom}:{start1}-{end1}",
                            f"{second.chrom}:{start2}-{end2}",
                            strand1,
                            strand2,
                            clip1,
                            clip2,
                        ]
                    )
                    writer.writerow(
                        [
                            readname, cell, first.chrom, start1, end1, second.chrom, start2, end2,
                            strand1, strand2, clip1, clip2, point1, point2,
                            first.cigar, second.cigar,
                            aligned_query_length(first), aligned_query_length(second),
                            aligned_edge_sequence(first, "LEFT", edge_sequence_length),
                            aligned_edge_sequence(first, "RIGHT", edge_sequence_length),
                            aligned_edge_sequence(second, "LEFT", edge_sequence_length),
                            aligned_edge_sequence(second, "RIGHT", edge_sequence_length),
                            evidence,
                        ]
                    )
                    group_stats["split_evidence"] += 1
    return interval_counts, profile_counts, breakpoint_counts, source_paths, group_stats

def build_modes(interval_counts: Counter, tight_window: int, min_reads: int):
    by_chrom: Dict[str, List[IntervalKey]] = defaultdict(list)
    for key in interval_counts:
        by_chrom[key[0]].append(key)
    modes: List[Mode] = []
    next_id = 1
    width = tight_window
    for chrom, keys in by_chrom.items():
        grid: Dict[Tuple[int, int], List[IntervalKey]] = defaultdict(list)
        for key in keys:
            grid[(key[1] // width, key[2] // width)].append(key)
        unassigned = set(keys)
        for center in sorted(keys, key=lambda key: (-interval_counts[key], key[1], key[2])):
            if center not in unassigned:
                continue
            cell = (center[1] // width, center[2] // width)
            nearby_keys = set()
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for key in grid.get((cell[0] + dx, cell[1] + dy), []):
                        if (
                            key in unassigned
                            and abs(key[1] - center[1]) < tight_window
                            and abs(key[2] - center[2]) < tight_window
                        ):
                            nearby_keys.add(key)

            member_keys = {center}
            for key in sorted(
                nearby_keys - {center},
                key=lambda item: (-interval_counts[item], item[1], item[2]),
            ):
                candidate_keys = member_keys | {key}
                starts = [item[1] for item in candidate_keys]
                ends = [item[2] for item in candidate_keys]
                if max(starts) - min(starts) < tight_window and max(ends) - min(ends) < tight_window:
                    member_keys.add(key)
            unassigned.difference_update(member_keys)
            count = sum(interval_counts[key] for key in member_keys)
            if count < min_reads:
                continue
            modes.append(
                Mode(
                    next_id,
                    chrom,
                    center[1],
                    center[2],
                    member_keys,
                    min(key[1] for key in member_keys),
                    max(key[2] for key in member_keys),
                    count,
                )
            )
            next_id += 1

    return modes

def score_modes_by_core_fraction(
    modes: List[Mode],
    profile_counts: Counter,
    breakpoint_counts: Counter,
    local_window: int,
    min_core_fraction: float,
):
    """计算严格核心在同染色体、同方向局部断点证据中的占比。"""
    key_to_mode = {}
    mode_by_id = {mode.provisional_id: mode for mode in modes}
    for mode in modes:
        for key in mode.member_keys:
            key_to_mode[key] = mode.provisional_id

    orientation_counts: Dict[int, Counter] = defaultdict(Counter)
    modal_profiles: Dict[int, Counter] = defaultdict(Counter)
    for profile, count in profile_counts.items():
        chrom, start, end, strand_value, clip, point = profile
        mode_id = key_to_mode.get((chrom, start, end))
        if mode_id is None:
            continue
        orientation_counts[mode_id][(strand_value, clip)] += count
        mode = mode_by_id[mode_id]
        if start == mode.modal_start and end == mode.modal_end:
            modal_profiles[mode_id][(strand_value, clip, point)] += count

    local_profiles: Dict[Tuple[str, str, str], Counter] = defaultdict(Counter)
    for (chrom, strand_value, clip, point), count in breakpoint_counts.items():
        local_profiles[(chrom, strand_value, clip)][point] += count
    local_coords = {}
    local_prefix = {}
    for signature, point_counts in local_profiles.items():
        coords = sorted(point_counts)
        prefix = [0]
        for point in coords:
            prefix.append(prefix[-1] + point_counts[point])
        local_coords[signature] = coords
        local_prefix[signature] = prefix

    passing = []
    for mode in modes:
        profiles = orientation_counts.get(mode.provisional_id)
        if not profiles:
            mode.decision = "REJECTED"
            mode.reject_reason = "NO_DIRECTION_PROFILE"
            continue
        dominant = max(profiles, key=lambda key: (profiles[key], key))
        mode.dominant_strand, mode.dominant_clip = dominant
        mode.core_observation_count = profiles[dominant]

        modal_counter = Counter()
        for (strand_value, clip, point), count in modal_profiles[mode.provisional_id].items():
            if (strand_value, clip) == dominant:
                modal_counter[point] += count
        if not modal_counter:
            mode.decision = "REJECTED"
            mode.reject_reason = "NO_MODAL_DIRECTION_PROFILE"
            continue
        mode.modal_breakpoint = max(modal_counter, key=lambda point: (modal_counter[point], -point))
        mode.modal_observation_count = sum(modal_counter.values())

        signature = (mode.chrom, mode.dominant_strand, mode.dominant_clip)
        coords = local_coords.get(signature, [])
        prefix = local_prefix.get(signature, [0])
        left = bisect.bisect_left(coords, mode.modal_breakpoint - local_window)
        right = bisect.bisect_right(coords, mode.modal_breakpoint + local_window)
        mode.local_observation_count = prefix[right] - prefix[left]
        if mode.local_observation_count:
            mode.modal_fraction = mode.modal_observation_count / mode.local_observation_count
            mode.core_fraction = mode.core_observation_count / mode.local_observation_count
        if mode.core_fraction < min_core_fraction:
            mode.decision = "REJECTED"
            mode.reject_reason = "CORE_FRACTION"
            continue
        passing.append(mode)
    return passing

def collect_members(
    evidence_path: Path,
    modes: List[Mode],
    provisional_path: Path,
    partner_bin: int,
    filter_breakpoint_window: int,
):
    key_to_mode = {}
    signature_modes: Dict[Tuple[str, str, str], List[Tuple[int, int]]] = defaultdict(list)
    for mode in modes:
        for key in mode.member_keys:
            key_to_mode[(key, mode.dominant_strand, mode.dominant_clip)] = mode.provisional_id
        signature_modes[(mode.chrom, mode.dominant_strand, mode.dominant_clip)].append(
            (mode.modal_breakpoint, mode.provisional_id)
        )
    signature_coords = {}
    mode_by_id = {mode.provisional_id: mode for mode in modes}
    for signature in signature_modes:
        signature_modes[signature].sort()
        signature_coords[signature] = [item[0] for item in signature_modes[signature]]

    def find_expanded_mode(chrom, strand_value, clip, point):
        signature = (chrom, strand_value, clip)
        candidates = signature_modes.get(signature, [])
        coords = signature_coords.get(signature, [])
        left = bisect.bisect_left(coords, point - filter_breakpoint_window)
        right = bisect.bisect_right(coords, point + filter_breakpoint_window)
        if left == right:
            return None
        return min(
            candidates[left:right],
            key=lambda item: (
                abs(item[0] - point),
                -mode_by_id[item[1]].core_observation_count,
                item[1],
            ),
        )[1]

    with evidence_path.open() as inp, provisional_path.open("w", newline="") as output:
        reader = csv.DictReader(inp, delimiter="\t")
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "ProvisionalModeID", "ReadName", "Cell", "Side", "Chrom", "Start", "End",
                "Strand", "SoftclipSide", "Breakpoint", "PartnerChrom", "PartnerBreakpoint",
                "PartnerBin", "Membership",
            ]
        )
        for row in reader:
            sides = [
                (
                    1,
                    (row["Chrom1"], int(row["Start1"]), int(row["End1"])),
                    row["Strand1"], row["SoftclipSide1"], int(row["Breakpoint1"]),
                    row["Chrom2"], int(row["Breakpoint2"]),
                ),
                (
                    2,
                    (row["Chrom2"], int(row["Start2"]), int(row["End2"])),
                    row["Strand2"], row["SoftclipSide2"], int(row["Breakpoint2"]),
                    row["Chrom1"], int(row["Breakpoint1"]),
                ),
            ]
            for side, key, strand_value, clip, point, partner_chrom, partner_point in sides:
                mode_id = key_to_mode.get((key, strand_value, clip))
                membership = "CORE"
                if mode_id is None:
                    mode_id = find_expanded_mode(key[0], strand_value, clip, point)
                    membership = "EXPANDED"
                    if mode_id is None:
                        continue
                writer.writerow(
                    [
                        mode_id, row["ReadName"], row["Cell"], side, key[0], key[1], key[2],
                        strand_value, clip, point, partner_chrom, partner_point,
                        partner_point // partner_bin, membership,
                    ]
                )

def evidence_sides(row):
    """把一行双侧 evidence 展开成两个结构相同、便于聚类的字典。"""
    sides = []
    for side, partner in ((1, 2), (2, 1)):
        sides.append(
            {
                "side": side,
                "chrom": row[f"Chrom{side}"],
                "start": int(row[f"Start{side}"]),
                "end": int(row[f"End{side}"]),
                "strand": row[f"Strand{side}"],
                "clip": row[f"SoftclipSide{side}"],
                "breakpoint": int(row[f"Breakpoint{side}"]),
                "anchor_length": int(row[f"AnchorLength{side}"]),
                "left_sequence": row[f"LeftEdgeSequence{side}"],
                "right_sequence": row[f"RightEdgeSequence{side}"],
                "partner_chrom": row[f"Chrom{partner}"],
                "partner_breakpoint": int(row[f"Breakpoint{partner}"]),
            }
        )
    return sides

def outer_edge(side):
    """返回远离soft-clip的区间边界；真实断点所在的clip边界不参与补充检测。"""
    genomic_clip = side["clip"]
    if side["strand"] == "-":
        genomic_clip = "+" if genomic_clip == "-" else "-"
    if genomic_clip == "+":
        return "LEFT", side["start"], side["left_sequence"]
    return "RIGHT", side["end"], side["right_sequence"]

def build_edge_families(
    evidence_path: Path,
    boundary_window: int,
    local_window: int,
    min_anchor: int,
    max_anchor: int,
    min_reads: int,
    min_boundary_fraction: float,
    min_sequence_fraction: float,
):
    """合并共享左边界或右边界、但匹配长度不同的嵌套短片段。"""
    edge_counts: Dict[EdgeSignature, Counter] = defaultdict(Counter)
    sequence_counts: Dict[Tuple[EdgeSignature, int], Counter] = defaultdict(Counter)

    with evidence_path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for side in evidence_sides(row):
                if not min_anchor <= side["anchor_length"] <= max_anchor:
                    continue
                edge, coordinate, sequence = outer_edge(side)
                if not sequence:
                    continue
                signature = (side["chrom"], side["strand"], side["clip"], edge)
                edge_counts[signature][coordinate] += 1
                sequence_counts[(signature, coordinate)][sequence] += 1

    families = []
    next_id = 1
    for signature, coordinate_counts in edge_counts.items():
        coordinates = sorted(coordinate_counts)
        prefix_counts = [0]
        for coordinate in coordinates:
            prefix_counts.append(prefix_counts[-1] + coordinate_counts[coordinate])
        unassigned = set(coordinate_counts)
        for center in sorted(coordinate_counts, key=lambda x: (-coordinate_counts[x], x)):
            if center not in unassigned:
                continue
            left = bisect.bisect_left(coordinates, center - boundary_window)
            right = bisect.bisect_right(coordinates, center + boundary_window)
            member_coordinates = {
                coordinate
                for coordinate in coordinates[left:right]
                if coordinate in unassigned
            }
            observation_count = sum(coordinate_counts[x] for x in member_coordinates)
            sequences = Counter()
            if observation_count >= min_reads:
                for coordinate in member_coordinates:
                    sequences.update(sequence_counts[(signature, coordinate)])
                leading = sequences.most_common(2)
                if (
                    len(member_coordinates) > 1
                    and len(leading) > 1
                    and leading[0][1] == leading[1][1]
                ):
                    # The original Counter tie order follows a comprehension
                    # over the mutable unassigned set, not coordinate order.
                    # Recreate it only when it can change the selected sequence.
                    member_coordinates = {
                        coordinate
                        for coordinate in unassigned
                        if abs(coordinate - center) <= boundary_window
                    }
                    sequences.clear()
                    for coordinate in member_coordinates:
                        sequences.update(sequence_counts[(signature, coordinate)])
            unassigned.difference_update(member_coordinates)
            if observation_count < min_reads:
                continue
            left = bisect.bisect_left(coordinates, center - local_window)
            right = bisect.bisect_right(coordinates, center + local_window)
            local_count = prefix_counts[right] - prefix_counts[left]
            dominant_sequence, sequence_count = sequences.most_common(1)[0]
            boundary_fraction = observation_count / local_count
            sequence_fraction = sequence_count / observation_count
            family = EdgeFamily(
                next_id, *signature, center, member_coordinates,
                observation_count, local_count, boundary_fraction,
                dominant_sequence, sequence_count, sequence_fraction,
            )
            reasons = []
            if boundary_fraction < min_boundary_fraction:
                reasons.append("BOUNDARY_FRACTION")
            if sequence_fraction < min_sequence_fraction:
                reasons.append("SEQUENCE_FRACTION")
            if reasons:
                family.decision = "REJECTED"
                family.reject_reason = ",".join(reasons)
            families.append(family)
            next_id += 1
    return families



def summarize_modes(
    modes, provisional_path, output_dir, min_reads=5, min_cells=2,
    min_partners=None, min_interval_length=10, max_interval_length=20, diagnostics=True,
):
    """Apply the original strict-mode membership thresholds."""
    summaries = {mode.provisional_id: ModeSummary() for mode in modes}
    with provisional_path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mode_id = int(row["ProvisionalModeID"])
            summary = summaries[mode_id]
            key = (row["Cell"], row["ReadName"])
            summary.reads.add(key)
            summary.cells.add(row["Cell"])
            if row["Membership"] == "CORE":
                summary.core_reads.add(key)
                summary.core_cells.add(row["Cell"])
            else:
                summary.expanded_reads.add(key)
            summary.partners.add((row["PartnerChrom"], int(row["PartnerBin"])))
            summary.strands.add(row["Strand"])
            summary.clips.add(row["SoftclipSide"])
            point = int(row["Breakpoint"])
            summary.breakpoint_min = min(summary.breakpoint_min, point)
            summary.breakpoint_max = max(summary.breakpoint_max, point)

    required = min_partners if min_partners is not None else 0
    valid = []
    for mode in modes:
        summary = summaries[mode.provisional_id]
        modal_length = mode.modal_end - mode.modal_start + 1
        reasons = []
        if len(summary.core_reads) < min_reads:
            reasons.append("MIN_CORE_READS")
        if len(summary.core_cells) < min_cells:
            reasons.append("MIN_CORE_CELLS")
        if len(summary.partners) < required:
            reasons.append("MIN_PARTNER_LOCI")
        if not min_interval_length <= modal_length <= max_interval_length:
            reasons.append("MODAL_INTERVAL_LENGTH")
        mode.decision = "REJECTED" if reasons else "ACCEPTED"
        mode.reject_reason = ",".join(reasons)
        if not reasons:
            valid.append(mode)
    valid.sort(key=lambda mode: (mode.chrom, mode.modal_start, mode.modal_end))
    filtered_names = defaultdict(set)
    destination = output_dir / "hotspot_clusters.tsv" if diagnostics else Path(os.devnull)
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "ClusterID", "Chrom", "ModalStart", "ModalEnd", "ModalIntervalLength",
            "IntervalStart", "IntervalEnd", "BreakpointMin", "BreakpointMax",
            "Strand", "SoftclipSide", "ReadCount", "CellCount", "CoreReadCount",
            "CoreCellCount", "ExpandedOnlyReadCount", "ModalBreakpoint",
            "ModalObservationCount", "CoreObservationCount", "LocalObservationCount",
            "ModalFraction", "CoreFraction", "PartnerLocusCount",
            "RequiredPartnerLocusCount", "LeftGap", "RightGap",
        ])
        for index, mode in enumerate(valid, 1):
            summary = summaries[mode.provisional_id]
            writer.writerow([
                index, mode.chrom, mode.modal_start, mode.modal_end,
                mode.modal_end - mode.modal_start + 1, mode.interval_start,
                mode.interval_end, summary.breakpoint_min, summary.breakpoint_max,
                ",".join(sorted(summary.strands)), ",".join(sorted(summary.clips)),
                len(summary.reads), len(summary.cells), len(summary.core_reads),
                len(summary.core_cells), len(summary.reads - summary.core_reads),
                mode.modal_breakpoint, mode.modal_observation_count,
                mode.core_observation_count, mode.local_observation_count,
                f"{mode.modal_fraction:.6f}", f"{mode.core_fraction:.6f}",
                len(summary.partners), required, mode.left_gap, mode.right_gap,
            ])
            for cell, readname in summary.reads:
                filtered_names[cell].add(readname)
    return valid, filtered_names


def summarize_edge_families(
    evidence_path, families, output_dir, filtered_names, min_reads=5,
    min_cells=2, partner_bin=150, min_anchor=10, max_anchor=30, diagnostics=True,
):
    """Collect the exact original edge-family members without writing read audits."""
    passing = [family for family in families if family.decision == "PENDING"]
    coordinate_to_family = {}
    family_by_id = {family.provisional_id: family for family in passing}
    for family in passing:
        signature = (family.chrom, family.strand, family.clip, family.fixed_edge)
        for coordinate in family.member_coordinates:
            coordinate_to_family[(signature, coordinate)] = family.provisional_id
    with evidence_path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for side in evidence_sides(row):
                if not min_anchor <= side["anchor_length"] <= max_anchor:
                    continue
                edge, coordinate, sequence = outer_edge(side)
                signature = (side["chrom"], side["strand"], side["clip"], edge)
                family = family_by_id.get(coordinate_to_family.get((signature, coordinate)))
                if family is None or sequence != family.dominant_sequence:
                    continue
                family.reads.add((row["Cell"], row["ReadName"]))
                family.cells.add(row["Cell"])
                family.partners.add((side["partner_chrom"], side["partner_breakpoint"] // partner_bin))
                family.anchor_lengths[side["anchor_length"]] += 1

    valid = []
    for family in passing:
        reasons = []
        if len(family.reads) < min_reads:
            reasons.append("MIN_READS")
        if len(family.cells) < min_cells:
            reasons.append("MIN_CELLS")
        family.decision = "REJECTED" if reasons else "ACCEPTED"
        family.reject_reason = ",".join(reasons)
        if not reasons:
            valid.append(family)
    destination = output_dir / "adapter_hotspot_families.tsv" if diagnostics else Path(os.devnull)
    with destination.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "FamilyID", "Chrom", "FixedEdge", "FixedCoordinate", "BoundaryMin",
            "BoundaryMax", "Strand", "SoftclipSide", "ReadCount", "CellCount",
            "ObservationCount", "LocalObservationCount", "BoundaryFraction",
            "AnchorLengthMin", "AnchorLengthMax", "DominantAnchorLength",
            "DominantEdgeSequence", "SequenceCount", "SequenceFraction",
            "PartnerLocusCount", "Decision", "RejectReason",
        ])
        for family in families:
            lengths = family.anchor_lengths
            writer.writerow([
                family.provisional_id, family.chrom, family.fixed_edge,
                family.fixed_coordinate, min(family.member_coordinates),
                max(family.member_coordinates), family.strand, family.clip,
                len(family.reads), len(family.cells), family.observation_count,
                family.local_observation_count, f"{family.boundary_fraction:.6f}",
                min(lengths) if lengths else "NA", max(lengths) if lengths else "NA",
                lengths.most_common(1)[0][0] if lengths else "NA",
                family.dominant_sequence, family.sequence_count,
                f"{family.sequence_fraction:.6f}", len(family.partners),
                family.decision, family.reject_reason,
            ])
    old_keys = {(cell, name) for cell, names in filtered_names.items() for name in names}
    edge_keys = set().union(*(family.reads for family in valid)) if valid else set()
    for cell, readname in edge_keys:
        filtered_names[cell].add(readname)
    return valid, edge_keys, edge_keys - old_keys, edge_keys & old_keys


def filter_sam(source, names):
    """Replace the new working SAM, never the STAR file or a reuse source."""
    if not names:
        return 0
    temporary = source.with_name(f".{source.name}.hotspot.tmp")
    removed = 0
    try:
        with source.open(errors="replace") as inp, temporary.open("w") as output:
            for line in inp:
                if not line.startswith("@") and line.split("\t", 1)[0] in names:
                    removed += 1
                else:
                    output.write(line)
        os.replace(temporary, source)
    finally:
        temporary.unlink(missing_ok=True)
    return removed


def write_summary(output_dir, values, diagnostics=True):
    if not diagnostics:
        return
    with (output_dir / "run_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Metric", "Value"])
        writer.writerows(values.items())


def run(args):
    """Global hotspot decision followed by bounded concurrent working-SAM filtering.

    Input: args.outdir/ChimericOut/{cell}.sam for begin <= cell <= end.
    Return: integer cell -> set of removed QNAMEs, for filtering gene annotation
    rows in the caller. Only aggregate hotspot tables persist on disk.
    """
    outdir = Path(args.outdir)
    output_dir = outdir / "HotSpotCluster"
    output_dir.mkdir(parents=True, exist_ok=True)
    cells = list(range(args.begin, args.end + 1))
    diagnostics = getattr(args, "diagnostics", True)
    enabled = getattr(args, "hotspot_filter", "on") not in (False, "off")
    if not enabled:
        for name in ("hotspot_clusters.tsv", "adapter_hotspot_families.tsv"):
            (output_dir / name).unlink(missing_ok=True)
        write_summary(output_dir, {"HotspotFilter": "off", "InputCells": len(cells), "FilteredReadNames": 0}, diagnostics)
        return {}
    # Constants match the v2 defaults. Decisions always use all cells together.
    with TemporaryDirectory(prefix=".hotspot-", dir=output_dir) as temporary:
        work = Path(temporary)
        evidence_path = work / "split_evidence.tsv"
        provisional_path = work / "members.tsv"
        intervals, profiles, breakpoints, sources, stats = extract_evidence(
            outdir, "ChimericOut/{cell}.sam", cells, evidence_path, 10,
        )
        all_modes = build_modes(intervals, 5, 5)
        modes = score_modes_by_core_fraction(all_modes, profiles, breakpoints, 150, 0.80)
        collect_members(evidence_path, modes, provisional_path, 150, 6)
        valid, filtered = summarize_modes(modes, provisional_path, output_dir, diagnostics=diagnostics)
        strict_count = sum(map(len, filtered.values()))
        families = build_edge_families(evidence_path, 2, 150, 10, 30, 5, 0.80, 0.80)
        valid_families, edge_keys, new_keys, overlap_keys = summarize_edge_families(
            evidence_path, families, output_dir, filtered, diagnostics=diagnostics,
        )
        jobs = [(sources[str(cell)], filtered.get(str(cell), set())) for cell in cells]
        with ThreadPoolExecutor(max_workers=max(1, args.threads)) as executor:
            removed_records = list(executor.map(lambda item: filter_sam(*item), jobs))
        write_summary(output_dir, {
            "HotspotFilter": "on",
            "InputCells": len(cells),
            "ThreeRecordGroups": stats["groups_3"],
            "OriginalSplitEvidence": stats["split_evidence"],
            "ExactIntervalKeys": len(intervals),
            "CoreModesBeforePurityFilter": len(all_modes),
            "CoreFractionPassingModes": len(modes),
            "FinalHotspotClusters": len(valid),
            "StrictRuleFilteredReadNames": strict_count,
            "EdgeFamilyCandidates": len(families),
            "FinalEdgeFamilies": len(valid_families),
            "EdgeFamilyReadNames": len(edge_keys),
            "EdgeFamilyNewReadNames": len(new_keys),
            "EdgeFamilyOverlappingReadNames": len(overlap_keys),
            "FilteredReadNames": sum(map(len, filtered.values())),
            "RemovedSAMRecords": sum(removed_records),
            "TightWindow": 5,
            "LocalEvidenceWindow": 150,
            "MinCoreFraction": 0.80,
            "FilterBreakpointWindow": 6,
            "MinReads": 5,
            "MinCells": 2,
            "MinPartnerLoci": "DISABLED",
            "MinModalIntervalLength": 10,
            "MaxModalIntervalLength": 20,
            "EdgeFamilyEnabled": True,
            "EdgeBoundaryWindow": 2,
            "EdgeMinAnchor": 10,
            "EdgeMaxAnchor": 30,
            "EdgeSequenceLength": 10,
            "MinEdgeBoundaryFraction": 0.80,
            "MinEdgeSequenceFraction": 0.80,
        }, diagnostics)
    return {int(cell): names for cell, names in filtered.items()}
