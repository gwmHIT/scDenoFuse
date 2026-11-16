"""Candidate clustering and essential per-cell counts, without read tracking."""
import ast
import csv
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import pysam

from .common import BIN, prefix_value, read_tsv, run_checked, write_tsv
from .orientation import annotate, read_gene_strands


HEADERS = (
    'Gene1', 'Gene2', 'CellSupport', 'SplitReadCountsByCell',
    'SpanningSupportingCellCount', 'SpanningReadCountsByCell',
    'Chrom1', 'Chrom2', 'Breakpoint1', 'Breakpoint2', 'ClusterIndex',
    'GC1', 'GC2', 'Gene1ExpressionBySplitCells', 'Gene2ExpressionBySplitCells',
    'CandidateID', 'SupportingCells', 'PairSpanningCells', 'PairSpanningCounts',
    'direction1', 'direction2', 'DirectionCountPP', 'DirectionCountPM',
    'DirectionCountMP', 'DirectionCountMM', 'DirectionAmbiguousCount', 'DirectionStatus',
    'strand1', 'strand2')


def background_rows(rows):
    """Reference FilterChiDist background selection and its stable ordering."""
    eligible = []
    for row in rows:
        counts = [int(value) for value in ast.literal_eval(row['SplitReadCountsByCell'])
                  if int(value) >= 1]
        if len(counts) <= 1 or sum(counts) / len(counts) > 5:
            continue
        item = dict(row, CellSupport=str(len(counts)), SplitReadCountsByCell=str(counts))
        eligible.append(item)
    eligible.sort(key=lambda row: len(row['SplitReadCountsByCell']), reverse=True)
    return eligible[min(10, len(eligible) // 50):]


def cell_count_rows(rows):
    """Sparse union of split-support cells and gene-pair spanning cells."""
    for row in rows:
        split = dict(zip(ast.literal_eval(row['SupportingCells']),
                         ast.literal_eval(row['SplitReadCountsByCell'])))
        spanning = dict(zip(ast.literal_eval(row['PairSpanningCells']),
                            ast.literal_eval(row['PairSpanningCounts'])))
        for cell in sorted(split.keys() | spanning.keys()):
            yield {'CandidateID': row['CandidateID'], 'Cell': cell,
                   'SplitCount': int(split.get(cell, 0)),
                   'SpanningCount': int(spanning.get(cell, 0))}


def run(args):
    outdir, genome_dir = Path(args.outdir), Path(args.genome_dir)
    prefix = prefix_value(args.prefix)
    directory = outdir / 'Candidates'
    directory.mkdir(parents=True, exist_ok=True)
    chimeric = outdir / 'ChimericOut'
    expression = outdir / 'Expr'
    env = {'SCFUSION_CHIMERIC_EXPR_DIR': str((outdir / 'ChimericExpression').resolve())}
    candidate_path = directory / (prefix + 'ChiDist.tsv')
    with TemporaryDirectory(prefix='gc_', dir=directory) as workspace:
        workspace = Path(workspace)
        score_path, gc_path = workspace / 'FusionScore.txt', workspace / 'GC.txt'
        strand_path = workspace / 'GeneStrands.tsv'
        reference = genome_dir / 'ref.fa'
        if not Path(str(reference) + '.fai').is_file():
            # faidx must not modify an external/reused GenomeDir implicitly.
            indexed = workspace / 'ref.fa'
            indexed.symlink_to(reference.resolve())
            pysam.faidx(str(indexed))
            reference = indexed
        with score_path.open('w') as handle:
            run_checked([sys.executable, BIN / 'FusionScore.py', chimeric,
                         args.begin, args.end, expression], stdout=handle, env=env)
        with gc_path.open('w') as handle:
            run_checked([sys.executable, BIN / 'FindHomoPattern_RAM.py', score_path,
                         reference, genome_dir / 'ref_annot.gtf.added', strand_path], stdout=handle)
        with candidate_path.open('w') as handle:
            run_checked([sys.executable, BIN / 'FindChiDist.py', chimeric,
                         args.begin, args.end, expression, gc_path, prefix],
                        stdout=handle, env=env)
        rows = annotate(read_tsv(candidate_path), read_gene_strands(strand_path))
    write_tsv(candidate_path, HEADERS, rows)
    write_tsv(directory / (prefix + 'ChiDist_background.tsv'), HEADERS, background_rows(rows))
    write_tsv(outdir / (prefix + 'candidate_cell_counts.tsv'),
              ('CandidateID', 'Cell', 'SplitCount', 'SpanningCount'), cell_count_rows(rows))
    return candidate_path
