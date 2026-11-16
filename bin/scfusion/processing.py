"""Bounded per-cell processing with a global hotspot barrier."""
from array import array
from bisect import bisect_right
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from functools import lru_cache
import multiprocessing
from pathlib import Path
import subprocess
import sys

from .common import BIN, executable, run_checked, unique_file

_mappability = None
_gtf = None
_annotation = None
_outdir = None
_removed = None


def load_mappability(path):
    intervals = {}
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            chrom, start, end, score = line.rstrip().split("\t")[:4]
            if float(score) >= 1:
                intervals.setdefault(chrom, {})[int(start)] = int(end)
    packed = {}
    for chrom in tuple(intervals):
        values = intervals.pop(chrom)
        starts = sorted(values)
        ends = [values[start] for start in starts]
        lower, upper = min(starts[0], min(ends)), max(starts[-1], max(ends))
        # Human genome coordinates fit uint32. Retain unusual signed or larger
        # integer inputs without changing the previous interval semantics.
        if 0 <= lower and upper < 2**32:
            typecode = "I"
        elif -2**63 <= lower and upper < 2**63:
            typecode = "q"
        else:
            typecode = None
        packed[chrom] = (array(typecode, starts), array(typecode, ends)) if typecode else (tuple(starts), tuple(ends))
    return packed


def position_passes(chrom, position, mappability):
    """Preserve the legacy POS and interval-boundary interpretation."""
    if chrom not in mappability:
        return True  # The original filter does not reject an unlisted contig.
    starts, ends = mappability[chrom]
    if position < starts[0] or ends[-1] <= position:
        return False
    if starts[0] <= position <= ends[0]:
        return True
    index = bisect_right(starts, position) - 1
    return position < ends[index]


def passes_group(lines, mappability):
    return _passes_fields([line.split("\t") for line in lines], mappability)


def _passes_fields(records, mappability):
    positions = [fields[3] for fields in records]
    if len(positions) <= 1 or (len(positions) == 2 and abs(int(positions[0]) - int(positions[1])) <= 10):
        return False
    for index, fields in enumerate(records):
        if index == 1 and positions[0] == positions[1]:
            continue
        if index == 2 and (positions[2] == positions[1] or positions[2] == positions[0]):
            continue
        chrom = fields[2] if fields[2].startswith("chr") else "chr" + fields[2]
        if not position_passes(chrom, int(fields[3]), mappability):
            return False
    return True


def filter_mappability(source, target, mappability):
    with Path(source).open(encoding="utf-8", errors="ignore") as handle, Path(target).open("w") as output:
        current, group = None, []
        for line in handle:
            fields = line.split("\t")
            if len(fields) < 5 or line.startswith("@") or "M" in fields[2]:
                continue
            if mappability is None:
                output.write(line)
                continue
            if fields[0] != current:
                if _passes_fields(group, mappability):
                    for values in group:
                        if not values[2].startswith("chr"):
                            values[2] = "chr" + values[2]
                        output.write("\t".join(values))
                current, group = fields[0], [fields]
            else:
                group.append(fields)
        # Historical mappability filtering commits a group when the next QNAME
        # arrives, not at EOF. Preserve this here rather than changing counts.


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def group_contributions(records):
    """Same pre-hotspot fractional expression counts as the v2 side table."""
    result = Counter()
    if len(records) == 2:
        if records[0][0] == records[1][0]:
            return result
        weights = (1, 1)
    elif len(records) == 3:
        if records[0][0] == records[1][0] == records[2][0]:
            return result
        pair = next(((a, b) for a, b in ((0, 1), (0, 2), (1, 2))
                     if records[a][10].upper() in
                     (records[b][10].upper(), reverse_complement(records[b][10].upper()))), None)
        if pair is None:
            return result
        weights = tuple(0.5 if i in pair else 1 for i in range(3))
    else:
        return result
    for record, weight in zip(records, weights):
        genes = dict.fromkeys(gene.strip().split(".")[0] for gene in record[0].split(";") if gene.strip())
        for gene in genes:
            result[gene] += weight
    return result


def count_expression(path, target):
    current, records, counts = None, [], Counter()
    with Path(path).open() as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 11 or line.startswith(("@", "#")):
                continue
            if current is not None and current != fields[1]:
                counts.update(group_contributions(records))
                records = []
            current = fields[1]
            records.append(fields)
    counts.update(group_contributions(records))
    with Path(target).open("w") as output:
        output.write("Gene\tChimericCount\n")
        for gene in sorted(counts):
            output.write("{}\t{:.8g}\n".format(gene, counts[gene]))


def initialize_annotation():
    import pyensembl
    global _annotation
    _annotation = pyensembl.Genome(reference_name="GRCH37", annotation_name="my_genome_features",
                                  gtf_path_or_url=str(_gtf))
    annotation_at_locus.cache_clear()
    # The parent creates the index once. Workers only open the completed DB.


@lru_cache(maxsize=8192)
def annotation_at_locus(chrom, position):
    """Bounded worker-local cache, with the original annotation/fallback order."""
    genes = _annotation.gene_names_at_locus(contig=chrom, position=position)
    if not genes:
        genes = _annotation.gene_names_at_locus(contig=chrom, position=position + 20)
    if not genes:
        genes = _annotation.gene_names_at_locus(contig=chrom, position=position - 20)
    return ";".join(genes)


def prepare_cell(cell):
    source = unique_file(_outdir / "STARMapping" / str(cell), "*Chimeric.out.sam")
    sam = _outdir / "ChimericOut" / "{}.sam".format(cell)
    geneanno = sam.with_name("{}_geneanno.sam".format(cell))
    filter_mappability(source, sam, _mappability)
    with sam.open() as handle, geneanno.open("w") as output:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.split("\t")
            chrom, position = fields[2], int(fields[3])
            output.write(annotation_at_locus(chrom, position) + "\t" + line)
    count_expression(geneanno, _outdir / "ChimericExpression" / "{}.gene_counts.tsv".format(cell))
    return cell


def support_cell(cell):
    path = _outdir / "ChimericOut" / "{}_geneanno.sam".format(cell)
    removed = _removed.get(cell, set())
    if removed:
        temporary = path.with_name(path.name + ".hotspot.tmp")
        with path.open() as handle, temporary.open("w") as output:
            for line in handle:
                if line.split("\t", 2)[1] not in removed:
                    output.write(line)
        temporary.replace(path)
    run_checked([sys.executable, BIN / "FindFusionSupport.py", path,
                 _outdir / "ChimericOut" / "{}_FusionSupport.txt".format(cell)])
    return cell


def calculate_rpkm(cell, outdir, bed):
    bam = unique_file(outdir / "STARMapping" / str(cell), "*Aligned.sortedByCoord.out.bam")
    index = Path(str(bam) + ".bai")
    if not index.is_file():
        if index.is_symlink():
            raise ValueError("Cannot write through a missing reused BAM-index link: {}".format(index))
        run_checked([executable("samtools"), "index", bam, index])
    stats = run_checked([executable("samtools"), "idxstats", bam], stdout=subprocess.PIPE, text=True).stdout
    total = sum(int(line.split("\t")[2]) for line in stats.splitlines())
    if total == 0:
        raise ValueError("No mapped reads: {}".format(bam))
    command = [executable("bedtools"), "multicov", "-bams", str(bam), "-bed", str(bed)]
    target = outdir / "Expr" / "{}.rpkm.txt".format(cell)
    temporary = target.with_name(target.name + ".tmp")
    with temporary.open("w") as output, subprocess.Popen(command, stdout=subprocess.PIPE, text=True) as process:
        for line in process.stdout:
            fields = line.rstrip().split("\t")
            length, count = int(fields[2]) - int(fields[1]), int(fields[4])
            value = 1e9 * count / (length * total) if length >= 1 else 0
            output.write("{}\t{}\t{:.15g}\n".format(fields[3], count, value))
        if process.wait():
            raise subprocess.CalledProcessError(process.returncode, command)
    temporary.replace(target)
    return cell


def run(args):
    import pyensembl
    from . import hotspot
    global _mappability, _gtf, _outdir, _removed
    _outdir, _gtf = args.outdir, args.genome_dir / "ref_annot.gtf.added"
    cells = list(range(args.begin, args.end + 1))
    for directory in ("ChimericOut", "ChimericExpression", "Expr"):
        (_outdir / directory).mkdir(parents=True, exist_ok=True)
    if args.reuse_from and not getattr(args, "_reuse_prepared", False):
        from .reuse import prepare
        prepare(args)
    _mappability = None if args.no_mappability else load_mappability(args.mappability)
    data = pyensembl.Genome(reference_name="GRCH37", annotation_name="my_genome_features",
                           gtf_path_or_url=str(_gtf))
    data.index()
    context = multiprocessing.get_context("fork")
    workers = min(args.threads, len(cells))
    print("ReadProcessing: {} cells, {} workers; annotation and expression before hotspot".format(len(cells), workers), flush=True)
    with ProcessPoolExecutor(max_workers=workers, mp_context=context, initializer=initialize_annotation) as pool:
        list(pool.map(prepare_cell, cells))
    _mappability = None
    del data
    _removed = hotspot.run(args)
    print("ReadProcessing: fusion support after global hotspot filtering", flush=True)
    with ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
        list(pool.map(support_cell, cells))
    _removed = None
    if not args.reuse_from:
        bed = args.genome_dir / "exon_probe.hg19.gene.new.bed"
        print("RPKM: {} concurrent bedtools jobs".format(workers), flush=True)
        with ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(lambda cell: calculate_rpkm(cell, _outdir, bed), cells))
