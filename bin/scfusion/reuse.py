"""Read-only STAR/RPKM reuse, independent of detection and filtering."""
import hashlib
import math
from pathlib import Path
import shlex

from .common import unique_file, write_tsv


def link(source, target):
    source, target = Path(source).resolve(), Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.is_symlink() and target.resolve() == source:
        return
    if target.exists() or target.is_symlink():
        raise FileExistsError("Reuse target already exists: {}".format(target))
    target.symlink_to(source)


def validate_expression(path, genes):
    with Path(path).open() as handle:
        rows = [line.rstrip("\r\n").split("\t") for line in handle]
    if len(rows) != len(genes):
        raise ValueError("RPKM row count differs from exon BED: {}".format(path))
    for gene, row in zip(genes, rows):
        if len(row) != 3 or row[0] != gene:
            raise ValueError("RPKM gene order differs from exon BED: {}".format(path))
        if int(row[1]) < 0 or not math.isfinite(float(row[2])) or float(row[2]) < 0:
            raise ValueError("Invalid expression value in {}".format(path))


def validate_cell_identity(directory, cell, file_dir):
    """An optional FASTQ input must identify the same cell as STAR's log."""
    source_reads = None
    for log in sorted(directory.glob("*Log.out")):
        with log.open() as handle:
            for line in handle:
                values = shlex.split(line)
                if values and values[0] == "readFilesIn" and len(values) >= 3:
                    source_reads = tuple(Path(value).resolve() for value in values[1:3])
                    break
        if source_reads:
            break
    if file_dir is None:
        return str(source_reads[0]) if source_reads else "", "original_numeric_cell_id"
    if source_reads is None:
        raise ValueError("Cannot verify supplied FASTQ identity without STAR Log.out: cell {}".format(cell))
    for mate, source in enumerate(source_reads, 1):
        current = next((file_dir / "{}_{}{}".format(cell, mate, suffix)
                        for suffix in (".fastq.gz", ".fastq")
                        if (file_dir / "{}_{}{}".format(cell, mate, suffix)).is_file()), None)
        if current is None or source != current.resolve():
            raise ValueError("FASTQ identity differs from reused STAR cell {} mate {}".format(cell, mate))
    return str(source_reads[0]), "exact_resolved_fastq_pair"


def prepare(args):
    """Validate the complete selection, then link only STAR files and Expr."""
    source, output = args.reuse_from.resolve(), args.outdir.resolve()
    if source == output:
        raise ValueError("Reuse source and output must be different directories")
    bed = args.genome_dir / "exon_probe.hg19.gene.new.bed"
    with bed.open() as handle:
        genes = [line.rstrip().split("\t")[3] for line in handle if line.strip()]
    provenance = source / "Expr/reuse_provenance.tsv"
    if provenance.is_file():
        metadata = dict(line.rstrip().split("\t", 1) for line in provenance.read_text().splitlines()[1:])
        expected = metadata.get("ExonBED_SHA256")
        if expected and hashlib.sha256(bed.read_bytes()).hexdigest() != expected:
            raise ValueError("Reuse source was generated with a different exon BED")
    selected = []
    for cell in range(args.begin, args.end + 1):
        directory = source / "STARMapping" / str(cell)
        bam = unique_file(directory, "*Aligned.sortedByCoord.out.bam")
        sam = unique_file(directory, "*Chimeric.out.sam")
        expr = source / "Expr" / "{}.rpkm.txt".format(cell)
        validate_expression(expr, genes)
        read1, identity = validate_cell_identity(directory, cell, getattr(args, "file_dir", None))
        selected.append((cell, directory, bam, sam, expr, read1, identity))
    rows = []
    for cell, directory, bam, sam, expr, read1, identity in selected:
        for original in sorted(directory.iterdir()):
            if original.is_file():
                link(original, output / "STARMapping" / str(cell) / original.name)
        link(expr, output / "Expr" / expr.name)
        rows.append(dict(Cell=cell, SourceBAM=str(bam.resolve()), SourceChimericSAM=str(sam.resolve()),
                         SourceExpression=str(expr.resolve()), SourceR1=read1, CellIdentityCheck=identity))
    write_tsv(output / "reuse_manifest.tsv",
              ["Cell", "SourceBAM", "SourceChimericSAM", "SourceExpression", "SourceR1", "CellIdentityCheck"], rows)
    args._reuse_prepared = True
    print("Reused STAR and RPKM for {} cells from {}".format(len(rows), source), flush=True)
