"""Reference preparation and explicitly requested FASTQ name conversion."""
from pathlib import Path
import shutil

from .common import executable, run_checked


def attributes(text):
    result = {}
    for part in text.split(";"):
        if part.strip():
            key, value = part.strip().split(None, 1)
            result[key] = value.strip('"')
    return result


def build_star(args):
    args.star_index.mkdir(parents=True, exist_ok=True)
    run_checked([executable("STAR"), "--runMode", "genomeGenerate", "--genomeDir", args.star_index,
                 "--runThreadN", args.threads, "--genomeFastaFiles", args.genome,
                 "--sjdbGTFfile", args.annotation])


def run(args):
    import pysam
    import pyensembl
    directory = args.genome_dir
    directory.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.genome, directory / "ref.fa")
    shutil.copy2(args.annotation, directory / "ref_annot.gtf")
    gtf = directory / "ref_annot.gtf.added"
    with args.annotation.open() as source, gtf.open("w") as output:
        for line in source:
            output.write(line if line.startswith(("#", "chr")) else "chr" + line)
    with gtf.open() as handle, (directory / "exon_probe.hg19.gene.new.bed").open("w") as bed, \
            (directory / "GenePos.txt").open("w") as positions:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields[2] != "gene":
                continue
            info = attributes(fields[8])
            gene = info.get("gene_name") or info.get("gene_id")
            if gene:
                # Preserve the original scFusion BED coordinates and row order.
                bed.write("\t".join([fields[0], fields[3], fields[4], gene]) + "\n")
                positions.write("\t".join([gene, fields[0], fields[3], fields[4],
                                           info.get("gene_type", info.get("gene_biotype", ""))]) + "\n")
    pysam.faidx(str(directory / "ref.fa"))
    pyensembl.Genome(reference_name="GRCH37", annotation_name="my_genome_features",
                    gtf_path_or_url=str(gtf)).index()


def rename(args, restore=False):
    directory = args.file_dir
    mapping = directory / "RenameList.txt"
    changes = []
    if restore:
        with mapping.open() as handle:
            names = [line.rstrip().split("\t") for line in handle]
        for number, name in names:
            for mate in (1, 2):
                for suffix in (".fastq", ".fastq.gz"):
                    source = directory / "{}_{}{}".format(number, mate, suffix)
                    if source.is_file():
                        changes.append((source, directory / "{}_{}{}".format(name, mate, suffix)))
    else:
        if mapping.exists():
            raise FileExistsError(mapping)
        firsts = sorted(path for path in directory.iterdir()
                        if path.name.endswith(("_1.fastq", "_1.fastq.gz")))
        names = []
        for number, first in enumerate(firsts, 1):
            suffix = ".fastq.gz" if first.name.endswith(".gz") else ".fastq"
            name = first.name[:-(len(suffix) + 2)]
            names.append((str(number), name))
            for mate in (1, 2):
                changes.append((directory / "{}_{}{}".format(name, mate, suffix),
                                directory / "{}_{}{}".format(number, mate, suffix)))
    sources = {source for source, _ in changes}
    for source, target in changes:
        if not source.is_file():
            raise FileNotFoundError(source)
        if source != target and (target.exists() or target in sources):
            raise FileExistsError("Name conversion would overwrite {}".format(target))
    for source, target in changes:
        if source != target:
            source.rename(target)
    if not restore:
        with mapping.open("w") as output:
            for number, name in names:
                output.write(number + "\t" + name + "\n")
