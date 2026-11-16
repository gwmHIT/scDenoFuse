"""Original STAR settings, with a single shared thread budget."""
from pathlib import Path

from .common import executable, run_checked


def run(args):
    if args.reuse_from:
        from .reuse import prepare
        return prepare(args)
    for cell in range(args.begin, args.end + 1):
        first = next((args.file_dir / (str(cell) + "_1" + suffix)
                      for suffix in (".fastq.gz", ".fastq")
                      if (args.file_dir / (str(cell) + "_1" + suffix)).is_file()), None)
        if first is None:
            raise FileNotFoundError("Missing R1 FASTQ for cell {}".format(cell))
        second = first.with_name(first.name.replace("_1.fastq", "_2.fastq"))
        if not second.is_file():
            raise FileNotFoundError(second)
        directory = args.outdir / "STARMapping" / str(cell)
        if directory.is_symlink() or (directory.exists() and any(path.is_symlink() for path in directory.iterdir())):
            raise ValueError("Fresh STAR mapping cannot overwrite reused links: {}".format(directory))
        directory.mkdir(parents=True, exist_ok=True)
        command = [executable("STAR"), "--runThreadN", args.threads,
                   "--genomeDir", args.star_index, "--readFilesIn", first, second,
                   "--outSAMtype", "BAM", "SortedByCoordinate",
                   "--chimOutType", "SeparateSAMold", "--outSAMunmapped", "Within", "KeepPairs",
                   "--quantMode", "GeneCounts", "--outFileNamePrefix", str(directory / "human"),
                   "--chimSegmentMin", 12, "--chimJunctionOverhangMin", 8,
                   "--alignSJDBoverhangMin", 10, "--alignMatesGapMax", 100000,
                   "--alignIntronMax", 100000, "--chimSegmentReadGapMax", 3,
                   "--alignSJstitchMismatchNmax", 5, -1, 5, 5,
                   "--alignSplicedMateMapLminOverLmate", 0, "--alignSplicedMateMapLmin", 30,
                   "--chimMultimapScoreRange", 3, "--chimScoreJunctionNonGTAG", -4,
                   "--chimNonchimScoreDropMin", 10, "--peOverlapMMp", 0.1]
        if first.name.endswith(".gz"):
            command += ["--readFilesCommand", "zcat"]
        print("STAR: cell {}".format(cell), flush=True)
        run_checked(command)
