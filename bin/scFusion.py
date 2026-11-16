#!/usr/bin/env python3
"""The sole public entry point for lightweight scFusion v3."""
import argparse
from pathlib import Path
import time


def build_parser():
    parser = argparse.ArgumentParser(description="scFusion v3: hotspot + original ZINB + dynamic threshold + baseline1 + baseline3 (R2/R3) + orientation + recurrent background filters")
    commands = parser.add_subparsers(dest="command", required=True)

    def command(name, description, output=False, genome=False, cells=False, prefix=False):
        p = commands.add_parser(name, help=description)
        if output:
            p.add_argument("-o", "--OutDir", "--outdir", dest="outdir", type=Path, required=True)
        if genome:
            p.add_argument("-d", "--GenomeDir", "--genome-dir", dest="genome_dir", type=Path, required=True)
        if cells:
            p.add_argument("-b", "--Begin", "--begin", dest="begin", type=int, required=True)
            p.add_argument("-e", "--End", "--end", dest="end", type=int, required=True)
        if prefix:
            p.add_argument("-p", "--Prefix", "--prefix", dest="prefix", default=".")
        return p

    def threads(p):
        p.add_argument("-t", "--Thread", "--threads", dest="threads", type=int, default=16)

    def mapping(p):
        p.add_argument("-f", "--FileDir", dest="file_dir", type=Path)
        p.add_argument("-s", "--STARIndex", dest="star_index", type=Path)
        p.add_argument("--reuse-from", type=Path, help="Reuse only STAR and RPKM from this result root")

    def processing(p):
        group = p.add_mutually_exclusive_group()
        group.add_argument("-m", "--Mappability", dest="mappability", type=Path)
        group.add_argument("-M", "--NoMappabilityFilter", dest="no_mappability", action="store_true")
        p.add_argument("--hotspot-filter", "--HotSpotFilter", choices=("on", "off"), default="on")

    def reporting(p):
        from scfusion.background import DEFAULT_BLACKLIST
        p.add_argument("-v", "--PvalueCutoff", "--fdr", dest="fdr_cutoff", type=float, default=0.05)
        p.add_argument("--exon-window", type=int, default=100)
        p.add_argument("--rscript", default="Rscript")
        p.add_argument("--baseline3", choices=("on", "off"), default="on",
                       help="R2 reference-side homology and R3 lncRNA readthrough filters after baseline1 (default: on)")
        p.add_argument("--orientation-filter", choices=("on", "off"), default="on",
                       help="Remove known 5prime--5prime / 3prime--3prime retained portions after baseline3; keep unknowns (default: on)")
        p.add_argument("--background-filter", choices=("on", "off"), default="on",
                       help="Exact recurrent background breakpoint-pair filter after orientation filtering (default: on)")
        p.add_argument("--background-blacklist", type=Path, default=DEFAULT_BLACKLIST,
                       help="Arriba-format resource with a recurrent breakpoints section (default: bundled GRCh38 v2.5.1)")

    p = command("BuildSTARIndex", "Build a STAR genome index")
    p.add_argument("-s", "--STARIndex", dest="star_index", type=Path, required=True)
    p.add_argument("-g", "--Genome", dest="genome", type=Path, required=True)
    p.add_argument("-a", "--Annotation", dest="annotation", type=Path, required=True)
    threads(p)
    p = command("Index", "Prepare scFusion reference files", genome=True)
    p.add_argument("-g", "--Genome", dest="genome", type=Path, required=True)
    p.add_argument("-a", "--Annotation", dest="annotation", type=Path, required=True)
    for name in ("Rename", "RestoreName"):
        p = command(name, "Explicitly rename or restore paired FASTQ files")
        p.add_argument("-f", "--FileDir", dest="file_dir", type=Path, required=True)
    p = command("ReadMapping", "Original paired-end STAR mapping or explicit reuse", output=True, cells=True)
    p.add_argument("-d", "--GenomeDir", dest="genome_dir", type=Path, help="Required with --reuse-from")
    mapping(p)
    threads(p)
    p = command("ReadProcessing", "Mappability, global hotspots, support and expression", output=True, genome=True, cells=True)
    p.add_argument("--reuse-from", type=Path)
    p.add_argument("-f", "--FileDir", dest="file_dir", type=Path, help="Optional FASTQ identity check with reuse")
    processing(p)
    threads(p)
    p = command("FusionCandidate", "Candidate clustering and essential per-cell counts", output=True, genome=True, cells=True, prefix=True)
    threads(p)
    p = command("FusionReport", "Original ZINB, dynamic threshold and post-filters including orientation/background", output=True, genome=True, cells=True, prefix=True)
    p.add_argument("-f", "--FileDir", dest="file_dir", type=Path, help="Accepted for old command compatibility; not needed")
    reporting(p)
    threads(p)
    p = command("Filter", "Reapply threshold and post-filters to existing v3 Allresult without refitting", output=True, genome=True, cells=True, prefix=True)
    reporting(p)
    threads(p)
    p = command("Run", "Run the complete lightweight pipeline", output=True, genome=True, cells=True, prefix=True)
    mapping(p)
    processing(p)
    reporting(p)
    threads(p)
    return parser


def stage(name, function, args):
    started = time.monotonic()
    print("[{}] start".format(name), flush=True)
    result = function(args)
    print("[{}] complete ({:.1f}s)".format(name, time.monotonic() - started), flush=True)
    return result


def main(argv=None, *, diagnostics=True):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.diagnostics = diagnostics
    for key, value in vars(args).items():
        if isinstance(value, Path):
            setattr(args, key, value.expanduser().resolve())
    if hasattr(args, "begin") and not 1 <= args.begin <= args.end:
        parser.error("Require 1 <= begin <= end")
    if hasattr(args, "threads") and args.threads < 1:
        parser.error("threads must be positive")
    if args.command in ("ReadMapping", "Run"):
        if args.reuse_from:
            if args.genome_dir is None:
                parser.error("--reuse-from requires -d GenomeDir to validate expression gene order")
        elif args.file_dir is None or args.star_index is None:
            parser.error("Mapping without reuse requires -f FileDir and -s STARIndex")
    if args.command in ("ReadProcessing", "Run") and not args.no_mappability and args.mappability is None:
        parser.error("Provide -m Mappability for this genome, or explicitly disable it with -M")
    if hasattr(args, "fdr_cutoff") and not 0 < args.fdr_cutoff < 1:
        parser.error("Require 0 < fdr < 1")
    if hasattr(args, "exon_window") and args.exon_window < 0:
        parser.error("exon-window must be nonnegative")
    if getattr(args, "background_filter", "off") == "on":
        from scfusion.background import validate_resource
        try:
            validate_resource(args.background_blacklist, args.genome_dir)
        except (OSError, ValueError) as error:
            parser.error(str(error))
    from scfusion import indexing
    if args.command == "BuildSTARIndex":
        return indexing.build_star(args)
    if args.command == "Index":
        return indexing.run(args)
    if args.command in ("Rename", "RestoreName"):
        return indexing.rename(args, restore=args.command == "RestoreName")
    if args.command == "ReadMapping":
        from scfusion import mapping
        return stage(args.command, mapping.run, args)
    if args.command == "ReadProcessing":
        from scfusion import processing
        return stage(args.command, processing.run, args)
    if args.command == "FusionCandidate":
        from scfusion import candidates
        return stage(args.command, candidates.run, args)
    if args.command == "Filter":
        from scfusion import filtering
        return stage(args.command, filtering.run, args)
    from scfusion import model, filtering
    if args.command == "Run":
        from scfusion import mapping, processing, candidates
        stage("ReadMapping", mapping.run, args)
        stage("ReadProcessing", processing.run, args)
        stage("FusionCandidate", candidates.run, args)
    stage("ZINB", model.run, args)
    return stage("Filter", filtering.run, args)


if __name__ == "__main__":
    main()
