#!/usr/bin/env python3
"""Public scDenoFuse interface; scientific computation is unchanged."""
import argparse
from pathlib import Path
import re
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent / 'bin'))
from scfusion import release


def cell_range(directory):
    mates = {}
    for path in directory.iterdir():
        match = re.fullmatch(r'([1-9][0-9]*)_([12])\.fastq(\.gz)?', path.name)
        if match and path.is_file():
            cell, mate = map(int, match.group(1, 2))
            if (cell, mate) in mates:
                raise ValueError('Multiple FASTQ files for cell {}'.format(cell))
            mates[cell, mate] = path
    cells = sorted({cell for cell, mate in mates})
    if not cells:
        raise ValueError('Expected paired numbered FASTQs such as 1_1.fastq.gz and 1_2.fastq.gz')
    if cells != list(range(cells[0], cells[-1] + 1)):
        raise ValueError('Cell numbers must be contiguous; identities are not renumbered')
    for cell in cells:
        if (cell, 1) not in mates or (cell, 2) not in mates:
            raise ValueError('Missing FASTQ mate for cell {}'.format(cell))
        if mates[cell, 1].suffixes != mates[cell, 2].suffixes:
            raise ValueError('FASTQ mates must use the same compression')
    return cells[0], cells[-1]


def build_parser():
    parser = argparse.ArgumentParser(description='scDenoFuse: compact results and per-cell evidence',
        epilog='Full run: scDenoFuse.py -i FASTQ -s STAR_INDEX -d REFERENCE -m MAPPABILITY -o OUTPUT')
    parser.add_argument('--version', action='version', version='scDenoFuse ' + release.VERSION)
    commands = parser.add_subparsers(dest='command', required=True)
    for name in ('Run', 'ReadMapping', 'ReadProcessing', 'FusionCandidate', 'FusionReport',
                 'Filter', 'Clean', 'BuildSTARIndex', 'Index'):
        stage = commands.add_parser(name, help={
            'Run': 'Full pipeline; removes work/STARMapping after successful export',
            'Clean': 'Remove local work/STARMapping; preserve shared symlink targets',
        }.get(name, name))
        if name not in ('BuildSTARIndex', 'Index'):
            stage.add_argument('-o', '--output', type=Path, required=True)
        if name in ('Run', 'ReadMapping'):
            stage.add_argument('-i', '-f', '--fastq', type=Path, required=True)
        if name in ('Run', 'ReadMapping', 'BuildSTARIndex'):
            stage.add_argument('-s', '--star-index', type=Path, required=True)
        if name not in ('ReadMapping', 'BuildSTARIndex', 'Clean'):
            stage.add_argument('-d', '--reference', type=Path, required=True)
        if name in ('Run', 'ReadProcessing'):
            stage.add_argument('-m', '--mappability', type=Path, required=True)
        if name in ('BuildSTARIndex', 'Index'):
            stage.add_argument('-g', '--genome', type=Path, required=True)
            stage.add_argument('-a', '--annotation', type=Path, required=True)
        if name not in ('Index', 'Clean'):
            stage.add_argument('-t', '--threads', type=int, help='Thread budget (default: 16)')
        if name == 'Run':
            stage.add_argument('--keep-star', action='store_true',
                help='Retain STARMapping for IGV, reuse or rerunning ReadProcessing')
    return parser


def main(argv=None):
    values = list(sys.argv[1:] if argv is None else argv)
    if values and values[0].startswith('-') and values[0] not in ('-h', '--help', '--version'):
        values.insert(0, 'Run')
    parser = build_parser()
    args = parser.parse_args(values)
    if getattr(args, 'threads', None) is not None and args.threads < 1:
        parser.error('threads must be positive')
    try:
        bounds = cell_range(args.fastq.expanduser()) if hasattr(args, 'fastq') else None
        return release.run(args, bounds)
    except (OSError, ValueError, RuntimeError) as error:
        parser.exit(1, 'scDenoFuse: {}\n'.format(error))


if __name__ == '__main__':
    main()
