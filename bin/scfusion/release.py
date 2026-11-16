"""Output layout, presentation and STAR cleanup; no scientific logic."""
from contextlib import contextmanager
import datetime
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

VERSION = '2026.11.16'
SCHEMA = 'scDenoFuse-output-v1'
PACKAGE = Path(__file__).resolve().parents[2]
LOG_NAMES = ('threshold_summary.tsv', 'reuse_manifest.tsv', 'Allresult.selection_status.tsv')


def timestamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def write_json(path, data):
    if path.is_symlink():
        raise ValueError('Refusing symlink metadata: {}'.format(path))
    fd, tmp = tempfile.mkstemp(prefix='.pending-', dir=str(path.parent))
    try:
        with os.fdopen(fd, 'w') as handle:
            json.dump(data, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


@contextmanager
def output_context(path, create=False):
    """Manage only an empty directory or a previously registered release output."""
    root = Path(os.path.abspath(str(path.expanduser())))
    if root.is_symlink():
        raise ValueError('Output must be a real directory, not a symlink')
    if not root.exists():
        if not create:
            raise ValueError('Output does not exist; run ReadMapping or Run first')
        root.mkdir(parents=True)
    if not root.is_dir():
        raise ValueError('Output is not a directory')
    info_path = root / 'logs/run_info.json'
    if not info_path.is_file() and any(root.iterdir()):
        raise ValueError('Not a managed scDenoFuse output. Use a new empty output directory; '
                         'historical outputs are not moved or cleaned automatically.')
    for name in ('work', 'logs'):
        directory = root / name
        if directory.is_symlink():
            raise ValueError('Refusing a symlink {} directory'.format(name))
        directory.mkdir(exist_ok=True)
    lock_path = root / 'logs/.run.lock'
    if lock_path.is_symlink():
        raise ValueError('Refusing a symlink lock')
    with lock_path.open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Another scDenoFuse stage is using this output')
        if info_path.is_symlink():
            raise ValueError('Refusing symlink run metadata')
        if info_path.exists():
            info = json.loads(info_path.read_text())
            if info.get('schema') != SCHEMA:
                raise ValueError('Unrecognized output layout version')
        else:
            if not create:
                raise ValueError('Missing managed run metadata')
            info = {'schema': SCHEMA, 'version': VERSION, 'created_utc': timestamp(), 'history': []}
            write_json(info_path, info)
        yield root, info


def save_info(root, info):
    info['updated_utc'] = timestamp()
    write_json(root / 'logs/run_info.json', info)


def bounds_from_info(info):
    begin, end = info.get('begin'), info.get('end')
    if (type(begin) is not int or type(end) is not int or begin < 1 or end < begin):
        raise ValueError('Missing or invalid stored cell range; run ReadMapping first')
    return begin, end


def core_command(args, work=None, bounds=None):
    command = [args.command]
    for name, option in (('fastq', '-f'), ('star_index', '-s'), ('reference', '-d'),
                         ('mappability', '-m'), ('threads', '-t'), ('genome', '-g'),
                         ('annotation', '-a')):
        value = getattr(args, name, None)
        if value is not None:
            command += [option, str(value.expanduser() if isinstance(value, Path) else value)]
    if work is not None:
        command += ['-o', str(work), '-b', str(bounds[0]), '-e', str(bounds[1])]
    return command


def invoke_core(command, log_path=None):
    launch = [sys.executable, '-B', '-u', str(PACKAGE / 'scripts/run_core.py')] + command
    if log_path is None:
        result = subprocess.run(launch)
        if result.returncode:
            raise RuntimeError('Core stage failed with exit code {}'.format(result.returncode))
        return
    if log_path.is_symlink():
        raise ValueError('Refusing a symlink pipeline log')
    with log_path.open('a', encoding='utf-8') as log:
        log.write('\n[{}] {}\n'.format(timestamp(), json.dumps(launch, ensure_ascii=False)))
        log.flush()
        with subprocess.Popen(launch, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              universal_newlines=True, encoding='utf-8', errors='replace',
                              bufsize=1, start_new_session=True) as proc:
            try:
                for line in proc.stdout:
                    log.write(line)
                    log.flush()
                    print(line, end='', flush=True)
                code = proc.wait()
            except BaseException:
                import signal
                try:
                    os.killpg(proc.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                proc.wait()
                raise
        if code:
            raise RuntimeError('Core stage failed with exit code {}; see {}'.format(code, log_path))


def validate_write_paths(root, stage):
    """Core uses open(w): refuse writable symlink outputs before delegating."""
    work = root / 'work'
    directories = {
        'Run': ('STARMapping', 'ChimericOut', 'Expr', 'ChimericExpression', 'Candidates', 'Model', 'HotSpotCluster'),
        'ReadMapping': ('STARMapping',),
        # calculate_rpkm may build missing BAM indexes inside STARMapping.
        'ReadProcessing': ('STARMapping', 'ChimericOut', 'Expr', 'ChimericExpression', 'HotSpotCluster'),
        'FusionCandidate': ('Candidates',),
        'FusionReport': ('Model',),
        'Filter': (),
    }[stage]
    files = list(LOG_NAMES)
    if stage in ('Run', 'FusionReport', 'Filter'):
        files += ['Result.full.txt']
    if stage in ('Run', 'FusionReport'):
        files += ['Allresult.txt']
    if stage in ('Run', 'FusionCandidate'):
        files += ['candidate_cell_counts.tsv']
    for name in files:
        if (work / name).is_symlink():
            raise ValueError('Refusing a symlink writable output: {}'.format(work / name))
    for name in directories:
        directory = work / name
        if directory.is_symlink():
            raise ValueError('Refusing a symlink writable directory: {}'.format(directory))
        if directory.exists():
            for base, dirs, names in os.walk(directory, followlinks=False):
                for entry in dirs + names:
                    if (Path(base) / entry).is_symlink():
                        raise ValueError('Writable stage contains shared symlinks: {}'.format(Path(base) / entry))
    for name in ('Result.full.txt', 'Result.abridged.txt'):
        if (root / name).is_symlink():
            raise ValueError('Refusing a symlink public result: {}'.format(root / name))


def relocate_logs(root):
    for name in LOG_NAMES:
        source, dest = root / 'work' / name, root / 'logs' / name
        if os.path.lexists(source):
            if source.is_symlink() or dest.is_symlink():
                raise ValueError('Refusing a symlink diagnostic file: {}'.format(name))
            os.replace(source, dest)


def threshold_status(root):
    import csv
    path = root / 'logs/threshold_summary.tsv'
    if not path.is_file():
        raise ValueError('Missing threshold summary after reporting')
    with path.open() as handle:
        values = {row['Metric']: row['Value'] for row in csv.DictReader(handle, delimiter='\t')}
    if values.get('ThresholdStatus') not in ('VALID', 'NO_VALID_THRESHOLD'):
        raise ValueError('Invalid threshold status')
    return values['ThresholdStatus']


def clean(root, info):
    """Remove only local STARMapping, after verifying the retained result export."""
    from .presentation import validate_export
    bounds_from_info(info)
    validate_export(root / 'work', root)
    work = root / 'work'
    if root.is_symlink() or work.is_symlink():
        raise ValueError('Refusing cleanup through symlink root/work')
    target = work / 'STARMapping'
    info['cleanup'] = {'status': 'RUNNING', 'target': 'work/STARMapping', 'started_utc': timestamp()}
    save_info(root, info)
    if target.is_symlink():
        target.unlink()
        action = 'UNLINKED_SHARED_DIRECTORY'
    elif target.exists():
        if not target.is_dir() or os.path.ismount(target):
            raise ValueError('STARMapping is not a regular, non-mounted directory')
        if not shutil.rmtree.avoids_symlink_attacks:
            raise RuntimeError('This platform lacks safe directory removal')
        shutil.rmtree(target)
        action = 'REMOVED_LOCAL_DIRECTORY'
    else:
        action = 'ALREADY_ABSENT'
    info['cleanup'].update(status='COMPLETE', action=action, completed_utc=timestamp())
    save_info(root, info)
    print('[Clean] {}: {} (shared source files preserved)'.format(action, target), flush=True)


def run(args, bounds=None):
    if args.command in ('BuildSTARIndex', 'Index'):
        return invoke_core(core_command(args))
    create = args.command in ('ReadMapping', 'Run')
    with output_context(args.output, create=create) as (root, info):
        if bounds is not None:
            if 'begin' in info and bounds_from_info(info) != bounds:
                raise ValueError('Cell range differs from this output; use a new output directory')
            info.update(begin=bounds[0], end=bounds[1], total_cells=bounds[1] - bounds[0] + 1)
        else:
            bounds = bounds_from_info(info)
        if args.command == 'ReadProcessing':
            mapped = root / 'work/STARMapping'
            if not mapped.is_dir() or any(not (mapped / str(i)).is_dir() for i in range(bounds[0], bounds[1] + 1)):
                raise ValueError('STARMapping is unavailable or incomplete. Rerun ReadMapping; '
                                 'use Run --keep-star to retain it next time.')
        recorded_args = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}
        record = {'stage': args.command, 'started_utc': timestamp(), 'status': 'RUNNING',
                  'arguments': recorded_args}
        info.setdefault('history', []).append(record)
        info['status'] = 'RUNNING'
        info['last_command'] = recorded_args
        info['version'] = VERSION
        save_info(root, info)
        phase = 'core'
        try:
            if args.command == 'Clean':
                # A failed earlier scientific stage must be resolved before deletion.
                earlier = [item for item in info['history'][:-1] if item.get('stage') != 'Clean']
                last = earlier[-1] if earlier else {}
                completed = last.get('status') == 'COMPLETE'
                cleanup_only_failure = (last.get('stage') == 'Run' and
                                        last.get('failure_phase') == 'cleanup' and
                                        last.get('export_complete') is True)
                if (last.get('stage') not in ('Run', 'FusionReport', 'Filter') or
                        not (completed or cleanup_only_failure)):
                    raise ValueError('Finish Run, FusionReport or Filter for the current inputs before Clean')
                phase = 'cleanup'
                clean(root, info)
            else:
                validate_write_paths(root, args.command)
                invoke_core(core_command(args, root / 'work', bounds), root / 'logs/pipeline.log')
                relocate_logs(root)
                if args.command in ('Run', 'FusionReport', 'Filter'):
                    phase = 'export'
                    from .presentation import export_results
                    exported = export_results(root / 'work', root)
                    status = threshold_status(root)
                    record['export_complete'] = True
                    info['selected_fusions'] = exported['rows']
                    info['threshold_status'] = status
                    info['last_export_utc'] = timestamp()
                    save_info(root, info)
                    if args.command == 'Run' and not args.keep_star:
                        if status == 'VALID':
                            phase = 'cleanup'
                            clean(root, info)
                        else:
                            print('[Clean] skipped: NO_VALID_THRESHOLD; STARMapping retained for diagnosis', flush=True)
            record.update(status='COMPLETE', finished_utc=timestamp())
            info['status'] = 'COMPLETE'
            save_info(root, info)
            return root / 'Result.full.txt' if args.command in ('Run', 'FusionReport', 'Filter') else None
        except BaseException as error:
            record.update(status='FAILED', error=str(error), failure_phase=phase, finished_utc=timestamp())
            info['status'] = 'FAILED'
            save_info(root, info)
            raise
