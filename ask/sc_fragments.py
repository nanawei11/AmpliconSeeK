"""Optional sinto fragments generation using the recorded CB/tagtoname workflow."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

import pysam


def report_missing_fragments(bam_path):
    """Explain the required input and show an explicit BAM preprocessing command."""
    import shlex
    bam = shlex.quote(str(bam_path)) if bam_path is not None else 'sample.bam'
    script = shlex.quote(str(Path(__file__).resolve()))
    print(
        '[single-cell] No fragments file was supplied. Cell assignment was not started.\n'
        'To enable cell assignment, derive fragments from a coordinate-sorted, '
        'paired-end BAM with CB barcode tags. If the tools are not installed, run:\n'
        '  conda create -n sinto-fragments -c conda-forge -c bioconda '
        'python=3.10 "setuptools<81" sinto samtools coreutils\n'
        'Generate fragments separately (this can take time):\n'
        f'  conda run -n sinto-fragments python {script} -i {bam} -o preprocessing/sample -p 4\n'
        'Then rerun ASK with --fragments preprocessing/sample/fragments.tsv.gz.\n'
        'No software is installed and no fragments are generated automatically. '
        'See README.md: Single-cell fragments from BAM. On Windows, use WSL2/Linux '
        'for this workflow.', file=sys.stderr, flush=True)


#------------------------------------------------------------------------------#
def add_fragment_generation_arguments(parser):
    parser.add_argument('--make-fragments', action='store_true',
                        help='Legacy inline preprocessing (can be slow). Prefer ask-make-fragments separately, then --fragments FILE.')
    parser.add_argument('--sinto-bin', default='sinto', help='sinto executable or absolute path.')
    parser.add_argument('--samtools-bin', default='samtools', help='samtools executable or absolute path.')
    parser.add_argument('--fragments-threads', type=int, default=4, help='sinto fragments workers (default: 4).')
    parser.add_argument('--fragments-mapq', type=int, default=20, help='sinto minimum MAPQ (default: 20).')
    parser.add_argument('--fragments-use-chrom', default='(?i)^chr', help='sinto chromosome regex; default matches historical chr-prefixed workflow.')
    parser.add_argument('--fragments-dir', default=None, help='Directory for generated fragments and provenance; completed matching output is validated and reused.')


#------------------------------------------------------------------------------#
def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


class FragmentDependencyError(RuntimeError):
    """Required fragment-generation executables are unavailable."""


def check_fragment_dependencies(sinto='sinto', samtools='samtools'):
    """Resolve every required executable before reading BAM or creating output."""
    executables, missing = {}, []
    for name, value in [('sinto', sinto), ('samtools', samtools), ('sort', 'sort')]:
        resolved = shutil.which(str(value))
        if resolved is None:
            missing.append(f"{name} (requested: {value})")
        else:
            executables[name] = str(Path(resolved).resolve())
    if missing:
        raise FragmentDependencyError(
            "Missing fragment-generation software: " + ", ".join(missing) + "\n"
            "No fragment generation was started.\n"
            "Example installation with Conda (run manually):\n"
            '  conda create -n sinto-fragments -c conda-forge -c bioconda python "setuptools<82" sinto samtools coreutils\n'
            "  conda activate sinto-fragments\n"
            "  command -v sinto samtools sort\n"
            "Keep running ASK with its existing Python environment. If the tools are installed "
            "elsewhere, add their bin directory to PATH or provide --sinto-bin and --samtools-bin.\n"
            "Installation documentation: https://timoast.github.io/sinto/installation.html"
        )
    return executables


#------------------------------------------------------------------------------#
def generate_fragments(bam_path, outdir, sinto='sinto', samtools='samtools', threads=4,
                       min_mapq=20, use_chrom='(?i)^chr'):
    """Generate and validate BGZF fragments; never modify BAM or overwrite output.

    Matches the 20260910 CB/tagtoname profile, not the older filtered-BAM,
    MAPQ30/collapse_within author-union profile. Temporary BAMs are removed.
    """
    if threads < 1 or min_mapq < 0:
        raise ValueError('fragments threads must be >=1 and MAPQ must be >=0')
    executables = check_fragment_dependencies(sinto, samtools)
    bam_path = Path(bam_path).resolve(strict=True)
    with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
        if bam.header.to_dict().get('HD', {}).get('SO') != 'coordinate':
            raise ValueError('sinto generation requires a coordinate-sorted BAM')
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=False)
    final = outdir / 'fragments.tsv.gz'
    stat = bam_path.stat()
    started_clock = time.monotonic()
    record = dict(status='running', source=str(bam_path), source_size=stat.st_size,
                  source_mtime_ns=stat.st_mtime_ns, output=str(final), executables=executables,
                  profile='CB_tagtoname_20260910', threads=threads, min_mapq=min_mapq,
                  min_distance=10, max_distance=5000, shift_plus=4, shift_minus=-5,
                  use_chrom=use_chrom, chunksize=500000, collapse_within=False,
                  pysam_version=pysam.__version__, commands=[], started=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
    method = outdir / '\u505a\u6cd5.md'
    method.write_text(f'# sinto fragment generation\n\nStatus: running; validation is not complete.\nGenerator: {Path(__file__).resolve()}\n'
                      f'Input BAM: {bam_path}\nOutput: {final} and .tbi\n'
                      f'Full parameters, commands, executables, and input size/mtime: {outdir}/manifest.json\n'
                      f'Log: {outdir}/commands.log\n'
                      'Method: retain CB-tagged reads, run sinto tagtoname then sinto fragments; sort by chromosome/start with LC_ALL=C, '
                      'then use pysam BGZF compression and tabix indexing. Column 5 contains sinto collapsed support counts. No additional whitelist or blacklist filtering. '
                      'Uses the 20260910 CB/tagtoname convention, which differs from the legacy author-union MAPQ30/collapse_within convention.\n'
                      'Validation: this generator checks five columns, coordinates, barcodes, counts, sorting, and the tabix index; no independent validator. '
                      f'Downstream: {Path(__file__).resolve().parent / "sc_assignment.py"} reads the output for coverage/Q1. Temporary BAM/BED files are removed on exit.\n')
    def save():
        record['elapsed_seconds'] = time.monotonic() - started_clock
        (outdir / 'manifest.json').write_text(json.dumps(record, indent=2))
    env = dict(os.environ, LC_ALL='C', PYTHONNOUSERSITE='1')
    def run(command, stdout=None):
        command = [str(x) for x in command]
        record['commands'].append(command)
        save()
        with (outdir / 'commands.log').open('a') as log:
            log.write(json.dumps(command) + '\n'); log.flush()
            subprocess.run(command, stdout=stdout if stdout is not None else log, stderr=log,
                           check=True, env=env, cwd=outdir)
    save()
    try:
        with tempfile.TemporaryDirectory(prefix='work_', dir=outdir) as temp:
            work = Path(temp)
            cb_bam, tagged = work/'CB.bam', work/'tagged.bam'
            unsorted, sorted_bed = work/'unsorted.bed', work/'sorted.bed'
            run([executables['sinto'], '--version'])
            run([executables['samtools'], '--version'])
            run([executables['samtools'], 'quickcheck', bam_path])
            run([executables['samtools'], 'view', '-h', '-d', 'CB', '-b', '-o', cb_bam, bam_path])
            run([executables['sinto'], 'tagtoname', '-b', cb_bam, '-O', 'b', '--tag', 'CB', '-o', tagged])
            run([executables['samtools'], 'quickcheck', tagged])
            run([executables['samtools'], 'index', tagged])
            run([executables['sinto'], 'fragments', '-b', tagged, '-p', threads, '-f', unsorted,
                 '--barcode_regex', '[^:]*', '--min_mapq', min_mapq, '--min_distance', 10,
                 '--max_distance', 5000, '--shift_plus', 4, '--shift_minus', -5,
                 '--chunksize', 500000, '--use_chrom', use_chrom])
            run([executables['sort'], '-T', work, '-S', '2G', '-k1,1', '-k2,2n', '-o', sorted_bed, unsorted])
            rows, support, previous, barcodes = 0, 0, None, set()
            with sorted_bed.open() as handle:
                for line in handle:
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) != 5:
                        raise ValueError('Expected five-column sinto fragments')
                    chrom, start, end, barcode, count = fields
                    start, end, count = int(start), int(end), int(count)
                    key = (chrom, start)
                    if not chrom or not barcode or ':' in barcode or start < 0 or end <= start or count < 1:
                        raise ValueError('Invalid fragment coordinates/barcode/count')
                    if previous is not None and key < previous:
                        raise ValueError('Fragments are not sorted')
                    rows += 1; support += count; barcodes.add(barcode); previous = key
            if rows == 0:
                raise ValueError('No fragments generated; check CB tags, pairing, MAPQ and chromosome naming')
            pysam.tabix_compress(str(sorted_bed), str(final), force=False)
            pysam.tabix_index(str(final), preset='bed', force=False)
            with gzip.open(final, 'rb') as handle:
                while handle.read(8 * 1024 * 1024):
                    pass
            with pysam.TabixFile(str(final)) as indexed:
                contigs = list(indexed.contigs)
            current = bam_path.stat()
            if (current.st_size, current.st_mtime_ns) != (stat.st_size, stat.st_mtime_ns):
                raise RuntimeError('Input BAM changed during generation')
            record.update(status='completed', fragment_rows=rows, cells_with_fragments=len(barcodes),
                          sum_column5=support, contigs=contigs, sha256=file_sha256(final),
                          index_sha256=file_sha256(str(final) + '.tbi'))
        save()
        with method.open('a') as handle:
            handle.write(f'\nCompleted: {rows} rows, {len(barcodes)} barcodes; format, compression, and index checks passed. '
                         'Exact equivalence to Cell Ranger or historical author-union fragments has not been verified.\n')
        return str(final)
    except Exception as error:
        record.update(status='failed', error=str(error)); save()
        with method.open('a') as handle:
            handle.write(f'\nFailed; incomplete: {error}\n')
        raise


#------------------------------------------------------------------------------#
class FragmentReuseError(RuntimeError):
    """Existing fragments cannot safely be reused."""


#------------------------------------------------------------------------------#
def reuse_fragments(bam_path, directory, min_mapq, use_chrom):
    """Reuse only a completed output matching the input, settings and hashes."""
    outdir = Path(directory).resolve()
    manifest = outdir / 'manifest.json'
    final = outdir / 'fragments.tsv.gz'
    index = Path(str(final) + '.tbi')
    if not outdir.is_dir():
        raise FragmentReuseError(f'Not a directory: {outdir}')
    try:
        record = json.loads(manifest.read_text())
    except (OSError, ValueError) as error:
        raise FragmentReuseError(f'Cannot read a valid manifest.json: {error}') from error
    if not isinstance(record, dict):
        raise FragmentReuseError('Invalid manifest.json: the top-level value must be an object')
    source = Path(bam_path).resolve(strict=True)
    stat = source.stat()
    expected = {
        'status': 'completed',
        'source': str(source),
        'source_size': stat.st_size,
        'source_mtime_ns': stat.st_mtime_ns,
        'output': str(final),
        'profile': 'CB_tagtoname_20260910',
        'min_mapq': min_mapq,
        'min_distance': 10,
        'max_distance': 5000,
        'shift_plus': 4,
        'shift_minus': -5,
        'use_chrom': use_chrom,
        'chunksize': 500000,
        'collapse_within': False,
    }
    for key, value in expected.items():
        if record.get(key) != value:
            raise FragmentReuseError(
                f'manifest.json mismatch for {key}: found {record.get(key)!r}, expected {value!r}'
            )
    for path, key in ((final, 'sha256'), (index, 'index_sha256')):
        if not path.is_file():
            raise FragmentReuseError(f'Missing file: {path}')
        if record.get(key) != file_sha256(path):
            raise FragmentReuseError(f'File checksum mismatch: {path}')
    try:
        with pysam.TabixFile(str(final)) as indexed:
            if list(indexed.contigs) != record.get('contigs'):
                raise FragmentReuseError('Tabix index contigs do not match manifest.json')
    except (OSError, ValueError) as error:
        raise FragmentReuseError(f'Cannot read the tabix index: {error}') from error
    print(f'[INFO] Reusing completed and validated fragments: {final}')
    return str(final)


#------------------------------------------------------------------------------#
def prepare_fragments(args, bam_path, outprefix):
    if not getattr(args, 'make_fragments', False):
        if not args.fragments:
            report_missing_fragments(bam_path)
            raise SystemExit('Fragments are required for single-cell assignment. Supply --fragments or request --make-fragments.')
        return args.fragments
    if args.fragments:
        raise ValueError('Use either --fragments or --make-fragments, not both')
    print('[preprocessing] --make-fragments runs external BAM preprocessing inside this invocation. '
          'For separate runtime and reusable output, run ask-make-fragments first, '
          'then pass --fragments FILE to ASK.', file=sys.stderr, flush=True)
    directory = args.fragments_dir or str(outprefix) + '_sinto_fragments'
    if Path(directory).exists():
        try:
            return reuse_fragments(bam_path, directory, args.fragments_mapq,
                                   args.fragments_use_chrom)
        except FragmentReuseError as error:
            raise FragmentReuseError(f'{directory}: {error}') from error
    try:
        return generate_fragments(bam_path, directory, sinto=args.sinto_bin, samtools=args.samtools_bin,
                                  threads=args.fragments_threads, min_mapq=args.fragments_mapq,
                                  use_chrom=args.fragments_use_chrom)
    except FragmentDependencyError as error:
        raise SystemExit(f"[ERROR] {error}") from None


def main(argv=None):
    """Standalone preprocessing entry point; does not import ASK analysis modules."""
    import argparse
    import shlex
    parser = argparse.ArgumentParser(
        description='Prepare reusable scATAC fragments with sinto, separately from ASK analysis.')
    parser.add_argument('-i', '--bam', required=True, help='Coordinate-sorted BAM with CB tags.')
    parser.add_argument('-o', '--outdir', required=True, help='Output directory; completed matching output is validated and reused.')
    parser.add_argument('--sinto-bin', default='sinto')
    parser.add_argument('--samtools-bin', default='samtools')
    parser.add_argument('-p', '--threads', type=int, default=4)
    parser.add_argument('--min-mapq', type=int, default=20)
    parser.add_argument('--use-chrom', default='(?i)^chr',
                        help='Chromosome regex; use .* to include all contigs.')
    args = parser.parse_args(argv)
    print('[preprocessing] Preparing fragments (generate or validate reuse); ASK analysis is not running.',
          file=sys.stderr, flush=True)
    started = time.monotonic()
    try:
        if Path(args.outdir).exists():
            final = reuse_fragments(args.bam, args.outdir, args.min_mapq, args.use_chrom)
        else:
            final = generate_fragments(args.bam, args.outdir, sinto=args.sinto_bin,
                                       samtools=args.samtools_bin, threads=args.threads,
                                       min_mapq=args.min_mapq, use_chrom=args.use_chrom)
    except FragmentDependencyError as error:
        raise SystemExit(f"[ERROR] {error}") from None
    except (FragmentReuseError, FileExistsError) as error:
        raise SystemExit(f"[ERROR] Cannot generate or reuse fragments: {error}\n"
                         "Check the existing output or choose a directory that does not exist.") from None
    print(f'[preprocessing] Completed in {time.monotonic() - started:.1f} seconds. '
          f'Reuse in ASK with --fragments {shlex.quote(final)}', file=sys.stderr, flush=True)
    print(final)


if __name__ == '__main__':
    main()
