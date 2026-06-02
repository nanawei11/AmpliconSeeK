from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Iterable


def run_command(cmd: list[str], dry_run: bool = False) -> None:
    """Print and run one command."""
    print(" ".join(str(x) for x in cmd))
    if not dry_run:
        subprocess.run([str(x) for x in cmd], check=True)


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def download_sra(
    sample_id: str,
    srr_ids: str | Iterable[str],
    workdir: str | Path = ".",
    threads: int = 12,
    dry_run: bool = False,
) -> tuple[Path, Path | None]:
    """Download SRA runs and return FASTQ paths.

    Requires ``parallel-fastq-dump`` or SRA Toolkit commands in PATH.
    Multiple SRR IDs can be provided as a list or as ``SRR1|SRR2``.
    """
    workdir = Path(workdir)
    fastq_dir = ensure_dir(workdir / sample_id / "fastq")
    if isinstance(srr_ids, str):
        runs = [x for x in srr_ids.replace(",", "|").split("|") if x]
    else:
        runs = list(srr_ids)

    for srr in runs:
        run_command(
            [
                "parallel-fastq-dump",
                "-t",
                str(threads),
                "--split-files",
                "--gzip",
                "-O",
                fastq_dir,
                "-s",
                srr,
            ],
            dry_run=dry_run,
        )

    r1_files = sorted(fastq_dir.glob("*_1.fastq.gz")) + sorted(fastq_dir.glob("*_R1.fastq.gz"))
    r2_files = sorted(fastq_dir.glob("*_2.fastq.gz")) + sorted(fastq_dir.glob("*_R2.fastq.gz"))
    if not r1_files and dry_run:
        return fastq_dir / f"{sample_id}_R1.fastq.gz", fastq_dir / f"{sample_id}_R2.fastq.gz"
    if not r1_files:
        raise FileNotFoundError(f"No R1 FASTQ files found in {fastq_dir}")

    return r1_files[0], r2_files[0] if r2_files else None


def map_fastq_to_bam(
    sample_id: str,
    fastq1: str | Path,
    fastq2: str | Path | None,
    bwa_index: str | Path,
    workdir: str | Path = ".",
    threads: int = 5,
    dry_run: bool = False,
) -> Path:
    """Map FASTQ files with BWA and create a sorted, marked-duplicate BAM."""
    workdir = Path(workdir)
    bwa_dir = ensure_dir(workdir / sample_id / "bwa")
    unsorted_bam = bwa_dir / f"{sample_id}.presort.bam"
    bam = bwa_dir / f"{sample_id}.bam"

    fastqs = [str(fastq1)]
    if fastq2:
        fastqs.append(str(fastq2))

    map_cmd = (
        f"bwa mem -t {threads} {bwa_index} {' '.join(fastqs)} "
        f"| samtools view -Shb - > {unsorted_bam}"
    )
    sort_cmd = (
        f"samtools fixmate --threads {threads} -m {unsorted_bam} - "
        f"| samtools sort --threads {threads} -T {bwa_dir / sample_id} - "
        f"| samtools markdup --threads {threads} -T {bwa_dir / sample_id} -S -s - {bam}"
    )

    print(map_cmd)
    if not dry_run:
        subprocess.run(map_cmd, shell=True, check=True)
    print(sort_cmd)
    if not dry_run:
        subprocess.run(sort_cmd, shell=True, check=True)

    run_command(["samtools", "index", bam], dry_run=dry_run)
    return bam


def run_ask_from_bam(
    bam: str | Path,
    sample_id: str | None = None,
    workdir: str | Path = ".",
    genome: str = "hg38",
    extra_args: Iterable[str] | None = None,
    dry_run: bool = False,
) -> Path:
    """Run ASK on an existing sorted and indexed BAM file."""
    bam = Path(bam)
    if sample_id is None:
        sample_id = bam.name.removesuffix(".bam")
    ask_dir = ensure_dir(Path(workdir) / sample_id / "ask")
    outprefix = ask_dir / sample_id
    cmd = [
        "ask",
        "-i",
        bam,
        "-o",
        outprefix,
        "-g",
        genome,
        "--knn",
        "3",
        "--subseg",
        "--method",
        "both",
        "--juncread",
        "5",
        "--SA_with_nm"
    ]
    if extra_args:
        cmd.extend(extra_args)
    run_command(cmd, dry_run=dry_run)
    return outprefix


def run_from_sra(
    sample_id: str,
    srr_ids: str | Iterable[str],
    bwa_index: str | Path,
    workdir: str | Path = ".",
    genome: str = "hg38",
    threads: int = 5,
    dry_run: bool = False,
) -> Path:
    """Download SRA data, map FASTQ files, and run ASK."""
    fastq1, fastq2 = download_sra(sample_id, srr_ids, workdir, threads=threads, dry_run=dry_run)
    bam = map_fastq_to_bam(sample_id, fastq1, fastq2, bwa_index, workdir, threads=threads, dry_run=dry_run)
    return run_ask_from_bam(bam, sample_id, workdir, genome, dry_run=dry_run)


def run_from_fastq(
    sample_id: str,
    fastq1: str | Path,
    fastq2: str | Path | None,
    bwa_index: str | Path,
    workdir: str | Path = ".",
    genome: str = "hg38",
    threads: int = 5,
    dry_run: bool = False,
) -> Path:
    """Map FASTQ files and run ASK."""
    bam = map_fastq_to_bam(sample_id, fastq1, fastq2, bwa_index, workdir, threads=threads, dry_run=dry_run)
    return run_ask_from_bam(bam, sample_id, workdir, genome, dry_run=dry_run)


def run_from_bam(
    sample_id: str,
    bam: str | Path,
    workdir: str | Path = ".",
    genome: str = "hg38",
    dry_run: bool = False,
) -> Path:
    """Run ASK from a prepared BAM file."""
    return run_ask_from_bam(bam, sample_id, workdir, genome, dry_run=dry_run)
