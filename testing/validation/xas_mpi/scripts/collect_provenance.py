#!/usr/bin/env python3
"""Generate a machine-readable manifest for completed XAS MPI validation runs."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path


RUNS = ("np1", "np2_pure", "np4_pure", "np4_shared")
EXPECTED_COMMIT = "b3818dad42aa56fa3a835a9cf8e972e23ea2dd87"
EXPECTED_EXECUTABLE_SHA256 = "df4af8292c337614a5755910c7e59634fe638c7cf1ff076127dd501d577e65c4"
EXPECTED_LAYOUTS = {
    "np1": (1, 1),
    "np2_pure": (2, 1),
    "np4_pure": (4, 1),
    "np4_shared": (4, 2),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_sha256s(path: Path) -> dict[str, str]:
    checksums = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        checksum, filename = line.split(maxsplit=1)
        checksums[filename.strip()] = checksum
    return checksums


def verify_sha256s(directory: Path, manifest: Path) -> dict[str, str]:
    checksums = read_sha256s(manifest)
    for filename, expected in checksums.items():
        actual = sha256(directory / filename)
        if actual != expected:
            raise RuntimeError(f"checksum mismatch for {directory / filename}: {actual} != {expected}")
    return checksums


def log_section(text: str, start: str, end: str) -> list[str]:
    if start not in text or end not in text:
        raise RuntimeError(f"missing provenance section {start!r} ... {end!r}")
    return [line for line in text.split(start, 1)[1].split(end, 1)[0].strip().splitlines() if line]


def required_files(run_dir: Path) -> list[Path]:
    paths = [
        run_dir / "inp.xml",
        run_dir / "kpts.xml",
        run_dir / "sym.xml",
        run_dir / "cdn.hdf",
        run_dir / "fleur_stdout.log",
        run_dir / "fleur_stderr.log",
        run_dir / "out.xml",
        run_dir / "result_sha256.txt",
        run_dir / "RUN_COMPLETE",
    ]
    paths.extend(sorted(run_dir.glob("xas_mpi_L3_*_eta0p030.dat")))
    paths.extend(sorted(run_dir.glob("xas_mpi_L3_*_transitions_rank*.dat")))
    missing = [path for path in paths if not path.exists()]
    if missing:
        raise RuntimeError(f"missing required files in {run_dir}: {missing}")
    return paths


def xml_metadata(path: Path) -> dict[str, object]:
    root = ET.parse(path).getroot()
    program = root.find("programVersion")
    compilation = program.find("compilationInfo")
    git_info = program.find("gitInfo")
    execution = program.find("programExecution")
    parallel = root.find("parallelSetup")
    mpi = parallel.find("mpi")
    omp = parallel.find("openMP")
    start = root.find("startDateAndTime")
    end = root.find("endDateAndTime")
    return {
        "program_version": program.attrib,
        "compilation": compilation.attrib,
        "git": git_info.attrib,
        "execution": execution.attrib,
        "mpi": mpi.attrib if mpi is not None else None,
        "openmp": omp.attrib if omp is not None else None,
        "start": start.attrib,
        "end": end.attrib,
    }


def run_manifest(run_dir: Path, label: str) -> dict[str, object]:
    environment_logs = sorted(run_dir.glob("job_environment_*.log"))
    if len(environment_logs) != 1:
        raise RuntimeError(f"expected one environment log in {run_dir}, found {len(environment_logs)}")
    environment = environment_logs[0]
    text = environment.read_text(encoding="utf-8", errors="replace")
    job_match = re.search(r"SLURM_JOB_ID=([0-9]+)", text)
    commit_match = re.search(r"Source commit:\s*([0-9a-f]{40})", text)
    executable_match = re.search(r"FLEUR executable:\s*(.+)", text)
    executable_hash_match = re.search(r"FLEUR executable SHA256:\s*([0-9a-f]{64})", text)
    executable_build_match = re.search(r"^Modify:\s*([^\r\n]+)", text, re.MULTILINE)
    launcher_match = re.search(r"Exact launcher command:(.+)", text)
    pe_match = re.search(r"Requested PEs per k-point:\s*([0-9]+)", text)
    ranks_match = re.search(r"Expected MPI ranks:\s*([0-9]+)", text)
    clean_match = re.search(r"Source status --short:\s*(\S+)", text)
    host_match = re.search(r"Host:\s*(\S+)", text)
    start_match = re.search(r"Started at\s*(\S+)", text)
    finish_match = re.search(r"Finished at\s*(\S+)", text)
    exit_match = re.search(r"Exit status:\s*([0-9]+)", text)
    mklroot_match = re.search(r"MKLROOT=(\S+)", text)
    hdf5_match = re.search(r"HDF5 Version:\s*(\S+)", text)
    cmake_match = re.search(r"CMake version:\s*\ncmake version\s+(\S+)", text)
    libxml_match = re.search(r"libxml2 version:\s*\n(\S+)", text)
    post_restart_match = re.search(r"Restart checksum after execution:\s*([0-9a-f]{64})", text)
    required_matches = (
        job_match,
        commit_match,
        executable_match,
        executable_hash_match,
        executable_build_match,
        launcher_match,
        pe_match,
        ranks_match,
        clean_match,
        host_match,
        start_match,
        finish_match,
        exit_match,
        mklroot_match,
        hdf5_match,
        cmake_match,
        libxml_match,
        post_restart_match,
    )
    if not all(required_matches):
        raise RuntimeError(f"incomplete provenance log: {environment}")
    pre_run_block = log_section(text, "Input and original restart checksums before execution:\n", "Expected MPI ranks:")
    pre_run_sha256s = {filename: checksum for checksum, filename in (line.split(maxsplit=1) for line in pre_run_block)}
    seed_sha256s = read_sha256s(run_dir / "SEED_SHA256")
    if pre_run_sha256s != seed_sha256s:
        raise RuntimeError(f"pre-run checksums and SEED_SHA256 differ in {run_dir}")
    result_sha256s = verify_sha256s(run_dir, run_dir / "result_sha256.txt")
    mpi_ranks = int(ranks_match.group(1))
    pe_per_kpt = int(pe_match.group(1))
    if (mpi_ranks, pe_per_kpt) != EXPECTED_LAYOUTS[label]:
        raise RuntimeError(f"unexpected MPI layout for {label}: {(mpi_ranks, pe_per_kpt)}")
    if commit_match.group(1) != EXPECTED_COMMIT or clean_match.group(1) != "CLEAN":
        raise RuntimeError(f"source provenance is not the approved clean commit in {environment}")
    if executable_hash_match.group(1) != EXPECTED_EXECUTABLE_SHA256:
        raise RuntimeError(f"unexpected executable checksum in {environment}")
    if int(exit_match.group(1)) != 0 or "Run completed successfully." not in text:
        raise RuntimeError(f"run did not complete successfully: {environment}")
    files = required_files(run_dir)
    files.append(environment)
    return {
        "slurm_job_id": job_match.group(1) if job_match else None,
        "host": host_match.group(1),
        "started": start_match.group(1),
        "finished": finish_match.group(1),
        "exit_status": int(exit_match.group(1)),
        "source": {
            "commit": commit_match.group(1),
            "status_short": clean_match.group(1),
        },
        "executable": {
            "path": executable_match.group(1).strip(),
            "sha256": executable_hash_match.group(1),
            "filesystem_modify_time": executable_build_match.group(1).strip(),
        },
        "launcher": launcher_match.group(1).strip(),
        "mpi_layout": {
            "ranks": mpi_ranks,
            "openmp_threads": 1,
            "pe_per_kpt": pe_per_kpt,
        },
        "module_environment": log_section(text, "Loaded modules:\n", "Compiler versions:\n"),
        "compiler_versions": log_section(text, "Compiler versions:\n", "MPI wrapper configuration:\n"),
        "mpi_wrapper_configuration": log_section(text, "MPI wrapper configuration:\n", "MPI and launcher versions:\n"),
        "mpi_and_launcher_versions": log_section(text, "MPI and launcher versions:\n", "MKL environment:\n"),
        "libraries": {
            "mklroot": mklroot_match.group(1),
            "hdf5_version": hdf5_match.group(1),
            "cmake_version": cmake_match.group(1),
            "libxml2_version": libxml_match.group(1),
            "dynamic_linkage": log_section(text, "Relevant dynamic libraries:\n", "Exact launcher command:"),
        },
        "seed_sha256s": seed_sha256s,
        "post_run_restart_sha256": post_restart_match.group(1),
        "result_sha256s": result_sha256s,
        "environment_log": environment.name,
        "environment_log_sha256": sha256(environment),
        "out_xml": xml_metadata(run_dir / "out.xml"),
        "files": {path.name: sha256(path) for path in files},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--output",
        type=Path,
        help="default: <root>/provenance/clean_jureca_manifest.json",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output or root / "provenance" / "clean_jureca_manifest.json"
    if output.exists():
        print(f"Refusing to overwrite existing manifest: {output}", file=sys.stderr)
        return 2
    try:
        fixture_sha256s = verify_sha256s(root / "fixtures", root / "fixtures" / "SHA256SUMS")
        canonical_sha256s = verify_sha256s(
            root / "reference" / "canonical", root / "reference" / "canonical" / "SHA256SUMS"
        )
        runs = {label: run_manifest(root / "runs" / label, label) for label in RUNS}
        for label, run in runs.items():
            if run["seed_sha256s"] != fixture_sha256s:
                raise RuntimeError(f"{label} did not use the tracked fixture checksums")
        for polarization in ("x", "y", "z"):
            filename = f"xas_mpi_L3_{polarization}_eta0p030.dat"
            if canonical_sha256s[filename] != runs["np1"]["files"][filename]:
                raise RuntimeError(f"canonical {polarization} spectrum is not the approved np1 result")
        manifest = {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "phase1_commit": (root / "provenance" / "PHASE1_COMMIT").read_text(encoding="utf-8").strip(),
            "approved_executable_sha256": EXPECTED_EXECUTABLE_SHA256,
            "fixture_sha256s": fixture_sha256s,
            "canonical_reference_sha256s": canonical_sha256s,
            "runs": runs,
        }
        output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    except (OSError, RuntimeError, ET.ParseError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
