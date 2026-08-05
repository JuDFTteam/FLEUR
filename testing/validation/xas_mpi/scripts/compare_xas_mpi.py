#!/usr/bin/env python3
"""Compare serial, pure-k-point, and shared-k-point XAS MPI runs."""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import Counter
from pathlib import Path


RUNS = {
    "np1": ("runs/np1", "xas_mpi", 1, 1),
    "np2_pure": ("runs/np2_pure", "xas_mpi", 2, 1),
    "np4_pure": ("runs/np4_pure", "xas_mpi", 4, 1),
    "np4_shared": ("runs/np4_shared", "xas_mpi", 4, 2),
}
POLS = ("x", "y", "z")
EXPECTED_TRANSITIONS = 2810
EXPECTED_FULL_KPTS = set(range(1, 28))
EXPECTED_PARENT_KPTS = set(range(1, 7))
ABS_TOL = 1.0e-12
REL_TOL = 1.0e-10
GRID_TOL = 1.0e-12
ARITH_ABS_TOL = 1.0e-18
ARITH_REL_TOL = 1.0e-12


class ValidationError(RuntimeError):
    pass


def close_enough(abs_diff: float, scale: float, abs_tol: float = ABS_TOL, rel_tol: float = REL_TOL) -> bool:
    rel_diff = abs_diff / scale if scale > 0.0 else (0.0 if abs_diff == 0.0 else math.inf)
    return abs_diff < abs_tol or rel_diff < rel_tol


def read_spectrum(path: Path) -> tuple[list[float], list[float]]:
    grid: list[float] = []
    intensity: list[float] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 2:
                raise ValidationError(f"Malformed spectrum row in {path}: {line.rstrip()}")
            grid.append(float(fields[0]))
            intensity.append(float(fields[1]))
    if len(grid) < 2:
        raise ValidationError(f"Spectrum has fewer than two rows: {path}")
    return grid, intensity


def spectrum_stats(grid: list[float], values: list[float]) -> dict[str, float | int]:
    spacing = grid[1] - grid[0]
    trap = sum(0.5 * (values[i] + values[i + 1]) * (grid[i + 1] - grid[i]) for i in range(len(grid) - 1))
    peak_index = max(range(len(values)), key=values.__getitem__)
    return {
        "rows": len(grid),
        "emin": grid[0],
        "emax": grid[-1],
        "spacing": spacing,
        "trap": trap,
        "rect": sum(values) * spacing,
        "min": min(values),
        "max": values[peak_index],
        "peak_energy": grid[peak_index],
        "negative": sum(value < -1.0e-14 for value in values),
        "nonfinite": sum(not math.isfinite(value) for value in grid + values),
    }


def rank_from_path(path: Path) -> int:
    match = re.search(r"_rank([0-9]{4})\.dat$", path.name)
    if not match:
        raise ValidationError(f"Cannot parse MPI rank from {path}")
    return int(match.group(1))


def read_transition_file(path: Path) -> list[tuple[tuple[int, ...], tuple[float, ...]]]:
    rows = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 13:
                raise ValidationError(f"Expected 13 transition columns in {path}: {line.rstrip()}")
            key = tuple(int(value) for value in fields[:6])
            values = tuple(float(value) for value in fields[6:])
            rows.append((key, values))
    return rows


def transition_stats(paths: list[Path]) -> dict[str, object]:
    by_rank: dict[int, list[tuple[tuple[int, ...], tuple[float, ...]]]] = {}
    all_rows = []
    for path in paths:
        rank = rank_from_path(path)
        if rank in by_rank:
            raise ValidationError(f"Duplicate transition file for rank {rank}")
        by_rank[rank] = read_transition_file(path)
        all_rows.extend(by_rank[rank])

    keys = [row[0] for row in all_rows]
    duplicate_keys = [key for key, count in Counter(keys).items() if count > 1]
    nonfinite = 0
    negative_strength = 0
    max_arith_abs = 0.0
    max_arith_rel = 0.0
    weighted_sum = 0.0
    for _, values in all_rows:
        transition_energy, transition_ev, occupation, k_weight, one_minus_occ, abs_m2, weighted = values
        nonfinite += sum(not math.isfinite(value) for value in values)
        negative_strength += int(abs_m2 < -1.0e-14 or weighted < -1.0e-14)
        expected = k_weight * one_minus_occ * abs_m2
        error = abs(weighted - expected)
        scale = max(abs(weighted), abs(expected))
        max_arith_abs = max(max_arith_abs, error)
        max_arith_rel = max(max_arith_rel, error / scale if scale else 0.0)
        weighted_sum += weighted

    return {
        "by_rank": by_rank,
        "rows": len(all_rows),
        "keys": set(keys),
        "duplicates": duplicate_keys,
        "full_kpts": {key[0] for key in keys},
        "parent_kpts": {key[1] for key in keys},
        "weighted_sum": weighted_sum,
        "nonfinite": nonfinite,
        "negative_strength": negative_strength,
        "max_arith_abs": max_arith_abs,
        "max_arith_rel": max_arith_rel,
    }


def read_log_metadata(run_dir: Path) -> dict[str, object]:
    stdout_path = run_dir / "fleur_stdout.log"
    stderr_path = run_dir / "fleur_stderr.log"
    if not stdout_path.exists() or not stderr_path.exists():
        raise ValidationError(f"Missing FLEUR stdout/stderr in {run_dir}")
    stdout = stdout_path.read_text(encoding="utf-8", errors="replace")
    stderr = stderr_path.read_text(encoding="utf-8", errors="replace")
    environment_files = sorted(run_dir.glob("job_environment_*.log"))
    if len(environment_files) != 1:
        raise ValidationError(f"Expected one job-environment log in {run_dir}, found {len(environment_files)}")
    environment = environment_files[0].read_text(encoding="utf-8", errors="replace")
    match = re.search(r"Requested PEs per k-point:\s*([0-9]+)", environment, re.IGNORECASE)
    commit_match = re.search(r"Source commit:\s*([0-9a-f]{40})", environment)
    if not match or not commit_match:
        raise ValidationError(f"Missing requested layout or source commit in {environment_files[0]}")
    pe_node_match = re.search(r"Number of PE/node\s*:\s*([0-9]+)", stdout, re.IGNORECASE)
    debug_files = sorted(run_dir.glob("xas_debug_*.txt"))
    if len(debug_files) != 1:
        raise ValidationError(f"Expected one XAS debug file in {run_dir}, found {len(debug_files)}")
    debug = debug_files[0].read_text(encoding="utf-8", errors="replace")
    weight_match = re.search(
        r"star weight sums: selected=\s*([^\s]+)\s+expanded=\s*([^\s]+)\s+diff=\s*([^\s]+)", debug
    )
    if not weight_match:
        raise ValidationError(f"Missing parent/star weight summary in {debug_files[0]}")
    strength_match = re.search(
        r"strength summary total.*?\s([+-]?[0-9.]+E[+-][0-9]+)\s+"
        r"([+-]?[0-9.]+E[+-][0-9]+)\s+([+-]?[0-9.]+E[+-][0-9]+)",
        debug,
    )
    if not strength_match:
        raise ValidationError(f"Missing additive strength summary in {debug_files[0]}")
    bad_runtime = re.search(r"\b(?:abort|segmentation fault|nan|[+-]?inf(?:inity)?)\b", stdout + "\n" + stderr, re.IGNORECASE)
    underflow_match = re.search(r"underflow occurred.*?([0-9]+) time", stdout, re.IGNORECASE)
    return {
        "pe_per_kpt": int(match.group(1)) if match else None,
        "pe_per_node": int(pe_node_match.group(1)) if pe_node_match else None,
        "weights": tuple(float(weight_match.group(i)) for i in range(1, 4)),
        "strengths": tuple(float(strength_match.group(i)) for i in range(1, 4)),
        "underflow": int(underflow_match.group(1)) if underflow_match else 0,
        "completion": (run_dir / "RUN_COMPLETE").exists() and "XAS wrote spectrum to" in stdout,
        "bad_runtime_text": bad_runtime.group(0) if bad_runtime else None,
        "source_commit": commit_match.group(1),
        "environment_path": environment_files[0],
    }


def load_run(root: Path, label: str, spec: tuple[str, str, int, int]) -> dict[str, object]:
    directory, prefix, ranks, pe_per_kpt = spec
    run_dir = root / directory
    if not run_dir.is_dir():
        raise ValidationError(f"Missing run directory: {run_dir}")
    data: dict[str, object] = {
        "dir": run_dir,
        "prefix": prefix,
        "ranks": ranks,
        "pe_per_kpt": pe_per_kpt,
    }
    data["log"] = read_log_metadata(run_dir)
    spectra = {}
    transitions = {}
    for pol in POLS:
        spectrum_path = run_dir / f"{prefix}_L3_{pol}_eta0p030.dat"
        if not spectrum_path.exists():
            raise ValidationError(f"Missing spectrum: {spectrum_path}")
        grid, intensity = read_spectrum(spectrum_path)
        spectra[pol] = {"path": spectrum_path, "grid": grid, "intensity": intensity, "stats": spectrum_stats(grid, intensity)}

        paths = sorted(run_dir.glob(f"{prefix}_L3_{pol}_transitions_rank*.dat"))
        if len(paths) != ranks:
            raise ValidationError(f"{label}/{pol}: expected {ranks} transition files, found {len(paths)}")
        stats = transition_stats(paths)
        if set(stats["by_rank"]) != set(range(ranks)):
            raise ValidationError(f"{label}/{pol}: rank-file set is {sorted(stats['by_rank'])}")
        transitions[pol] = stats
    data["spectra"] = spectra
    data["transitions"] = transitions
    return data


def compare(root: Path, skip_reference: bool = False) -> bool:
    ok = True
    runs = {label: load_run(root, label, spec) for label, spec in RUNS.items()}
    baseline = runs["np1"]
    expected_commit = (root / "provenance" / "PHASE1_COMMIT").read_text(encoding="utf-8").strip()

    for label, run in runs.items():
        log = run["log"]
        pe_per_kpt = log["pe_per_kpt"]
        expected_pe = run["pe_per_kpt"]
        layout_ok = pe_per_kpt == expected_pe
        weights_ok = (
            abs(log["weights"][0] - 1.0) < ABS_TOL
            and abs(log["weights"][1] - 1.0) < ABS_TOL
            and abs(log["weights"][2]) < ABS_TOL
        )
        run_ok = (
            bool(log["completion"])
            and not log["bad_runtime_text"]
            and layout_ok
            and weights_ok
            and log["source_commit"] == expected_commit
        )
        ok &= run_ok
        print(
            f"{label}: ranks={run['ranks']} pe_per_node={log['pe_per_node']} "
            f"pe_per_kpt={pe_per_kpt} completion={log['completion']} "
            f"runtime_issue={log['bad_runtime_text']} layout_status={'PASS' if layout_ok else 'FAIL'}"
        )
        print(
            f"  source_commit={log['source_commit']} parent/star weights: "
            f"selected={log['weights'][0]:.16e} expanded={log['weights'][1]:.16e} "
            f"diff={log['weights'][2]:.3e}; strengths={log['strengths']}; underflow={log['underflow']}"
        )

        expected_nonempty = run["ranks"] // expected_pe
        for pol in POLS:
            spectrum = run["spectra"][pol]
            stats = spectrum["stats"]
            transitions = run["transitions"][pol]
            row_counts = {rank: len(rows) for rank, rows in transitions["by_rank"].items()}
            rank_kpts = {
                rank: sorted({row[0][0] for row in rows})
                for rank, rows in transitions["by_rank"].items()
            }
            nonempty = sum(count > 0 for count in row_counts.values())
            coverage_ok = (
                transitions["full_kpts"] == EXPECTED_FULL_KPTS
                and transitions["parent_kpts"] == EXPECTED_PARENT_KPTS
            )
            arithmetic_ok = (
                transitions["max_arith_abs"] < ARITH_ABS_TOL
                or transitions["max_arith_rel"] < ARITH_REL_TOL
            )
            local_ok = (
                stats["negative"] == 0
                and stats["nonfinite"] == 0
                and transitions["nonfinite"] == 0
                and transitions["negative_strength"] == 0
                and not transitions["duplicates"]
                and transitions["rows"] == EXPECTED_TRANSITIONS
                and coverage_ok
                and nonempty == expected_nonempty
                and arithmetic_ok
            )
            ok &= local_ok
            print(
                f"  {pol}: rows={stats['rows']} range=[{stats['emin']:.12e},{stats['emax']:.12e}] "
                f"step={stats['spacing']:.12e} trap={stats['trap']:.16e} rect={stats['rect']:.16e} "
                f"max={stats['max']:.16e}@{stats['peak_energy']:.12e} min={stats['min']:.3e}"
            )
            print(
                f"     transitions={transitions['rows']} per_rank={row_counts} "
                f"rank_full_k={rank_kpts} "
                f"weighted_sum={transitions['weighted_sum']:.16e} duplicates={len(transitions['duplicates'])} "
                f"full_k={sorted(transitions['full_kpts'])} parents={sorted(transitions['parent_kpts'])} "
                f"arith_abs={transitions['max_arith_abs']:.3e} arith_rel={transitions['max_arith_rel']:.3e} "
                f"status={'PASS' if local_ok else 'FAIL'}"
            )

    print("Cross-layout comparison against np1:")
    for label, run in runs.items():
        if label == "np1":
            continue
        for pol in POLS:
            ref_spec = baseline["spectra"][pol]
            cur_spec = run["spectra"][pol]
            if len(ref_spec["grid"]) != len(cur_spec["grid"]):
                raise ValidationError(f"{label}/{pol}: spectrum row count differs from np1")
            grid_diff = max(abs(a - b) for a, b in zip(ref_spec["grid"], cur_spec["grid"]))
            spec_diff = max(abs(a - b) for a, b in zip(ref_spec["intensity"], cur_spec["intensity"]))
            common_max = max(max(map(abs, ref_spec["intensity"])), max(map(abs, cur_spec["intensity"])))
            rel_diff = spec_diff / common_max if common_max else 0.0
            ref_trans = baseline["transitions"][pol]
            cur_trans = run["transitions"][pol]
            row_keys_ok = ref_trans["keys"] == cur_trans["keys"]
            row_count_ok = ref_trans["rows"] == cur_trans["rows"]
            weighted_diff = abs(ref_trans["weighted_sum"] - cur_trans["weighted_sum"])
            weighted_scale = max(abs(ref_trans["weighted_sum"]), abs(cur_trans["weighted_sum"]))
            strength_diff = max(
                abs(a - b) for a, b in zip(baseline["log"]["strengths"], run["log"]["strengths"])
            )
            strength_scale = max(
                max(map(abs, baseline["log"]["strengths"])), max(map(abs, run["log"]["strengths"]))
            )
            weight_diff = max(abs(a - b) for a, b in zip(baseline["log"]["weights"], run["log"]["weights"]))
            underflow_equal = baseline["log"]["underflow"] == run["log"]["underflow"]
            comparison_ok = (
                grid_diff < GRID_TOL
                and close_enough(spec_diff, common_max)
                and row_keys_ok
                and row_count_ok
                and close_enough(weighted_diff, weighted_scale)
                and close_enough(strength_diff, strength_scale)
                and weight_diff < ABS_TOL
                and underflow_equal
            )
            ok &= comparison_ok
            print(
                f"  {label}/{pol}: grid_max_abs={grid_diff:.3e} spectrum_max_abs={spec_diff:.3e} "
                f"spectrum_rel={rel_diff:.3e} rows_equal={row_count_ok} keys_equal={row_keys_ok} "
                f"weighted_sum_abs={weighted_diff:.3e} strength_abs={strength_diff:.3e} "
                f"weight_abs={weight_diff:.3e} underflow_equal={underflow_equal} "
                f"status={'PASS' if comparison_ok else 'FAIL'}"
            )

    if skip_reference:
        print("Canonical serial-reference comparison: SKIPPED by explicit request")
    else:
        print("Canonical serial-reference comparison:")
        for pol in POLS:
            reference_path = root / "reference" / "canonical" / f"xas_mpi_L3_{pol}_eta0p030.dat"
            if not reference_path.exists():
                raise ValidationError(
                    "Canonical JURECA reference is absent; use --skip-reference only for the pre-promotion checkpoint"
                )
            ref_grid, ref_values = read_spectrum(reference_path)
            current = baseline["spectra"][pol]
            if len(ref_grid) != len(current["grid"]):
                raise ValidationError(f"Canonical reference row count differs for polarization {pol}")
            grid_diff = max(abs(a - b) for a, b in zip(ref_grid, current["grid"]))
            spec_diff = max(abs(a - b) for a, b in zip(ref_values, current["intensity"]))
            scale = max(max(map(abs, ref_values)), max(map(abs, current["intensity"])))
            ref_ok = grid_diff < GRID_TOL and close_enough(spec_diff, scale)
            ok &= ref_ok
            print(
                f"  {pol}: grid_max_abs={grid_diff:.3e} spectrum_max_abs={spec_diff:.3e} "
                f"spectrum_rel={(spec_diff / scale if scale else 0.0):.3e} "
                f"status={'PASS' if ref_ok else 'FAIL'}"
            )

    print(f"OVERALL STATUS: {'PASS' if ok else 'FAIL'}")
    return ok


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--skip-reference",
        action="store_true",
        help="validate layouts before the clean np1 result is promoted to the canonical reference",
    )
    args = parser.parse_args()
    try:
        return 0 if compare(args.root.resolve(), skip_reference=args.skip_reference) else 1
    except (OSError, ValueError, ValidationError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
