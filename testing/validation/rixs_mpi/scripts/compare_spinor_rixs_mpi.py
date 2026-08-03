#!/usr/bin/env python3
"""Compare 1-, 2-, and 4-rank spinor RIXS contribution validations."""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


SPECTRUM_ABS_TOL = 1.0e-12
SPECTRUM_REL_TOL = 1.0e-10
ROW_ABS_TOL = 1.0e-18
ROW_REL_TOL = 1.0e-12
CONTRIBUTION_ENERGY_ABS_TOL = 1.0e-12
MANIFOLD_DEGENERACY_TOL = 1.0e-10
GRID_TOL = 1.0e-13
RUNS = (("np1", "SQAy_np1", 1), ("np2", "SQAy_np2", 2), ("np4", "SQAy_np4", 4))


@dataclass(frozen=True)
class Spectrum:
    path: Path
    loss: tuple[float, ...]
    intensity: tuple[float, ...]


@dataclass(frozen=True)
class ContributionRow:
    ikpt: int
    band_v: int
    band_n: int
    absorber_atom: int
    absorber_type: int
    eps_v: float
    eps_n: float
    core_energy: float
    occupation_v: float
    occupation_n: float
    k_weight: float
    loss_energy: float
    loss_energy_ev: float
    denominator_real: float
    denominator_imag: float
    denominator_abs2: float
    amplitude_real: float
    amplitude_imag: float
    amplitude_abs2: float
    weighted_strength: float

    @property
    def transition_key(self) -> tuple[int, int, int, int, int]:
        # Contribution files are polarization-specific, so these are all of the
        # discrete identifiers needed to distinguish transitions within a file.
        return (
            self.ikpt,
            self.band_v,
            self.band_n,
            self.absorber_atom,
            self.absorber_type,
        )

    @property
    def numeric_values(self) -> tuple[float, ...]:
        return (
            self.eps_v,
            self.eps_n,
            self.core_energy,
            self.occupation_v,
            self.occupation_n,
            self.k_weight,
            self.loss_energy,
            self.loss_energy_ev,
            self.denominator_real,
            self.denominator_imag,
            self.denominator_abs2,
            self.amplitude_real,
            self.amplitude_imag,
            self.amplitude_abs2,
            self.weighted_strength,
        )


@dataclass
class ContributionAnalysis:
    files: list[Path]
    rows_by_rank: dict[int, list[ContributionRow]]
    rows: list[ContributionRow]
    ikpts_by_rank: dict[int, set[int]]
    combined_ikpts: set[int]
    missing_ikpts: set[int]
    unexpected_ikpts: set[int]
    duplicate_keys: list[tuple[int, int, int, int, int]]
    nonfinite_rows: int
    amplitude_identity_abs: float
    amplitude_identity_rel: float
    weighted_arithmetic_abs: float
    weighted_arithmetic_rel: float
    minimum_amplitude_abs2: float
    minimum_weighted_strength: float
    sum_weighted_strength: float


@dataclass
class LogAnalysis:
    files: list[Path]
    has_spinor_summary: bool
    has_completion: bool
    contribution_status: str
    reconstruction_abs: float | None
    reconstruction_rel: float | None
    abort_hits: list[str]
    nonfinite_hits: list[str]
    distribution_warning_hits: list[str]


@dataclass
class RunAnalysis:
    label: str
    directory: Path
    expected_ranks: int
    spectrum: Spectrum
    contribution: ContributionAnalysis
    logs: LogAnalysis


def numeric_lines(path: Path) -> Iterable[list[str]]:
    with path.open("r", encoding="utf-8", errors="strict") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if not fields:
                continue
            yield fields


def locate_one(directory: Path, pattern: str, description: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise ValueError(
            f"{directory}: expected one {description} matching {pattern!r}, "
            f"found {len(matches)}"
        )
    return matches[0]


def read_spectrum(path: Path) -> Spectrum:
    loss: list[float] = []
    intensity: list[float] = []
    for fields in numeric_lines(path):
        if len(fields) < 3:
            raise ValueError(f"{path}: spectrum row has fewer than three columns")
        loss.append(float(fields[0]))
        intensity.append(float(fields[2]))
    if len(loss) < 2:
        raise ValueError(f"{path}: spectrum needs at least two data rows")
    return Spectrum(path=path, loss=tuple(loss), intensity=tuple(intensity))


def spectrum_metrics(spectrum: Spectrum) -> dict[str, float | int]:
    loss = spectrum.loss
    intensity = spectrum.intensity
    spacings = [loss[index + 1] - loss[index] for index in range(len(loss) - 1)]
    spacing = (loss[-1] - loss[0]) / (len(loss) - 1)
    trapezoid = math.fsum(
        spacings[index] * (intensity[index + 1] + intensity[index]) / 2.0
        for index in range(len(spacings))
    )
    rectangle = math.fsum(intensity) * spacing
    finite_values = [value for value in intensity if math.isfinite(value)]
    if not finite_values:
        minimum = maximum = peak_loss = math.nan
    else:
        minimum = min(finite_values)
        maximum = max(finite_values)
        peak_loss = loss[intensity.index(maximum)]
    return {
        "rows": len(loss),
        "loss_min": loss[0],
        "loss_max": loss[-1],
        "spacing": spacing,
        "max_spacing_deviation": max(abs(value - spacing) for value in spacings),
        "trapezoidal_integral": trapezoid,
        "rectangular_sum": rectangle,
        "endpoint_correction": spacing * (intensity[0] + intensity[-1]) / 2.0,
        "minimum_intensity": minimum,
        "maximum_intensity": maximum,
        "peak_loss": peak_loss,
        "below_negative_tolerance": sum(value < -1.0e-14 for value in finite_values),
        "nonfinite_values": sum(
            not math.isfinite(value) for value in (*loss, *intensity)
        ),
    }


def max_abs_and_relative(
    left: Sequence[float], right: Sequence[float]
) -> tuple[float, float]:
    if len(left) != len(right):
        return math.inf, math.inf
    differences = [abs(a - b) for a, b in zip(left, right)]
    if any(not math.isfinite(value) for value in differences):
        return math.inf, math.inf
    maximum_absolute = max(differences, default=0.0)
    common_maximum = max(
        (abs(value) for value in (*left, *right) if math.isfinite(value)),
        default=0.0,
    )
    relative = maximum_absolute / max(common_maximum, sys.float_info.min)
    return maximum_absolute, relative


def passes_tolerance(maximum_absolute: float, maximum_relative: float) -> bool:
    return (
        maximum_absolute < SPECTRUM_ABS_TOL
        or maximum_relative < SPECTRUM_REL_TOL
    )


def parse_contribution_row(fields: list[str], path: Path) -> ContributionRow:
    if len(fields) != 20:
        raise ValueError(f"{path}: expected 20 contribution columns, found {len(fields)}")
    integers = [int(value) for value in fields[:5]]
    values = [float(value) for value in fields[5:]]
    return ContributionRow(*integers, *values)


def rank_from_filename(path: Path) -> int:
    match = re.search(r"_contrib_rank(\d+)\.dat$", path.name)
    if not match:
        raise ValueError(f"Cannot extract MPI rank from {path.name}")
    return int(match.group(1))


def update_error(
    actual: float,
    expected: float,
    maximum_absolute: float,
    maximum_relative: float,
) -> tuple[float, float]:
    difference = abs(actual - expected)
    scale = max(abs(actual), abs(expected), sys.float_info.min)
    return max(maximum_absolute, difference), max(
        maximum_relative, difference / scale
    )


def analyse_contributions(
    directory: Path, expected_ranks: int, expected_ikpts: set[int]
) -> ContributionAnalysis:
    files = sorted(directory.glob("*_contrib_rank*.dat"))
    rows_by_rank: dict[int, list[ContributionRow]] = {}
    amplitude_identity_abs = 0.0
    amplitude_identity_rel = 0.0
    weighted_arithmetic_abs = 0.0
    weighted_arithmetic_rel = 0.0
    nonfinite_rows = 0

    for path in files:
        rank = rank_from_filename(path)
        if rank in rows_by_rank:
            raise ValueError(f"{directory}: multiple contribution files for rank {rank}")
        rank_rows = [parse_contribution_row(fields, path) for fields in numeric_lines(path)]
        rows_by_rank[rank] = rank_rows
        for row in rank_rows:
            if any(not math.isfinite(value) for value in row.numeric_values):
                nonfinite_rows += 1
                continue
            amplitude_expected = (
                row.amplitude_real * row.amplitude_real
                + row.amplitude_imag * row.amplitude_imag
            )
            amplitude_identity_abs, amplitude_identity_rel = update_error(
                row.amplitude_abs2,
                amplitude_expected,
                amplitude_identity_abs,
                amplitude_identity_rel,
            )
            weighted_expected = (
                row.k_weight
                * row.occupation_v
                * (1.0 - row.occupation_n)
                * row.amplitude_abs2
            )
            weighted_arithmetic_abs, weighted_arithmetic_rel = update_error(
                row.weighted_strength,
                weighted_expected,
                weighted_arithmetic_abs,
                weighted_arithmetic_rel,
            )

    rows = [row for rank in sorted(rows_by_rank) for row in rows_by_rank[rank]]
    key_counts = Counter(row.transition_key for row in rows)
    duplicate_keys = sorted(key for key, count in key_counts.items() if count > 1)
    ikpts_by_rank = {
        rank: {row.ikpt for row in rank_rows}
        for rank, rank_rows in sorted(rows_by_rank.items())
    }
    combined_ikpts = set().union(*ikpts_by_rank.values()) if ikpts_by_rank else set()
    finite_rows = [
        row for row in rows if all(math.isfinite(value) for value in row.numeric_values)
    ]
    return ContributionAnalysis(
        files=files,
        rows_by_rank=rows_by_rank,
        rows=rows,
        ikpts_by_rank=ikpts_by_rank,
        combined_ikpts=combined_ikpts,
        missing_ikpts=expected_ikpts - combined_ikpts,
        unexpected_ikpts=combined_ikpts - expected_ikpts,
        duplicate_keys=duplicate_keys,
        nonfinite_rows=nonfinite_rows,
        amplitude_identity_abs=amplitude_identity_abs,
        amplitude_identity_rel=amplitude_identity_rel,
        weighted_arithmetic_abs=weighted_arithmetic_abs,
        weighted_arithmetic_rel=weighted_arithmetic_rel,
        minimum_amplitude_abs2=min(
            (row.amplitude_abs2 for row in finite_rows), default=math.nan
        ),
        minimum_weighted_strength=min(
            (row.weighted_strength for row in finite_rows), default=math.nan
        ),
        sum_weighted_strength=math.fsum(
            row.weighted_strength for row in finite_rows
        ),
    )


def analyse_logs(directory: Path) -> LogAnalysis:
    candidates = [
        directory / "fleur_stdout.log",
        directory / "fleur_stderr.log",
        *sorted(directory.glob("job_environment_*.log")),
    ]
    files = [path for path in candidates if path.is_file()]
    text = "\n".join(
        path.read_text(encoding="utf-8", errors="replace") for path in files
    )
    reconstruction = re.search(
        r"RIXS\s+xx\s+spinor contribution-spectrum check:\s*"
        r"max abs diff\s*=\s*([+\-0-9.Ee]+)\s*"
        r"max rel diff\s*=\s*([+\-0-9.Ee]+)\s*"
        r"status\s*=\s*(PASS|FAIL)",
        text,
        flags=re.MULTILINE,
    )
    abort_pattern = re.compile(
        r"\bMPI_ABORT\b|\bMPI error\b|\bFLEUR-Error\b|"
        r"\bsegmentation fault\b|^\s*STOP\s+[1-9]\d*",
        flags=re.IGNORECASE | re.MULTILINE,
    )
    nonfinite_pattern = re.compile(
        r"(?<![A-Za-z])(?:NaN|[+\-]?Inf(?:inity)?)(?![A-Za-z])",
        flags=re.IGNORECASE,
    )
    distribution_pattern = re.compile(
        r".*(?:duplicated|duplicate|missing).*(?:k-point|k point).*",
        flags=re.IGNORECASE,
    )
    return LogAnalysis(
        files=files,
        has_spinor_summary=(
            "Spinor treatment         : first-variation coherent core-mj amplitude"
            in text
        ),
        has_completion="RIXS wrote spectrum to" in text,
        contribution_status=reconstruction.group(3) if reconstruction else "MISSING",
        reconstruction_abs=float(reconstruction.group(1)) if reconstruction else None,
        reconstruction_rel=float(reconstruction.group(2)) if reconstruction else None,
        abort_hits=abort_pattern.findall(text),
        nonfinite_hits=nonfinite_pattern.findall(text),
        distribution_warning_hits=distribution_pattern.findall(text),
    )


def analyse_run(
    label: str,
    directory: Path,
    expected_ranks: int,
    expected_ikpts: set[int],
) -> RunAnalysis:
    spectrum_path = locate_one(directory, "*_rixs.dat", "production RIXS spectrum")
    return RunAnalysis(
        label=label,
        directory=directory,
        expected_ranks=expected_ranks,
        spectrum=read_spectrum(spectrum_path),
        contribution=analyse_contributions(
            directory, expected_ranks, expected_ikpts
        ),
        logs=analyse_logs(directory),
    )


def format_set(values: set[int]) -> str:
    return "{" + ", ".join(str(value) for value in sorted(values)) + "}"


def close_row_error(maximum_absolute: float, maximum_relative: float) -> bool:
    return maximum_absolute < ROW_ABS_TOL or maximum_relative < ROW_REL_TOL


ManifoldContext = tuple[int, int, int]
BandEnergyKey = tuple[ManifoldContext, int]
ManifoldBlockKey = tuple[ManifoldContext, tuple[int, ...], tuple[int, ...]]


def band_energy_maps(
    rows: Sequence[ContributionRow],
) -> tuple[dict[BandEnergyKey, float], dict[BandEnergyKey, float]]:
    valence: dict[BandEnergyKey, float] = {}
    intermediate: dict[BandEnergyKey, float] = {}
    for row in rows:
        context = (row.ikpt, row.absorber_atom, row.absorber_type)
        valence.setdefault((context, row.band_v), row.eps_v)
        intermediate.setdefault((context, row.band_n), row.eps_n)
    return valence, intermediate


def construct_manifolds(
    energies: dict[BandEnergyKey, float],
) -> dict[BandEnergyKey, tuple[int, ...]]:
    membership: dict[BandEnergyKey, tuple[int, ...]] = {}
    contexts = sorted(context for context, _ in energies)
    for context in dict.fromkeys(contexts):
        ordered = sorted(
            (
                (band, energy)
                for (band_context, band), energy in energies.items()
                if band_context == context
            ),
            key=lambda item: (item[1], item[0]),
        )
        groups: list[list[tuple[int, float]]] = []
        for band, energy in ordered:
            if groups and energy - groups[-1][0][1] <= MANIFOLD_DEGENERACY_TOL:
                groups[-1].append((band, energy))
            else:
                groups.append([(band, energy)])
        for group in groups:
            bands = tuple(sorted(band for band, _ in group))
            for band in bands:
                membership[(context, band)] = bands
    return membership


def construct_manifold_blocks(
    rows: Sequence[ContributionRow],
    valence_membership: dict[BandEnergyKey, tuple[int, ...]],
    intermediate_membership: dict[BandEnergyKey, tuple[int, ...]],
) -> tuple[dict[ManifoldBlockKey, set[tuple[int, int, int, int, int]]], int]:
    blocks: dict[ManifoldBlockKey, set[tuple[int, int, int, int, int]]] = {}
    for row in rows:
        context = (row.ikpt, row.absorber_atom, row.absorber_type)
        block = (
            context,
            valence_membership[(context, row.band_v)],
            intermediate_membership[(context, row.band_n)],
        )
        blocks.setdefault(block, set()).add(row.transition_key)

    incomplete = 0
    for block, keys in blocks.items():
        context, valence_bands, intermediate_bands = block
        ikpt, absorber_atom, absorber_type = context
        expected = {
            (ikpt, band_v, band_n, absorber_atom, absorber_type)
            for band_v in valence_bands
            for band_n in intermediate_bands
        }
        incomplete += keys != expected
    return blocks, incomplete


def format_transition_key(key: tuple[int, int, int, int, int]) -> str:
    ikpt, band_v, band_n, absorber_atom, absorber_type = key
    return (
        f"(ikpt={ikpt}, band_v={band_v}, band_n={band_n}, "
        f"absorber_atom={absorber_atom}, absorber_type={absorber_type})"
    )


def compare_row_tables(
    reference: ContributionAnalysis,
    candidate: ContributionAnalysis,
    candidate_label: str,
) -> tuple[list[str], list[str]]:
    failures: list[str] = []
    warnings: list[str] = []
    if reference.duplicate_keys or candidate.duplicate_keys:
        failures.append(
            f"np1/{candidate_label} cannot compare contribution manifolds "
            "with duplicate transition identities"
        )
        return failures, warnings

    reference_rows = sorted(reference.rows, key=lambda row: row.transition_key)
    candidate_rows = sorted(candidate.rows, key=lambda row: row.transition_key)
    reference_keys = [row.transition_key for row in reference_rows]
    candidate_keys = [row.transition_key for row in candidate_rows]
    keys_equal = reference_keys == candidate_keys
    print(f"  np1 versus {candidate_label} row keys equal: {keys_equal}")
    if not keys_equal:
        failures.append(f"np1/{candidate_label} contribution row keys differ")
        first_mismatch = next(
            (
                (left, right)
                for left, right in zip(reference_keys, candidate_keys)
                if left != right
            ),
            None,
        )
        print(f"    first key mismatch: {first_mismatch}")
        return failures, warnings

    maximum_energy_difference = max(
        (
            abs(left.loss_energy - right.loss_energy)
            for left, right in zip(reference_rows, candidate_rows)
        ),
        default=0.0,
    )
    print(
        "    contribution energy comparison: "
        f"max_abs={maximum_energy_difference:.17e} "
        f"tolerance={CONTRIBUTION_ENERGY_ABS_TOL:.1e} (loss_energy_Ha)"
    )
    if (
        not math.isfinite(maximum_energy_difference)
        or maximum_energy_difference > CONTRIBUTION_ENERGY_ABS_TOL
    ):
        failures.append(
            f"np1/{candidate_label} contribution energies differ beyond tolerance"
        )

    fields = (
        "amplitude_real",
        "amplitude_imag",
        "amplitude_abs2",
        "weighted_strength",
    )
    metrics: dict[str, tuple[float, float]] = {}
    print("    raw band-labelled row comparisons (diagnostic only):")
    for field in fields:
        left = [getattr(row, field) for row in reference_rows]
        right = [getattr(row, field) for row in candidate_rows]
        metrics[field] = max_abs_and_relative(left, right)
        maximum_absolute, maximum_relative = metrics[field]
        worst_index = max(
            range(len(reference_rows)),
            key=lambda index: abs(left[index] - right[index]),
            default=None,
        )
        worst_key = (
            reference_rows[worst_index].transition_key
            if worst_index is not None
            else (0, 0, 0, 0, 0)
        )
        print(
            f"      {field}: max_abs={maximum_absolute:.17e} "
            f"max_rel_common={maximum_relative:.17e} "
            f"worst={format_transition_key(worst_key)}"
        )

    reference_valence, reference_intermediate = band_energy_maps(reference_rows)
    candidate_valence, candidate_intermediate = band_energy_maps(candidate_rows)
    reference_valence_membership = construct_manifolds(reference_valence)
    reference_intermediate_membership = construct_manifolds(reference_intermediate)
    candidate_valence_membership = construct_manifolds(candidate_valence)
    candidate_intermediate_membership = construct_manifolds(candidate_intermediate)
    membership_equal = (
        reference_valence_membership == candidate_valence_membership
        and reference_intermediate_membership == candidate_intermediate_membership
    )
    print(f"    degeneracy threshold: {MANIFOLD_DEGENERACY_TOL:.1e} Ha")
    print(f"    manifold membership equal: {membership_equal}")
    if not membership_equal:
        failures.append(f"np1/{candidate_label} manifold membership differs")
        return failures, warnings

    reference_blocks, reference_incomplete = construct_manifold_blocks(
        reference_rows,
        reference_valence_membership,
        reference_intermediate_membership,
    )
    candidate_blocks, candidate_incomplete = construct_manifold_blocks(
        candidate_rows,
        candidate_valence_membership,
        candidate_intermediate_membership,
    )
    print(
        "    incomplete Cartesian products: "
        f"np1={reference_incomplete} {candidate_label}={candidate_incomplete}"
    )
    if reference_incomplete or candidate_incomplete:
        failures.append(
            f"np1/{candidate_label} has incomplete manifold Cartesian products"
        )
        return failures, warnings
    if reference_blocks.keys() != candidate_blocks.keys():
        failures.append(f"np1/{candidate_label} manifold block identities differ")
        return failures, warnings

    reference_by_key = {row.transition_key: row for row in reference_rows}
    candidate_by_key = {row.transition_key: row for row in candidate_rows}
    manifold_fields = (
        ("amplitude_abs2", "manifold amplitude comparison"),
        ("weighted_strength", "manifold weighted comparison"),
    )
    for field, report_label in manifold_fields:
        block_keys = sorted(reference_blocks)
        left = [
            math.fsum(
                getattr(reference_by_key[key], field)
                for key in reference_blocks[block]
            )
            for block in block_keys
        ]
        right = [
            math.fsum(
                getattr(candidate_by_key[key], field)
                for key in candidate_blocks[block]
            )
            for block in block_keys
        ]
        maximum_absolute, maximum_relative = max_abs_and_relative(left, right)
        print(
            f"    {report_label}: max_abs={maximum_absolute:.17e} "
            f"max_rel_common={maximum_relative:.17e}"
        )
        if not passes_tolerance(maximum_absolute, maximum_relative):
            failures.append(
                f"np1/{candidate_label} manifold-summed {field} "
                "differs beyond tolerance"
            )
    return failures, warnings


def print_run(run: RunAnalysis) -> tuple[list[str], list[str]]:
    failures: list[str] = []
    warnings: list[str] = []
    spectrum = spectrum_metrics(run.spectrum)
    contribution = run.contribution
    logs = run.logs

    print(f"\n[{run.label}] {run.directory}")
    print(f"  spectrum: {run.spectrum.path.name}")
    print(
        "  spectrum metrics: "
        f"rows={spectrum['rows']} grid={spectrum['loss_min']:.17g}..."
        f"{spectrum['loss_max']:.17g} spacing={spectrum['spacing']:.17g} "
        f"trap={spectrum['trapezoidal_integral']:.17e} "
        f"rect={spectrum['rectangular_sum']:.17e}"
    )
    print(
        "  intensity: "
        f"min={spectrum['minimum_intensity']:.17e} "
        f"max={spectrum['maximum_intensity']:.17e} "
        f"peak_loss={spectrum['peak_loss']:.17e} "
        f"below_-1e-14={spectrum['below_negative_tolerance']} "
        f"nonfinite={spectrum['nonfinite_values']}"
    )
    print(
        f"  contribution files: found={len(contribution.files)} "
        f"expected={run.expected_ranks}"
    )
    print(
        "  rows per rank: "
        + ", ".join(
            f"{rank}:{len(rows)}"
            for rank, rows in sorted(contribution.rows_by_rank.items())
        )
    )
    print(f"  total contribution rows: {len(contribution.rows)}")
    print(
        f"  sum(weighted_strength): "
        f"{contribution.sum_weighted_strength:.17e}"
    )
    print(
        "  amplitude identity: "
        f"max_abs={contribution.amplitude_identity_abs:.17e} "
        f"max_rel={contribution.amplitude_identity_rel:.17e}"
    )
    print(
        "  weighted arithmetic: "
        f"max_abs={contribution.weighted_arithmetic_abs:.17e} "
        f"max_rel={contribution.weighted_arithmetic_rel:.17e}"
    )
    print(
        "  contribution minima: "
        f"amplitude_abs2={contribution.minimum_amplitude_abs2:.17e} "
        f"weighted_strength={contribution.minimum_weighted_strength:.17e} "
        f"nonfinite_rows={contribution.nonfinite_rows}"
    )
    for rank, ikpts in sorted(contribution.ikpts_by_rank.items()):
        print(f"  rank {rank:04d} ikpt set: {format_set(ikpts)}")
    print(f"  combined ikpt set: {format_set(contribution.combined_ikpts)}")
    print(f"  missing ikpts: {format_set(contribution.missing_ikpts)}")
    print(f"  unexpected ikpts: {format_set(contribution.unexpected_ikpts)}")
    print(f"  duplicate transition keys: {len(contribution.duplicate_keys)}")
    print(
        "  log checks: "
        f"spinor_summary={logs.has_spinor_summary} "
        f"completion={logs.has_completion} "
        f"contribution_status={logs.contribution_status} "
        f"reconstruction_abs={logs.reconstruction_abs} "
        f"reconstruction_rel={logs.reconstruction_rel}"
    )
    print(
        "  log diagnostics: "
        f"abort_hits={len(logs.abort_hits)} "
        f"nonfinite_hits={len(logs.nonfinite_hits)} "
        f"distribution_warnings={len(logs.distribution_warning_hits)}"
    )

    expected_rank_ids = set(range(run.expected_ranks))
    actual_rank_ids = set(contribution.rows_by_rank)
    if len(contribution.files) != run.expected_ranks:
        failures.append(f"{run.label}: wrong number of contribution files")
    if actual_rank_ids != expected_rank_ids:
        failures.append(
            f"{run.label}: contribution rank IDs {sorted(actual_rank_ids)} "
            f"do not match {sorted(expected_rank_ids)}"
        )
    if spectrum["max_spacing_deviation"] > GRID_TOL:
        failures.append(f"{run.label}: loss grid is not uniform")
    if spectrum["below_negative_tolerance"]:
        failures.append(f"{run.label}: spectrum has values below -1e-14")
    if spectrum["nonfinite_values"]:
        failures.append(f"{run.label}: spectrum has non-finite values")
    if contribution.nonfinite_rows:
        failures.append(f"{run.label}: contribution rows have non-finite values")
    if contribution.minimum_amplitude_abs2 < 0.0:
        failures.append(f"{run.label}: negative amplitude_abs2")
    if contribution.minimum_weighted_strength < 0.0:
        failures.append(f"{run.label}: negative weighted_strength")
    if not close_row_error(
        contribution.amplitude_identity_abs, contribution.amplitude_identity_rel
    ):
        failures.append(f"{run.label}: amplitude identity exceeds tolerance")
    if not close_row_error(
        contribution.weighted_arithmetic_abs,
        contribution.weighted_arithmetic_rel,
    ):
        failures.append(f"{run.label}: weighted arithmetic exceeds tolerance")
    if contribution.missing_ikpts:
        failures.append(f"{run.label}: missing expected k points")
    if contribution.unexpected_ikpts:
        failures.append(f"{run.label}: unexpected k points")
    if contribution.duplicate_keys:
        failures.append(f"{run.label}: duplicated transition keys")
    if not logs.has_spinor_summary:
        failures.append(f"{run.label}: spinor contraction summary missing")
    if not logs.has_completion:
        failures.append(f"{run.label}: RIXS completion marker missing")
    if logs.contribution_status != "PASS":
        failures.append(f"{run.label}: contribution-spectrum check is not PASS")
    if logs.abort_hits:
        failures.append(f"{run.label}: abort/error text found in logs")
    if logs.nonfinite_hits:
        failures.append(f"{run.label}: NaN/Inf text found in logs")
    if logs.distribution_warning_hits:
        warnings.append(f"{run.label}: k-point distribution warning text found")
    return failures, warnings


def compare_spectra(
    reference: RunAnalysis, candidate: RunAnalysis
) -> tuple[list[str], tuple[float, float]]:
    failures: list[str] = []
    grid_abs, _ = max_abs_and_relative(
        reference.spectrum.loss, candidate.spectrum.loss
    )
    maximum_absolute, maximum_relative = max_abs_and_relative(
        reference.spectrum.intensity, candidate.spectrum.intensity
    )
    print(
        f"  np1 versus {candidate.label}: rows "
        f"{len(reference.spectrum.loss)}/{len(candidate.spectrum.loss)}, "
        f"grid_max_abs={grid_abs:.17e}, "
        f"spectrum_max_abs={maximum_absolute:.17e}, "
        f"spectrum_max_rel_common={maximum_relative:.17e}"
    )
    if len(reference.spectrum.loss) != len(candidate.spectrum.loss):
        failures.append(f"np1/{candidate.label} spectrum row counts differ")
    if grid_abs > GRID_TOL:
        failures.append(f"np1/{candidate.label} loss grids differ")
    if not passes_tolerance(maximum_absolute, maximum_relative):
        failures.append(f"np1/{candidate.label} spectra differ beyond tolerance")
    return failures, (maximum_absolute, maximum_relative)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Validate SQAy spinor RIXS spectra and rank-local contribution "
            "tables for 1, 2, and 4 MPI ranks."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
        help="RIXS_MPI_validation directory (default: package root)",
    )
    parser.add_argument(
        "--expected-kpoints",
        type=int,
        default=27,
        help="Expected full-k ikpt count (default: 27)",
    )
    args = parser.parse_args()

    root = args.root.resolve()
    expected_ikpts = set(range(1, args.expected_kpoints + 1))
    failures: list[str] = []
    warnings: list[str] = []
    analyses: dict[str, RunAnalysis] = {}

    try:
        for label, directory_name, expected_ranks in RUNS:
            directory = root / directory_name
            analyses[label] = analyse_run(
                label, directory, expected_ranks, expected_ikpts
            )
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    print("Spinor RIXS MPI validation report")
    print(f"Package root: {root}")
    print(
        "Acceptance tolerances: "
        f"spectrum abs<{SPECTRUM_ABS_TOL:g} OR rel<{SPECTRUM_REL_TOL:g}; "
        f"row arithmetic abs<{ROW_ABS_TOL:g} OR rel<{ROW_REL_TOL:g}; "
        f"contribution energy abs<={CONTRIBUTION_ENERGY_ABS_TOL:g}"
    )

    for label, _, _ in RUNS:
        run_failures, run_warnings = print_run(analyses[label])
        failures.extend(run_failures)
        warnings.extend(run_warnings)

    row_counts = {label: len(run.contribution.rows) for label, run in analyses.items()}
    if len(set(row_counts.values())) != 1:
        failures.append(f"combined contribution row counts differ: {row_counts}")

    print("\n[Production spectrum comparisons]")
    for label in ("np2", "np4"):
        comparison_failures, _ = compare_spectra(analyses["np1"], analyses[label])
        failures.extend(comparison_failures)

    print("\n[Contribution-table comparisons]")
    for label in ("np2", "np4"):
        comparison_failures, comparison_warnings = compare_row_tables(
            analyses["np1"].contribution,
            analyses[label].contribution,
            label,
        )
        failures.extend(comparison_failures)
        warnings.extend(comparison_warnings)

    reference_path = root / "reference" / "SQAy_stage2_serial_spectrum.dat"
    if reference_path.is_file():
        print("\n[Committed Stage 2 serial-reference comparisons]")
        serial_reference = read_spectrum(reference_path)
        for label in ("np1", "np2", "np4"):
            grid_abs, _ = max_abs_and_relative(
                serial_reference.loss, analyses[label].spectrum.loss
            )
            maximum_absolute, maximum_relative = max_abs_and_relative(
                serial_reference.intensity, analyses[label].spectrum.intensity
            )
            print(
                f"  reference versus {label}: grid_max_abs={grid_abs:.17e} "
                f"spectrum_max_abs={maximum_absolute:.17e} "
                f"spectrum_max_rel_common={maximum_relative:.17e}"
            )
            if grid_abs > GRID_TOL:
                failures.append(f"serial reference/{label} loss grids differ")
            if not passes_tolerance(maximum_absolute, maximum_relative):
                failures.append(
                    f"serial reference/{label} spectra differ beyond tolerance"
                )

    print("\n[Warnings]")
    if warnings:
        for warning in warnings:
            print(f"  WARNING: {warning}")
    else:
        print("  none")

    print("\n[Acceptance]")
    if failures:
        for failure in failures:
            print(f"  FAIL: {failure}")
        print("OVERALL STATUS: FAIL")
        return 1
    print("  all required checks passed")
    print("OVERALL STATUS: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
