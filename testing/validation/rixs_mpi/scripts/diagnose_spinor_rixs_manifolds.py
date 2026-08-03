#!/usr/bin/env python3
"""Diagnose MPI-layout sensitivity of band-labelled RIXS contributions."""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from compare_spinor_rixs_mpi import (
    RUNS,
    ContributionAnalysis,
    ContributionRow,
    analyse_contributions,
)


THRESHOLDS = (1.0e-12, 1.0e-10, 1.0e-8, 1.0e-6)
FIELDS = ("amplitude_abs2", "weighted_strength")
NEIGHBOR_COUNT = 2
NEAR_TRANSITION_TOL = 1.0e-8

TransitionKey = tuple[int, int, int, int, int]
BlockKey = tuple[int, int, int, tuple[int, ...], tuple[int, ...]]


def row_map(analysis: ContributionAnalysis) -> dict[TransitionKey, ContributionRow]:
    result: dict[TransitionKey, ContributionRow] = {}
    for row in analysis.rows:
        if row.transition_key in result:
            raise ValueError(f"duplicate transition key {row.transition_key}")
        result[row.transition_key] = row
    return result


def signed_delta(reference: ContributionRow, candidate: ContributionRow, field: str) -> float:
    return getattr(candidate, field) - getattr(reference, field)


def worst_difference(
    reference: dict[TransitionKey, ContributionRow],
    candidate: dict[TransitionKey, ContributionRow],
    field: str,
) -> tuple[float, TransitionKey]:
    return max(
        (abs(signed_delta(reference[key], candidate[key], field)), key)
        for key in reference
    )


def band_energy_maps(
    rows: Iterable[ContributionRow],
) -> tuple[dict[tuple[int, int], float], dict[tuple[int, int], float]]:
    valence_values: dict[tuple[int, int], list[float]] = defaultdict(list)
    intermediate_values: dict[tuple[int, int], list[float]] = defaultdict(list)
    for row in rows:
        valence_values[(row.ikpt, row.band_v)].append(row.eps_v)
        intermediate_values[(row.ikpt, row.band_n)].append(row.eps_n)

    def collapse(values: dict[tuple[int, int], list[float]], label: str) -> dict[tuple[int, int], float]:
        collapsed: dict[tuple[int, int], float] = {}
        for key, samples in values.items():
            spread = max(samples) - min(samples)
            if spread > 1.0e-14:
                raise ValueError(
                    f"inconsistent repeated {label} energy for {key}: spread={spread:.3e} Ha"
                )
            collapsed[key] = samples[0]
        return collapsed

    return collapse(valence_values, "valence"), collapse(intermediate_values, "intermediate")


def neighboring_bands(
    energies: dict[tuple[int, int], float], ikpt: int, target_band: int
) -> list[int]:
    bands = sorted(band for point, band in energies if point == ikpt)
    if target_band not in bands:
        return []
    position = bands.index(target_band)
    return bands[
        max(0, position - NEIGHBOR_COUNT) : position + NEIGHBOR_COUNT + 1
    ]


def report_fixed_band_neighborhood(
    reference: dict[TransitionKey, ContributionRow],
    candidate: dict[TransitionKey, ContributionRow],
    reference_energies: dict[tuple[int, int], float],
    candidate_energies: dict[tuple[int, int], float],
    key: TransitionKey,
    vary_valence: bool,
) -> None:
    ikpt, band_v, band_n, absorber_atom, absorber_type = key
    target_band = band_v if vary_valence else band_n
    varying_energies = reference_energies
    target_energy = reference_energies[(ikpt, target_band)]
    bands = neighboring_bands(varying_energies, ikpt, target_band)
    label = "valence" if vary_valence else "intermediate"
    print(f"      neighboring {label} bands (fixed other band):")
    amplitude_sum = 0.0
    weighted_sum = 0.0
    for band in bands:
        neighbor_key = (
            (ikpt, band, band_n, absorber_atom, absorber_type)
            if vary_valence
            else (ikpt, band_v, band, absorber_atom, absorber_type)
        )
        energy_reference = reference_energies[(ikpt, band)]
        energy_candidate = candidate_energies.get((ikpt, band), math.nan)
        if neighbor_key not in reference or neighbor_key not in candidate:
            print(
                f"        band={band:4d} eps_np1={energy_reference:.17e} "
                f"eps_candidate={energy_candidate:.17e} "
                f"gap_to_target={energy_reference - target_energy:+.3e} row=absent"
            )
            continue
        amplitude_delta = signed_delta(
            reference[neighbor_key], candidate[neighbor_key], "amplitude_abs2"
        )
        weighted_delta = signed_delta(
            reference[neighbor_key], candidate[neighbor_key], "weighted_strength"
        )
        amplitude_sum += amplitude_delta
        weighted_sum += weighted_delta
        print(
            f"        band={band:4d} eps_np1={energy_reference:.17e} "
            f"eps_candidate={energy_candidate:.17e} "
            f"gap_to_target={energy_reference - target_energy:+.3e} "
            f"delta_amp2={amplitude_delta:+.17e} "
            f"delta_weighted={weighted_delta:+.17e}"
        )
    print(
        f"        neighborhood signed sums: delta_amp2={amplitude_sum:+.17e} "
        f"delta_weighted={weighted_sum:+.17e}"
    )


def report_near_transition_rows(
    reference: dict[TransitionKey, ContributionRow],
    candidate: dict[TransitionKey, ContributionRow],
    target_key: TransitionKey,
) -> None:
    target = reference[target_key]
    nearby = [
        key
        for key, row in reference.items()
        if row.ikpt == target.ikpt
        and row.absorber_atom == target.absorber_atom
        and row.absorber_type == target.absorber_type
        and abs(row.loss_energy - target.loss_energy) <= NEAR_TRANSITION_TOL
    ]
    nearby.sort(
        key=lambda key: abs(
            signed_delta(reference[key], candidate[key], "amplitude_abs2")
        ),
        reverse=True,
    )
    amplitude_sum = math.fsum(
        signed_delta(reference[key], candidate[key], "amplitude_abs2")
        for key in nearby
    )
    weighted_sum = math.fsum(
        signed_delta(reference[key], candidate[key], "weighted_strength")
        for key in nearby
    )
    print(
        f"      near-equal transition-energy rows: count={len(nearby)} "
        f"window={NEAR_TRANSITION_TOL:.1e} Ha "
        f"signed_delta_amp2={amplitude_sum:+.17e} "
        f"signed_delta_weighted={weighted_sum:+.17e}"
    )
    for key in nearby[:12]:
        row = reference[key]
        print(
            f"        v={row.band_v:4d} n={row.band_n:4d} "
            f"loss_np1={row.loss_energy:.17e} "
            f"delta_loss={candidate[key].loss_energy - row.loss_energy:+.3e} "
            f"delta_amp2={signed_delta(row, candidate[key], 'amplitude_abs2'):+.17e} "
            f"delta_weighted={signed_delta(row, candidate[key], 'weighted_strength'):+.17e}"
        )


def report_worst_row(
    reference: dict[TransitionKey, ContributionRow],
    candidate: dict[TransitionKey, ContributionRow],
    reference_valence: dict[tuple[int, int], float],
    candidate_valence: dict[tuple[int, int], float],
    reference_intermediate: dict[tuple[int, int], float],
    candidate_intermediate: dict[tuple[int, int], float],
    candidate_label: str,
    field: str,
) -> None:
    maximum, key = worst_difference(reference, candidate, field)
    left = reference[key]
    right = candidate[key]
    print(f"\n    worst {field}: max_abs={maximum:.17e}")
    print(
        f"      identity: ikpt={left.ikpt} band_v={left.band_v} band_n={left.band_n} "
        f"absorber_atom={left.absorber_atom} absorber_type={left.absorber_type}"
    )
    print(
        f"      loss_energy_Ha: np1={left.loss_energy:.17e} "
        f"{candidate_label}={right.loss_energy:.17e} "
        f"delta={right.loss_energy - left.loss_energy:+.17e}"
    )
    print(
        f"      amplitude_abs2: np1={left.amplitude_abs2:.17e} "
        f"{candidate_label}={right.amplitude_abs2:.17e} "
        f"delta={right.amplitude_abs2 - left.amplitude_abs2:+.17e}"
    )
    print(
        f"      weighted_strength: np1={left.weighted_strength:.17e} "
        f"{candidate_label}={right.weighted_strength:.17e} "
        f"delta={right.weighted_strength - left.weighted_strength:+.17e}"
    )
    report_fixed_band_neighborhood(
        reference,
        candidate,
        reference_valence,
        candidate_valence,
        key,
        vary_valence=True,
    )
    report_fixed_band_neighborhood(
        reference,
        candidate,
        reference_intermediate,
        candidate_intermediate,
        key,
        vary_valence=False,
    )
    report_near_transition_rows(reference, candidate, key)


def cluster_bands(
    energies: dict[tuple[int, int], float], threshold: float
) -> dict[tuple[int, int], tuple[int, ...]]:
    groups_by_band: dict[tuple[int, int], tuple[int, ...]] = {}
    ikpts = sorted({ikpt for ikpt, _ in energies})
    for ikpt in ikpts:
        ordered = sorted(
            (
                (band, energy)
                for (point, band), energy in energies.items()
                if point == ikpt
            ),
            key=lambda item: (item[1], item[0]),
        )
        groups: list[list[tuple[int, float]]] = []
        for band, energy in ordered:
            if groups and energy - groups[-1][0][1] <= threshold:
                groups[-1].append((band, energy))
            else:
                groups.append([(band, energy)])
        for group in groups:
            bands = tuple(sorted(band for band, _ in group))
            for band in bands:
                groups_by_band[(ikpt, band)] = bands
    return groups_by_band


def manifold_blocks(
    keys: Iterable[TransitionKey],
    valence_groups: dict[tuple[int, int], tuple[int, ...]],
    intermediate_groups: dict[tuple[int, int], tuple[int, ...]],
) -> dict[BlockKey, set[TransitionKey]]:
    blocks: dict[BlockKey, set[TransitionKey]] = defaultdict(set)
    for key in keys:
        ikpt, band_v, band_n, absorber_atom, absorber_type = key
        block = (
            ikpt,
            absorber_atom,
            absorber_type,
            valence_groups[(ikpt, band_v)],
            intermediate_groups[(ikpt, band_n)],
        )
        blocks[block].add(key)
    return blocks


def group_statistics(
    groups: dict[tuple[int, int], tuple[int, ...]],
    energies: dict[tuple[int, int], float],
) -> tuple[int, int, float]:
    unique = {(ikpt, bands) for (ikpt, _), bands in groups.items()}
    nonsingleton = [(ikpt, bands) for ikpt, bands in unique if len(bands) > 1]
    maximum_spread = max(
        (
            max(energies[(ikpt, band)] for band in bands)
            - min(energies[(ikpt, band)] for band in bands)
            for ikpt, bands in nonsingleton
        ),
        default=0.0,
    )
    return len(unique), len(nonsingleton), maximum_spread


def maximum_energy_map_difference(
    reference: dict[tuple[int, int], float],
    candidate: dict[tuple[int, int], float],
) -> float:
    if reference.keys() != candidate.keys():
        return math.inf
    return max(
        (abs(reference[key] - candidate[key]) for key in reference),
        default=0.0,
    )


def report_manifold_comparison(
    reference: dict[TransitionKey, ContributionRow],
    candidate: dict[TransitionKey, ContributionRow],
    reference_valence: dict[tuple[int, int], float],
    reference_intermediate: dict[tuple[int, int], float],
    candidate_valence: dict[tuple[int, int], float],
    candidate_intermediate: dict[tuple[int, int], float],
    candidate_label: str,
    threshold: float,
) -> None:
    valence_groups = cluster_bands(reference_valence, threshold)
    intermediate_groups = cluster_bands(reference_intermediate, threshold)
    candidate_valence_groups = cluster_bands(candidate_valence, threshold)
    candidate_intermediate_groups = cluster_bands(candidate_intermediate, threshold)
    valence_membership_changes = sum(
        valence_groups[key] != candidate_valence_groups[key]
        for key in valence_groups
    )
    intermediate_membership_changes = sum(
        intermediate_groups[key] != candidate_intermediate_groups[key]
        for key in intermediate_groups
    )
    blocks = manifold_blocks(reference, valence_groups, intermediate_groups)
    complete_blocks: dict[BlockKey, set[TransitionKey]] = {}
    incomplete_blocks = 0
    for block, keys in blocks.items():
        ikpt, absorber_atom, absorber_type, valence_bands, intermediate_bands = block
        expected = {
            (ikpt, band_v, band_n, absorber_atom, absorber_type)
            for band_v in valence_bands
            for band_n in intermediate_bands
        }
        if keys == expected:
            complete_blocks[block] = keys
        else:
            incomplete_blocks += 1

    valence_stats = group_statistics(valence_groups, reference_valence)
    intermediate_stats = group_statistics(intermediate_groups, reference_intermediate)
    print(f"\n  threshold={threshold:.1e} Ha")
    print(
        "    valence groups: "
        f"total={valence_stats[0]} nonsingleton={valence_stats[1]} "
        f"max_spread={valence_stats[2]:.3e} Ha"
    )
    print(
        "    intermediate groups: "
        f"total={intermediate_stats[0]} nonsingleton={intermediate_stats[1]} "
        f"max_spread={intermediate_stats[2]:.3e} Ha"
    )
    print(
        "    cross-layout group-membership changes: "
        f"valence_states={valence_membership_changes} "
        f"intermediate_states={intermediate_membership_changes}"
    )
    print(
        f"    Cartesian blocks: complete={len(complete_blocks)} "
        f"incomplete_excluded={incomplete_blocks}"
    )

    for field in FIELDS:
        comparisons: list[tuple[float, BlockKey, float, float]] = []
        for block, keys in complete_blocks.items():
            reference_sum = math.fsum(getattr(reference[key], field) for key in keys)
            candidate_sum = math.fsum(getattr(candidate[key], field) for key in keys)
            comparisons.append(
                (abs(candidate_sum - reference_sum), block, reference_sum, candidate_sum)
            )
        maximum, worst_block, reference_sum, candidate_sum = max(
            comparisons, default=(0.0, (0, 0, 0, (), ()), 0.0, 0.0)
        )
        common_maximum = max(
            (
                max(abs(left_sum), abs(right_sum))
                for _, _, left_sum, right_sum in comparisons
            ),
            default=0.0,
        )
        relative_common = maximum / max(common_maximum, math.ulp(0.0))
        ikpt, absorber_atom, absorber_type, valence_bands, intermediate_bands = worst_block
        print(
            f"    {field}: max_abs={maximum:.17e} "
            f"max_rel_common={relative_common:.17e}"
        )
        print(
            f"      worst block: ikpt={ikpt} absorber_atom={absorber_atom} "
            f"absorber_type={absorber_type} valence={valence_bands} "
            f"intermediate={intermediate_bands} np1={reference_sum:.17e} "
            f"{candidate_label}={candidate_sum:.17e}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Diagnose degenerate-band redistribution in spinor RIXS contribution tables."
    )
    parser.add_argument(
        "--root",
        type=Path,
        required=True,
        help="RIXS validation result root containing SQAy_np1/np2/np4",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    expected_ikpts = set(range(1, 28))

    analyses = {
        label: analyse_contributions(root / directory, ranks, expected_ikpts)
        for label, directory, ranks in RUNS
    }
    maps = {label: row_map(analysis) for label, analysis in analyses.items()}
    if not (maps["np1"].keys() == maps["np2"].keys() == maps["np4"].keys()):
        raise ValueError("discrete transition keys differ across MPI layouts")
    energy_maps = {
        label: band_energy_maps(analysis.rows) for label, analysis in analyses.items()
    }

    print("Spinor RIXS degenerate-manifold diagnostic")
    print(f"Result root: {root}")
    print(f"Discrete transition rows: {len(maps['np1'])}")
    print("Available eigenvalue columns: eps_v_Ha and eps_n_Ha (sufficient for manifold construction)")

    for candidate_label in ("np2", "np4"):
        print(f"\n[np1 versus {candidate_label}: worst rows and local compensation]")
        print(
            "  corresponding eigenvalue drift: "
            f"valence_max_abs={maximum_energy_map_difference(energy_maps['np1'][0], energy_maps[candidate_label][0]):.17e} Ha "
            f"intermediate_max_abs={maximum_energy_map_difference(energy_maps['np1'][1], energy_maps[candidate_label][1]):.17e} Ha"
        )
        for field in FIELDS:
            report_worst_row(
                maps["np1"],
                maps[candidate_label],
                energy_maps["np1"][0],
                energy_maps[candidate_label][0],
                energy_maps["np1"][1],
                energy_maps[candidate_label][1],
                candidate_label,
                field,
            )

        print(f"\n[np1 versus {candidate_label}: Cartesian manifold sums]")
        for threshold in THRESHOLDS:
            report_manifold_comparison(
                maps["np1"],
                maps[candidate_label],
                energy_maps["np1"][0],
                energy_maps["np1"][1],
                energy_maps[candidate_label][0],
                energy_maps[candidate_label][1],
                candidate_label,
                threshold,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
