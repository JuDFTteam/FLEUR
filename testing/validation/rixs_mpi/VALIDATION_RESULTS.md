# Validated RIXS MPI checkpoint

Spinor RIXS is validated in serial and across 1-, 2-, and 4-rank pure k-point
MPI layouts, with additional successful validation of the shared-k-point
subgroup-root path. Spectra, transition coverage, reconstruction, and
degenerate-manifold-summed contributions are MPI invariant.

## Pure k-point validation

The 1-, 2-, and 4-rank runs used `-pe_per_kpt 1` and produced:

- 50,784 combined contribution rows in every layout;
- complete coverage of all 27 k points;
- zero duplicate transition identities;
- passing contribution-spectrum reconstruction and row-arithmetic checks;
- no abort, non-finite, or k-point-distribution diagnostics.

Production spectra agreed with np1 to:

| Comparison | Maximum absolute difference | Maximum relative difference |
|---|---:|---:|
| np1/np2 | 3.7947076036992655e-19 | 5.092959261543016e-15 |
| np1/np4 | 5.5565361339882102e-19 | 7.457547490116559e-15 |

The maximum absolute differences from the preserved serial reference were
`1.8702487475374952e-18`, `2.0870891820345960e-18`, and
`2.3852447794681098e-18` for np1, np2, and np4 respectively.

With the fixed `1e-10` Ha degeneracy threshold, manifold-summed contributions
agreed to:

| Comparison | `amplitude_abs2` max abs | `weighted_strength` max abs |
|---|---:|---:|
| np1/np2 | 5.0989671327600507e-13 | 1.2742852921207956e-14 |
| np1/np4 | 6.2514102255911046e-13 | 1.5622934682642996e-14 |

Manifold membership was identical across layouts and every manifold block had
a complete valence × intermediate Cartesian product.

## Shared-k-point validation

The shared-k-point runs assigned multiple MPI ranks to the same k-point
subgroup. Only `fmpi%n_rank == 0` evaluated RIXS transitions, while every rank
participated in the unchanged global reductions. The spectra and combined
transition rows remained invariant, and the former factor-of-N duplication was
absent.

## Interpretation

Raw band-labelled rows can differ because independently regenerated
eigenvectors may rotate within degenerate valence and intermediate subspaces.
These raw differences remain diagnostic output. Acceptance uses complete
manifold sums derived independently from `eps_v_Ha` and `eps_n_Ha`; it does not
group transitions by loss energy.
