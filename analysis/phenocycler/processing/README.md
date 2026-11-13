# Islet inflammation ordering

The input is the pancreas release that excludes EpCAM from protein typing. The outcome is an ordering of islets within each donor by observed local immune-cell density. Native islet membership and the established exterior 20 µm band are reused. The pancreas tissue mask excludes pLN, and the exterior band excludes every islet core and endocrine cluster. Bands from two different islets can overlap, as they do in the native membership table; each islet is measured separately.

For each islet, count admitted cells with final primary type **Immune** in the core and surrounding band. Add the counts and divide by the sum of the two assessed areas, in mm². This is an area-weighted density, not the sum of two densities. Both regions must have positive area and an admitted cell; an assessed region can have zero immune cells. Within a donor, average the ranks of exactly tied densities and rescale the smallest and largest rank to 0 and 1. All measured zeros share 0. An all-zero donor has rank 0 throughout; a donor with only one distinct positive density and no zero has no relative rank. The smallest positive value may be rank 0 if that donor has no measured zero.

The primary fields are `combined_immune_per_mm2`, `pseudotime_rank`, and `density_eligible` in the islet AnnData; `islet_pseudotime_rank` in cell AnnData inherits it only through a unique region link. This order is not an estimate of elapsed time. The original composition graph is available as `X_islet_graph`, `composition_dpt`, and `composition_dpt_rank`.

To rebuild from the private full-section features, memberships, native image calibration, and masks:

```sh
python islet_inflammation.py measure --islets "$original/pancreas_islets.h5ad" \
  --features-root "$native_features" --manifest "$image_manifest" --output "$staging"
python islet_inflammation.py attach --source "$original" \
  --islets "$staging/pancreas_islets.h5ad" --output "$staging"
```

The optional `--geometry-root` and `--tissue-root` relocate the mask files by basename. Use a new private staging directory; apply its contents only after the complete run succeeds. To attach an already verified islet result at a second location, copy the islet AnnData, donor summary CSV, islet ranks Parquet and islet method JSON to staging, then run only `attach` against the original cell files. The file manifest records pre/post cell hashes and linkage totals. The typing release records the hashes of the producer, islet data, method and cell-file manifest. Untouched assay tables and pLN files stay with the original typing release.

## Dedicated pLN compartment sensitivity

Use the Senior Python environment with its `core/multiplex5000` package on the import path. The Stage 03 roster identifies the native export for each dedicated section. The producer verifies each object ID and nucleus admission, reproduces the released whole-cell marker calls, and independently fits the established section gates to unshifted whole-cell and own-compartment means. EpCAM is excluded by the tissue rule manifest.

```bash
python pln_compartment_gates.py \
  --exports "$qupath_cells" --release "$dedicated_release" --calls "$epcam_excluded_calls" \
  --panel "$pln_manifest" --registry "$compartment_registry" \
  --out "$private_review" --scanpy "$scanpy_bundle"
```

The two compact Scanpy tables report marker transitions, measured cell counts, fitted thresholds and complete root qualification by donor. The per-object tables under `$private_review/unmatched` retain evidence for the separate cluster comparison. This sensitivity uses raw compartment intensities. The dedicated pLN release does not include compensated REDSEA intensities.

## Dedicated pLN candidate scoring and marker evidence

The continuous-score comparison uses SCIMAP 2.3.5 and the Senior dedicated pLN tissue manifest. It preserves the assay-specific unavailable branches and native cell IDs. Run `scimap_root_sensitivity.py` for the stratified comparison followed by `scimap_all_cells.py` for complete sections. Run `scimap_anchor_audit.py` to check the actual positivity of defining markers in formerly unmatched cells. Each command accepts explicit private input and output paths via `--help`. For paired donor/marker evidence, `pln_unmatched_review_markers.py` accepts `--mode review` or `--mode scimap`.

The separate October cluster- and cell-marker review provides broad immune evidence. These annotated review labels, strict rule outcomes and continuous-score candidates remain individually available in the private Scanpy data and notebook 05.
