# Trajectory parameters: every one, with where it came from

Audited 2026-09-15 on the maintainer's instruction that **no pre-existing trajectory parameter
is used without explicit approval**. 189 parameters were extracted from the five trajectory
sources and the app module, then every provenance claim was checked against the file, the git
history and the installed library rather than taken from a comment. **57 claims did not
survive that check.**

Provenance means: `MEASURED_HERE` requires a number, a cohort and a procedure. A comment saying
a value is 'chosen for stability' or 'works well' is not a measurement and was downgraded.

| bucket | n | meaning |
|---|---|---|
| live | 98 | active in a per-donor, blind-rooted, no-Harmony ordering |
| dormant | 67 | dies with Harmony, insulin-rooting or cohort pooling |
| dead | 24 | unreachable in this path |

Of the live ones: **53 ARBITRARY**, **24 INHERITED_DEFAULT**, **16 LIBRARY_DEFAULT**, **4 MEASURED_HERE**, **1 UNKNOWN**.


## Live — these need a ruling


### `apps/islet/R/00_globals.R`

**`EXCLUDED_MARKERS` = `c("DAPI", "IAPP", "Ker8_18", "Keratin_5")`** — apps/islet/R/00_globals.R:273 (applied at mod_trajectory_server.R:235, 261, 969, 1503) — *ARBITRARY* **[claim overturned]**  
Controls: Four markers are removed from every trajectory selector, the multi-marker heatmap default and the drill-down colour-by list.  
Checked: The cited authority does not exist and does not match: there is no scripts/senior/marker_taxonomy.py (the module is core/phenocycler/phenocycler/marker_taxonomy.py), and evaluating its EXCLUDED gives ['CD11b','DAPI','IAPP'] -- so the app hides Ker8_18 and Keratin_5 that upstream does not exclude, and still offers CD11b that upstream does; this is a hand-maintained list that has drifted from the source it claims to mirror.


### `apps/islet/R/mod_trajectory_server.R`

**`obs column whitelist (cols)` = `c("combined_islet_id","base_islet_id","donor_status","Case ID","case_id","imageid","dpt_pseudotime","dpt_pseudotime_combined","age","gender","total_cells")`** — apps/islet/R/mod_trajectory_server.R:126-131 — *ARBITRARY*  
Controls: Which obs columns survive into the app; anything not listed is invisible to every downstream card.  
Checked: The comment only describes the two pseudotime columns and no commit explains the list; measured against the deployed file it drops five columns that exist in obs -- acq_run, stain_group, endocrine_fraction, total_cells_peri and pseudotime_unavailable_reason -- so both finer batch variables and the peri cell count are unreachable from this tab.

**`UMAP feature colour scale` = `scale_color_viridis_c(option = "inferno", limits = range(df$value), na.value = "#bbbbbb")`** — apps/islet/R/mod_trajectory_server.R:1340-1342 — *ARBITRARY*  
Controls: Colour ramp for the selected feature on the UMAP, rescaled to the currently filtered data's min/max -- so the same islet changes colour when the cohort filter changes.  
Checked: "Continuous colormap scaled to data min/max" describes the behaviour without justifying either the palette or the rescaling, and the rescaling is real -- limits comes from the post-filter frame.

**`INTENSITY_FLOOR` = `0.50`** — apps/islet/R/mod_trajectory_server.R:1368 (used 1420-1423) — *ARBITRARY*  
Controls: Minimum colour saturation of a pseudotime strip bin, so the dimmest bin still reads as its donor colour -- compressing the displayed feature range into the upper half of the ramp.  
Checked: Commit 77a9601 introduced it and only restates the design ("5-95% clip, 50% floor"); no measurement, and the inline comment is a preference ("bars never go pale-white").

**`pseudotime strip bin count` = `n_bins <- 25 (breaks seq(0, 1, length.out = 26))`** — apps/islet/R/mod_trajectory_server.R:1375-1378 — *ARBITRARY*  
Controls: Resolution of the dominant-status strip under the scatter.  
Checked: Introduced in 77a9601 with no rationale for 25 in either the commit or the code.

**`dominant-status normalization by status totals` = `status_totals over c("ND","Aab+","T1D"); dens <- counts / status_totals; dominant <- names(dens)[which.max(dens)]`** — apps/islet/R/mod_trajectory_server.R:1385, 1399, 1401 — *ARBITRARY* **[claim overturned]**  
Controls: Each bin's hue is the status most enriched relative to that status's whole-cohort pool, so the strip is defined only across pooled donors.  
Checked: I found the source of the quoted counts and they are stale: ND 2300 / Aab+ 1668 / T1D 1246 sums to 5,214, the legacy CODEX trajectory file documented at README.md:168, whereas the deployed run-20 file has ND 6,303 / T1D 999 / Aab+ 0 (3,758 / 613 / 0 with finite pseudotime) -- status_totals is computed live so only the comment misleads, but the rule itself is a design choice with no measurement, and with ND outnumbering T1D six to one a bin holding one T1D islet can read as T1D-dominant.

**`strip intensity percentile clip` = `stats::quantile(hm$feat_mean, c(0.05, 0.95))`** — apps/islet/R/mod_trajectory_server.R:1412-1419 — *ARBITRARY*  
Controls: Bins below the 5th and above the 95th percentile of bin means are saturated to the ends of the tint ramp.  
Checked: "robust 5-95th percentile clip" names the rule without justifying the cut points; commit 77a9601 likewise only restates it.

**`multi-heatmap default marker selection` = `all core marker_mean_* columns except EXCLUDED_MARKERS, inline = TRUE`** — apps/islet/R/mod_trajectory_server.R:1499-1506 — *ARBITRARY* **[claim overturned]**  
Controls: Which markers the heatmap draws before any user choice -- every non-excluded core marker at once.  
Checked: Two different decisions are bundled: the core-only half is supported by the verified 116/58 column count, but selecting them ALL by default comes from commit c883cad ("Select all 31 markers by default (was 8)"), which cites only the resulting pixel height; measured on this file 55 boxes tick by default (58 core minus DAPI, Ker8_18, Keratin_5 -- IAPP is not present), giving a 2,040 px heatmap.

**`Case ID regex` = `stringr::str_extract(combined_islet_id, "^[0-9]{3,4}")`** — apps/islet/R/mod_trajectory_server.R:154 — *ARBITRARY* **[claim overturned]**  
Controls: Parses the donor id out of combined_islet_id when imageid is absent; a donor id outside 3-4 digits is silently not recovered.  
Checked: Not unknown so much as unrecorded: it predates the module split (first appears in 230af4d, the 5,397-line app.R extraction) with no commit rationale; I measured the deployed ids and they happen to be 3 or 4 digits ("115" plus nine 6xxx donors), and because imageid is present the coalesce at line 152 means this branch never actually fires today.

**`multi-heatmap rank normalization (second copy)` = `result$pt <- rank(result$pt, ties.method="average") / nrow(result)`** — apps/islet/R/mod_trajectory_server.R:1540-1541 — *ARBITRARY*  
Controls: Re-ranks pseudotime over the multi-heatmap's own frame, which is NOT cohort-filtered, so the heatmap x-axis and the scatter x-axis can disagree.  
Checked: "matches scatter x-axis" is the only note and it is not true: I read traj_multi_feature_data() (lines 1510-1545) end to end and it applies only traj_min_cells -- there is no cohort$donor_ids() filter -- so the two ranks are computed over different row sets whenever the cohort is narrowed.

**`multi-heatmap bin minimum support` = `if (length(bin_data) >= 1) mean(...) else NA`** — apps/islet/R/mod_trajectory_server.R:1573 — *ARBITRARY*  
Controls: A tile drawn from a single islet is rendered identically to one drawn from hundreds; only an empty bin becomes NA.  
Checked: No comment or commit addresses the floor, and nothing downstream carries the per-bin count into the tile (bin_counts is computed at 1572 and never used).

**`multi-heatmap minimum bins per marker` = `if (length(valid_means) < 2) next`** — apps/islet/R/mod_trajectory_server.R:1577 — *ARBITRARY*  
Controls: A marker with fewer than two occupied bins is dropped from the heatmap without notice.  
Checked: Unrecorded; 2 is simply the minimum for sd() to return a finite value on the next line, and the drop is silent -- the marker stays ticked in the selector but has no row.

**`multi-heatmap z-score clamp and scale limits` = `sd forced to 1 when zero; z <- pmax(-2.5, pmin(2.5, z)); limits = c(-2.5, 2.5)`** — apps/islet/R/mod_trajectory_server.R:1581-1582 and 1609 — *ARBITRARY*  
Controls: Row-wise z across bins is clipped at +/-2.5 SD, so a marker with a sharp excursion looks the same as one that just reaches the cap.  
Checked: The only annotation is the word "clamp"; no commit gives a basis for 2.5, and the sd <- 1 substitution silently turns a flat marker's zeros into z = 0 rather than marking it unavailable.

**`multi-heatmap marker row order` = `hormones c("INS","GCG","SST") first, then immune c("CD3e","CD4","CD8a","CD20","CD45","CD68","CD163","HLADR"), then the rest`** — apps/islet/R/mod_trajectory_server.R:1598-1603 — *ARBITRARY*  
Controls: Row ordering of the heatmap -- endocrine markers are placed at the top regardless of the question being asked.  
Checked: "Order markers: hormones first, then immune, then others" restates the code; no record of why, and measured against the panel three of the eight immune names (CD8a, CD45, HLADR) do not exist among the 58 core marker means, so CD8 and HLA_DR sort into "the rest".

**`donor_status recode map (hard-coded class list)` = `c("nd","non-diabetic","nondiabetic")->"ND"; c("aab+","aab","autoantibody")->"Aab+"; c("t1d","type1","diabetic")->"T1D"`** — apps/islet/R/mod_trajectory_server.R:162-169 — *ARBITRARY* **[claim overturned]**  
Controls: The only three disease labels normalised; any other spelling passes through unmapped and then falls outside every factor level below.  
Checked: Nothing is inherited here -- the app has a domain source for these levels (domain_group_levels(), used two hundred lines later at the legend) and this map bypasses it; no commit message mentions it and apps/islet/CLAUDE.md explicitly names hard-coded status literals as a domain leak.

**`remove_outliers default` = `reactive(FALSE)`** — apps/islet/R/mod_trajectory_server.R:17 — *ARBITRARY*  
Controls: Whether flagged points are excluded from the LOESS trend fits (they are always drawn, as grey x markers).  
Checked: Commit 04b1679 records the intent -- "Changed default outlier exclusion state to OFF for a more transparent initial display" -- a display judgement with no measurement behind it; app.R:410 overrides it with cohort$outlier_enabled.

**`default pseudotime column` = `obs2$pseudotime <- obs2$dpt_pseudotime (else NA_real_)`** — apps/islet/R/mod_trajectory_server.R:171-176 — *INHERITED_DEFAULT*  
Controls: The generic pseudotime used before any mode toggle.  
Checked: The column name is whatever the producer writes, so the app inherits it -- but the comment's "(default to PAGA dpt_pseudotime)" is wrong for the deployed file: trajectory_islet.py calls sc.tl.diffmap + sc.tl.dpt and never runs PAGA (no paga reference in that file at all).

**`outlier_threshold default` = `reactive(3.0)`** — apps/islet/R/mod_trajectory_server.R:18 — *INHERITED_DEFAULT* **[claim overturned]**  
Controls: |z| above which an islet's feature value is flagged and (when the filter is on) dropped from the trend fits.  
Checked: Not a library default -- apps/R/ieCore is first-party code in this repo; the 3 mirrors ieCore's own fallback (utils_stats.R:36) and the cohort panel's fallback (mod_cohort.R:155), both set by commit 04b1679 ("configurable z-threshold slider (default 3.0, range 0.5-10)"), and the underlying 3-sigma is convention, not measurement.

**`pseudotime mode fallback chain` = `combined if requested and present, else dpt_pseudotime, else dpt_pseudotime_combined`** — apps/islet/R/mod_trajectory_server.R:193-200 — *ARBITRARY*  
Controls: Silently substitutes the other pseudotime when the requested one is missing (a caption notes it only for the combined case).  
Checked: The comment states the behaviour without justifying it and no commit explains the ordering; the third branch is entirely silent -- asking for core on a file with only a combined column swaps in the peri-inclusive ordering with no caption at all (the note at 536-541 covers only the reverse case).

**`peri marker columns hidden from the feature selector` = `core <- nm[!grepl("^marker_mean_peri_", nm)]`** — apps/islet/R/mod_trajectory_server.R:223 — *MEASURED_HERE*  
Controls: Only core-region marker means can be put on the y-axis; the 20 um expansion means exist in the file but are unreachable here.  
Checked: I counted it in the deployed file: uns['marker_mean_names'] holds exactly 116 entries, 58 core and 58 peri, so the "116-box selector" in the comment at 1498-1499 is a real, checkable number on this dataset and the halving is what the line does.

**`hormone_markers category list` = `c("INS", "GCG", "SST")`** — apps/islet/R/mod_trajectory_server.R:238 (reused at 259-260) — *ARBITRARY*  
Controls: Grouping under "Hormone Markers"/"Hormone markers" in the feature dropdown, and which marker means head the list.  
Checked: Only "Organize features into categories" justifies it and it predates the module split (230af4d); measured, the line-238 intersect is empty on this data because var_names are the eight phenotype classes (Alpha_cell, Beta_cell, ...), so only the line-259 reuse against marker means is live -- and all three do exist among the 58 core marker means.

**`immune_markers category list` = `c("CD3e","CD4","CD8a","CD68","CD163","CD20","CD45","HLADR")`** — apps/islet/R/mod_trajectory_server.R:239 — *ARBITRARY*  
Controls: Grouping of the immune entries in the feature dropdown; a marker not listed falls into "Other Features".  
Checked: No comment or commit justifies the membership, and three of the eight names do not exist in this 60-plex at all -- the panel carries CD8, HLA_DR and no CD45 (checked against the 58 core marker-mean names) -- while the intersect target (var_names) is the eight phenotype classes, so this group never populates on the deployed file.

**`vascular_markers category list` = `c("CD31", "CD34", "SMA", "ColIV")`** — apps/islet/R/mod_trajectory_server.R:240 — *ARBITRARY*  
Controls: Grouping of vascular entries in the feature dropdown.  
Checked: Unrecorded (first appears in 230af4d); measured against the deployed panel, "ColIV" does not exist -- the marker is named Collagen_IV -- and like the other category lists it is intersected against var_names, which holds the eight phenotype classes.

**`neural_markers category list` = `c("B3TUBB", "GAP43", "PGP9.5")`** — apps/islet/R/mod_trajectory_server.R:241 — *ARBITRARY*  
Controls: Grouping of neural entries in the feature dropdown.  
Checked: Unrecorded, and two of the three (GAP43, PGP9.5) are absent from this panel's 58 core marker means -- only B3TUBB exists.

**`default y-axis feature` = `INS marker mean if present, else first marker mean, else first composition feature`** — apps/islet/R/mod_trajectory_server.R:270-272 — *ARBITRARY*  
Controls: The feature the tab opens on -- insulin by default, so the first thing every reader sees is INS versus pseudotime.  
Checked: The comment ("INS if the markers are there, else whatever leads") states the rule without a reason and no commit explains it; INS is present in the deployed marker means (median 1,693, max 29,608), so the tab does open on insulin.

**`min-cells source switch (cohort$spec()$min_cells vs input$min_cells)` = `if (!is.null(cohort)) traj_min_cells_value(cohort$spec()$min_cells) else traj_min_cells_value(input$min_cells)`** — apps/islet/R/mod_trajectory_server.R:28-31 — *INHERITED_DEFAULT* **[claim overturned]**  
Controls: Whether the tab honours the shared cohort panel's threshold or its own (now removed) numericInput; the threshold is compared against total_cells, which is core-only on this tab.  
Checked: The comment states a code-shape fact, not a measurement (no number, no cohort, no procedure); the value is inherited from mod_cohort.R:34/151 (numericInput value = 1) via commit c1a788e, and the core-only claim checks out -- the deployed obs carries total_cells and a separate total_cells_peri, and this tab never calls cohort$apply(cells_expr=).

**`pt_mode fallback / validated default` = `input$pt_mode %||% "core"; anything not in c("core","combined") -> "core"`** — apps/islet/R/mod_trajectory_server.R:33-39 — *ARBITRARY*  
Controls: Whether the x-axis reads dpt_pseudotime (core features) or dpt_pseudotime_combined (core + peri).  
Checked: The only comment describes the two columns and gives no basis for preferring core; commit 6b8c346 added the toggle with no rationale for the default, and both columns exist in the deployed file (4,371 finite values each).

**`expression source for a selected feature` = `marker -> adata$obsm[["X_marker_mean"]][, idx]; composition -> adata$X[, idx]`** — apps/islet/R/mod_trajectory_server.R:330-337 — *INHERITED_DEFAULT* **[claim overturned]**  
Controls: Which matrix the y value is read from, selected by whether the name is in uns['marker_mean_names']; a name colliding across the two vocabularies resolves to the marker.  
Checked: The comment's content is verifiably true (the deployed file carries obsm['X_marker_mean'] of shape 7302x116 with uns['marker_mean_names'], and X is the eight phenotype classes Alpha_cell...Neural) but that is the producer's file layout being followed, not a measured choice of value -- build_islet_marker_means.py decides it, this code inherits it.

**`total_cells NA fill` = `result_df$total_cells[is.na(...)] <- 1L`** — apps/islet/R/mod_trajectory_server.R:374-376 — *ARBITRARY*  
Controls: Islets with a missing cell count are given a count of 1, which becomes their LOESS weight (log1p(1)=0.69) and their point size (sqrt(1)=1).  
Checked: "Guard: ... (safe default weight)" is an assertion, not a measurement, and it substitutes a value for a missing measurement against the repo's own rule; measured on the deployed file the column has zero NAs (min 1, median 14), so the substitution is a no-op today and would only bite on a rebuild that leaves the count unavailable.

**`trajectory H5AD filename (claimed insulin-rooted)` = `"adata_ins_root.h5ad" (and domain_data_file("adata_ins_root"))`** — apps/islet/R/mod_trajectory_server.R:44-45 — *INHERITED_DEFAULT*  
Controls: The only file this tab will ever plot; the name is fixed, so a differently-named ordering cannot be displayed without renaming or re-pathing.  
Checked: The name is set by the producers -- scripts/reaggregate_islets.py:40 (legacy, root = medoid of top-K ND by INS) and core/phenocycler/scripts/typing/trajectory_islet.py:424 -- but the extractor's gloss is wrong for the deployed file: its uns says root='6476_Islet_551', root chosen as the medoid of the top-50 ND islets by endocrine_fraction (pick_root, line 130) with batch_correction='none', so the name outlived the insulin rooting it describes.

**`trajectory H5AD fallback search path` = `file.path("..","..","scripts","adata_ins_root.h5ad"); project_root/data/app_data/adata_ins_root.h5ad; project_root/scripts/adata_ins_root.h5ad`** — apps/islet/R/mod_trajectory_server.R:46-48 — *INHERITED_DEFAULT*  
Controls: Which file is loaded when the domain root has none -- these candidates escape the active domain's data root.  
Checked: Legacy pre-monorepo paths carried through commit b64ffd0 with no explanation, and the escape is not hypothetical: project_root resolves to the repo root (00_globals.R:146) and /home/smith6jt/IO60panc2nd/Islet-Explorer-Senior/data/app_data/adata_ins_root.h5ad is a real 13 MB CODEX-domain file dated Apr 25, so a missing or renamed run-20 file silently serves the other cohort under run-20 labels.

**`cohort scope applied by donor id` = `result_df filtered to cohort$donor_ids() when non-empty; skipped entirely when the list is empty`** — apps/islet/R/mod_trajectory_server.R:462-468 — *ARBITRARY* **[claim overturned]**  
Controls: Which donors reach the scatter; an empty donor list is treated as "no filter" rather than "no donors", and traj_multi_feature_data() never applies this filter at all.  
Checked: The comment records a history (the tab had no donor filter until c1a788e) rather than a measurement, and the "empty means unfiltered" semantics has no record; the extractor's second point is confirmed by reading -- traj_multi_feature_data (1510-1545) applies only the min-cells filter, so the heatmap keeps plotting donors the scatter has dropped.

**`pseudotime rank normalization` = `result_df$pt <- rank(result_df$pt, ties.method="average") / nrow(result_df) (raw kept as pt_raw)`** — apps/islet/R/mod_trajectory_server.R:478-482 — *ARBITRARY* **[claim overturned]**  
Controls: Replaces raw DPT with a within-view rank percentile, so every islet's x position depends on which other islets survived the cohort and min-cells filters.  
Checked: I tested the stated justification against the deployed file and it fails: only 13.1% of the 4,371 finite DPT values lie in 0.25-0.70 (quantiles 5%=0.047, 50%=0.127, 95%=0.560), not "95%" -- the numbers describe the legacy 5,214-islet CODEX file, so the transform is applied without a measurement that holds on this cohort.

**`ADATA_PATH environment override` = `Sys.getenv("ADATA_PATH", unset = "") (last candidate)`** — apps/islet/R/mod_trajectory_server.R:49 — *ARBITRARY* **[claim overturned]**  
Controls: Lets an operator point the tab at any H5AD, including one from the other domain.  
Checked: Nothing inherits this -- no config key, no other module reads ADATA_PATH, and no commit message explains it; the only trace is the user-facing hint at line 493, so it is an unrecorded escape hatch that bypasses domain_data_root() entirely.

**`validity caption gate and thresholds display` = `reads uns[["trajectory_validation"]][["core"]]; badge #228833 if hard_pass else #CC6633 "check gates"; prints chosen_theta, INS ρ, T1D>ND p, within-status batch η²`** — apps/islet/R/mod_trajectory_server.R:516-534 — *INHERITED_DEFAULT* **[claim overturned]**  
Controls: The tab's headline verdict on the ordering -- sourced from gates that consult INS and a T1D>ND direction, computed upstream and only displayed here.  
Checked: I opened the deployed H5AD: its uns['trajectory_validation'] is FLAT (gate_root_is_ND, gate_stable_under_reroot, hard_pass=True, mwu_T1D_gt_ND_p, spearman_pseudotime_beta_proportion, ...) with no 'core' member and no chosen_theta/spearman_INS/within_status_eta2_*, so this caption never renders on the run-20 domain -- the keys it reads are scripts/senior/compute_trajectory_senior.py's (the Harmony θ-sweep whose root is the top-K ND medoid by INS and whose hard gate requires spearman(pt,INS) < -0.3, GCG and an ND<Aab+<T1D ordering), which is a different producer from the one that wrote this file.

**`single-donor support caption threshold` = `s$donors == 1`** — apps/islet/R/mod_trajectory_server.R:554 (and the legend copy at 1180) — *ARBITRARY* **[claim overturned]**  
Controls: A status carried by exactly one donor gets an italic warning and is described as "excluded from the ordering statistic"; two donors get no warning.  
Checked: The comment documents a real fix (replacing hard-coded "Aab+ is a single donor (6521)" literals with status_support() computed from obs) but the threshold itself is unmeasured and disagrees with the producer: trajectory_islet.py sets DONOR_ETA2_MIN_DONORS = 3, so a two-donor status is also excluded from that statistic yet triggers no caption -- and on the deployed file (8 ND, 2 T1D, 0 Aab+) the warning never fires while the T1D η² of 0.46 rests on one pairwise contrast.

**`donor_status factor levels (hard-coded class list)` = `c("ND", "Aab+", "T1D")`** — apps/islet/R/mod_trajectory_server.R:592, 671, 688, 787, 793, 1281, 1384 — *ARBITRARY* **[claim overturned]**  
Controls: Legend order, colour mapping, trend-line loop and the heatmap's status set; a dataset with any other status label drops out of the factor and its points lose their colour.  
Checked: Nothing is inherited -- the app has domain_group_levels() and uses it in this very file at line 1170, so these seven literals bypass the DOMAIN, exactly the leak apps/islet/CLAUDE.md warns about ("choices = c(\"ND\", \"Aab+\", \"T1D\") ... every one a fact about the CODEX cohort"); measured, the deployed run-20 file contains zero Aab+ islets.

**`flag_outliers(group_col=)` = `group_col = "donor_status"`** — apps/islet/R/mod_trajectory_server.R:612 — *INHERITED_DEFAULT* **[claim overturned]**  
Controls: Outliers are judged within a disease-status pool of donors, not within a donor; with one donor loaded it collapses to a single group.  
Checked: Same first-party-helper point (utils_stats.R:27), not a library default; measured on the deployed run-20 file this pools 8 ND donors into one group and 2 T1D donors into another (Aab+ has zero islets), and groups with fewer than 3 finite values are skipped outright (utils_stats.R:86).

**`flag_outliers(method=)` = `method = "zscore"`** — apps/islet/R/mod_trajectory_server.R:613 — *INHERITED_DEFAULT* **[claim overturned]**  
Controls: Chooses the robust median/MAD modified z-score rule over the IQR rule for flagging.  
Checked: ieCore is in-repo, not an installed library, so this is the repo's own helper default (utils_stats.R:28) restated at the call site; I read the implementation and the in-tab wording is right -- "zscore" is Iglewicz-Hoaglin 0.6745*(x-median)/MAD with a classical mean/SD fallback when MAD=0 (utils_stats.R:97-108), while the stale file-header comment at utils_stats.R:6-7 still describes mean/SD.

**`jitter_amount_x` = `diff(range(df$pt)) * 0.012`** — apps/islet/R/mod_trajectory_server.R:652,654,657 — *ARBITRARY*  
Controls: Horizontal displacement applied to every plotted islet -- the displayed pseudotime is not the computed pseudotime; at the rank scale this is +/-1.2% of the axis.  
Checked: The comment explains only why the kept-rows range is used as the scale, never where 0.012 came from, and git -S finds it only in the pre-split monolith.

**`jitter_amount_y` = `diff(range(df$value)) * 0.02`** — apps/islet/R/mod_trajectory_server.R:653,655,658 — *ARBITRARY*  
Controls: Vertical displacement of every plotted value, +/-2% of the value range, including for marker intensities.  
Checked: No comment, no commit; on marker means this is +/-2% of a range that reaches 29,608, i.e. hundreds of intensity units of displacement on a value the axis label calls a mean intensity.

**`jitter clamp on pseudotime` = `pmax(0, pmin(1, pt + runif(...)))`** — apps/islet/R/mod_trajectory_server.R:654,657 — *ARBITRARY*  
Controls: Piles jittered points at exactly 0 and 1, so the two extremes of the ordering are visually denser than they are.  
Checked: Unexplained anywhere; the clamp is only coherent because of the rank normalization above it (raw DPT is not on 0-1 for this file, whose 95th percentile is 0.56), and it is applied to the ghost layer too.

**`LOESS weight transform` = `df$loess_weight <- log1p(total_cells)`** — apps/islet/R/mod_trajectory_server.R:661-662 — *ARBITRARY* **[claim overturned]**  
Controls: How much a large islet pulls the trend line relative to a small one.  
Checked: The cited measurement is real but belongs to another cohort and no longer holds: commit 4a1b5c9 (Mar 2026, pre-run-20) introduced it citing "60% have <=10 cells", and the comment's "median=9, max=1902" does not match the deployed file (median 14 / max 1928 over all islets; median 34 and only 4.0% at <=10 cells among the 4,371 actually plotted), and no procedure compares log1p against any alternative.

**`traj_min_cells_value() floor / default` = `v <- as.numeric(x %||% 1); if (length(v)==1 && is.finite(v) && v >= 1) v else 1`** — apps/islet/R/mod_trajectory_server.R:7-10 — *ARBITRARY*  
Controls: Islets with fewer than this many CORE cells are dropped before the scatter, the rank normalization and both heatmaps; the default 1 means no filtering.  
Checked: git log ea33878 ("Stop a cleared Min cells/islet box from blanking the Trajectory tab") shows the helper was added purely for NULL/NA coercion, and 1 is simply the off position -- no measurement anywhere; the real cell floor on this data is upstream (COMPOSITION_FLOOR_CELLS=10 in trajectory_islet.py, which already leaves 2,931 of 7,302 islets with no pseudotime).

**`geom_smooth(method=)` = `method = "loess"`** — apps/islet/R/mod_trajectory_server.R:781, 799, 806 — *ARBITRARY*  
Controls: The trend model drawn over the points, at all three call sites.  
Checked: No record anywhere, and it is a non-default at this size: ggplot2's method = "auto" picks loess only below 1,000 rows and gam above it, while the deployed file plots 4,371 islets.

**`geom_smooth(se=)` = `se = FALSE (all three call sites)`** — apps/islet/R/mod_trajectory_server.R:781, 799, 806 — *ARBITRARY*  
Controls: Suppresses the confidence band, so a curve fitted from three islets looks as firm as one fitted from three thousand.  
Checked: Confirmed non-default -- formals(ggplot2::geom_smooth)$se is TRUE in the installed ggplot2 -- and no comment or commit explains the override.

**`LOESS span` = `span = 0.75`** — apps/islet/R/mod_trajectory_server.R:782, 800, 807 — *LIBRARY_DEFAULT*  
Controls: Smoothing bandwidth of every trend line -- how much local structure along pseudotime survives into the curve.  
Checked: Confirmed against the installed R 4.5.1: args(stats::loess) shows span = 0.75, so all three call sites merely write out the package default.

**`y-axis metric label switch` = `marker -> "— mean aligned intensity per islet (core)"; else "— proportion of islet-core cells (0–1)"; x label "Pseudotime (rank-normalized)"`** — apps/islet/R/mod_trajectory_server.R:819-838 — *MEASURED_HERE*  
Controls: What the axis claims the numbers are; the branch is keyed purely on the "marker_mean_" prefix.  
Checked: The load-bearing number checks out on the deployed file: obsm['X_marker_mean'] runs 0 to 29,608 (INS median 1,693), so "run to five figures" is a real measurement of this dataset and the old 0-1 label would indeed have been false for every marker.


### `apps/islet/R/mod_trajectory_ui.R`

**`traj_show_trend default` = `choices none / overall / by_donor, selected = "by_donor"`** — apps/islet/R/mod_trajectory_ui.R:20-24 (server fallback identical, line 645) — *ARBITRARY*  
Controls: Draws one LOESS curve per disease status by default -- a cohort-pooled summary laid over the points.  
Checked: No comment, and git -S finds it only in the pre-split monolith (230af4d); the label "By Donor Status" makes it a disease-pooled default with nothing recorded behind the choice.

**`pt_mode selectInput default` = `choices = c("Core only"="core", "Core + Peri"="combined"), selected = "core"`** — apps/islet/R/mod_trajectory_ui.R:25-28 — *ARBITRARY*  
Controls: The pseudotime column the tab opens on.  
Checked: Introduced by commit 6b8c346 whose message says nothing about which mode should lead, and no comment justifies it.

**`traj_color_by choices and default` = `donor_status / donor_id / leiden_0.3 / 0.5 / 0.8 / 1.0, selected = "donor_status"`** — apps/islet/R/mod_trajectory_ui.R:31-38 (server fallback identical, line 643) — *ARBITRARY*  
Controls: What the scatter is coloured by on load; the default frames every islet by disease status.  
Checked: No comment and no commit message addresses the default; it predates the module split (230af4d).

**`multi_heatmap_nbins` = `min = 10, max = 40, value = 20, step = 5`** — apps/islet/R/mod_trajectory_ui.R:74 (server fallback 20 at line 1554) — *ARBITRARY*  
Controls: Number of pseudotime bins the multi-marker heatmap averages within -- and therefore how many islets back each tile.  
Checked: git -S traces it to 184f523 ("Add phenotype composition, demographic filters, multi-feature heatmap, embedded segmentation"), which records nothing about the bin count.


### `core/phenocycler/scripts/typing/trajectory_islet.py`

**`sc.pp.neighbors relied-on defaults (method, knn, random_state)` = `method="umap", knn=True, random_state=0 (not passed)`** — core/phenocycler/scripts/typing/trajectory_islet.py:123 — *LIBRARY_DEFAULT*  
Controls: The UMAP connectivity kernel and its fixed seed produce the transition matrix DPT integrates.  
Checked: Confirmed against the installed scanpy 1.12.3 at .venv/lib/python3.12/site-packages/scanpy/neighbors/__init__.py:83-95 (n_neighbors=15, knn=True, method="umap", random_state=0) - note metric there defaults to None, not "euclidean", so passing METRIC explicitly is what pins it.

**`sc.tl.diffmap(n_comps=...)` = `10`** — core/phenocycler/scripts/typing/trajectory_islet.py:125 — *INHERITED_DEFAULT*  
Controls: Number of diffusion components computed, capping the dimensionality DPT then integrates over.  
Checked: scanpy 1.12.3's own default is n_comps=15 (tools/_diffmap.py:17), so 10 is a deliberate override that matches the legacy lane verbatim (compute_trajectory_senior.py:141 and :336) with no comment or commit text in either file.

**`sc.tl.dpt() relied-on defaults` = `n_dcs=10, n_branchings=0, min_group_size=0.01 (none passed)`** — core/phenocycler/scripts/typing/trajectory_islet.py:126 — *LIBRARY_DEFAULT*  
Controls: n_dcs consumes exactly the 10 components diffmap computed; n_branchings=0 yields a single unbranched ordering.  
Checked: Confirmed in .venv scanpy/tools/_dpt.py:43-46 - n_dcs=10, n_branchings=0, min_group_size=0.01 - so the apparent agreement between diffmap's 10 and dpt's n_dcs is a coincidence of two independently set numbers, and lowering diffmap's n_comps alone would silently change what dpt integrates.

**`--donors (default None -> cohort.CLEAN)` = `None, resolving to the 10 CLEAN donors`** — core/phenocycler/scripts/typing/trajectory_islet.py:155 — *INHERITED_DEFAULT*  
Controls: Which donors are pooled into one feature matrix and one graph - the default is what makes this a cohort-pooled ordering.  
Checked: cohort.resolve(None) returns the 10 donors I confirmed against the data (7,780 islets, 7,302 with a core endocrine cell), and core/phenocycler/scripts/typing/cohort.py documents the six-setting missed-islet sweep behind that set - so the value is inherited here from a measurement made there, not made here.

**`--include-marginal` = `False (store_true)`** — core/phenocycler/scripts/typing/trajectory_islet.py:156 — *INHERITED_DEFAULT*  
Controls: Whether 6521, 6551 and 6380 join the cohort, including the only Aab+ donor outside batch 2.  
Checked: The store_true default mirrors cohort.resolve(include_marginal=False); cohort.py's MARGINAL dict prices each donor (6380 alone at 9.0% missed islets) and its COHORT_LIMITS calls including 6521 "a maintainer decision, not a default".

**`classes (derived from count_* columns)` = `sorted(c[len("count_"):] for c in frame.columns if c.startswith("count_"))`** — core/phenocycler/scripts/typing/trajectory_islet.py:161 — *INHERITED_DEFAULT* **[claim overturned]**  
Controls: Sets both the ALR dimensionality and adata.var; an upstream taxonomy change silently changes the feature space.  
Checked: Not UNKNOWN but inherited from upstream data: I read the CSV and it yields exactly eight classes today (Alpha_cell, Beta_cell, Blood_Vessel, Delta_cell, Exocrine, Immune, Lymphatic, Neural), which are precisely the DISPLAY values written by build_cells.py:41-47, so the feature space is pinned by that map and by nothing in this file.

**`defined = endocrine_core > 0` = `> 0`** — core/phenocycler/scripts/typing/trajectory_islet.py:167 — *ARBITRARY* **[claim overturned]**  
Controls: Islets with no core endocrine cell are dropped from the matrix entirely before the ALR is built.  
Checked: Downgraded from MEASURED_HERE: no measurement produced it and the stated rationale does not hold - line 113 adds ALR_PSEUDOCOUNT to the denominator too, so a zero-endocrine islet yields the perfectly finite log((n+0.5)/0.5), making this a modelling choice rather than an undefined ratio; it drops 478 of 7,780 islets, 351 of them from donor 6539 alone.

**`ALR denominator = endocrine_core for both the core and peri blocks` = `alr_matrix(frame, classes, "", endocrine_core) and alr_matrix(frame, classes, "peri_", endocrine_core)`** — core/phenocycler/scripts/typing/trajectory_islet.py:179-180 — *ARBITRARY* **[claim overturned]**  
Controls: Both region blocks are normalised by the islet's core endocrine count rather than by their own region total, so the peri axis reads per endocrine cell in the islet.  
Checked: Determinable rather than UNKNOWN: the docstring gives a design rationale and no measurement, `git log -S` shows the construction never had another form, and no ancestor supplies it - the legacy lane normalises a scVI latent and weights its second block explicitly (ALPHA_STRUCT = 0.15), so this is a declared choice made in 6b8c346.

**`adata.X scaling = counts / total_cells_core` = `np.divide(counts, core_total[:, None], out=np.full_like(counts, np.nan), where=core_total[:, None] > 0)`** — core/phenocycler/scripts/typing/trajectory_islet.py:187-190 — *ARBITRARY* **[claim overturned]**  
Controls: What the Trajectory tab plots against pseudotime and what the reported Beta_cell correlation at line 300 is computed on; it is not what the ordering uses (that is the obsm ALR blocks).  
Checked: The cited observation is real - I measured max per-class count 1,741 and max total_cells_core 1,928 in this cohort's CSV, so "counts run to ~1,500" checks out - but it motivated a label-correctness fix rather than selecting a value, and no alternative scaling was compared, which is below the bar for MEASURED_HERE.

**`primary ordering representation` = `"X_composition_alr" (core only) for dpt_pseudotime; "X_composition_alr_combined" (core hstack peri) for dpt_pseudotime_combined`** — core/phenocycler/scripts/typing/trajectory_islet.py:236-237 — *UNKNOWN*  
Controls: Which feature block the deployed pseudotime is built on, and that the combined block is an unweighted concatenation.  
Checked: Nothing records it - no comment, no commit text, and no ancestor to inherit from, since the legacy lane's combined mode weights its second block explicitly (ALPHA_STRUCT = 0.15 at compute_trajectory_senior.py:59) while line 218 here is a plain hstack that lets the peri block double the dimensionality at equal weight.

**`hardcoded reporting class "Beta_cell"` = `"Beta_cell"`** — core/phenocycler/scripts/typing/trajectory_islet.py:300 — *INHERITED_DEFAULT*  
Controls: Which class's per-islet proportion is correlated against pseudotime for spearman_pseudotime_beta_proportion; NaN if the taxonomy renames it.  
Checked: Nothing in this file justifies the choice; it is the composition analogue of the legacy INS check named in the docstring, the name resolves today because build_cells.py's DISPLAY map emits it, and the deployed run recorded spearman_pseudotime_beta_proportion = -0.099.

**`sc.tl.umap() relied-on defaults` = `min_dist=0.5, spread=1.0, n_components=2, init_pos="spectral", random_state=0 (none passed)`** — core/phenocycler/scripts/typing/trajectory_islet.py:338 — *LIBRARY_DEFAULT*  
Controls: The islet-level UMAP written to obsm and to islet_leiden.parquet as leiden_umap_1/2; display only, no effect on the ordering.  
Checked: Confirmed in .venv scanpy/tools/_umap.py:44-52 (min_dist=0.5, spread=1.0, n_components=2, init_pos="spectral", random_state=0); the Junior lane's different values (reaggregate_islets.py:162, init_pos='paga', min_dist=0.3, spread=1.5) and the legacy lane's spread=2.0 (compute_trajectory_senior.py:426) were both dropped without comment.

**`sc.tl.leiden(flavor, n_iterations, directed)` = `flavor="igraph", n_iterations=2, directed=False`** — core/phenocycler/scripts/typing/trajectory_islet.py:343-344 — *LIBRARY_DEFAULT* **[claim overturned]**  
Controls: Which Leiden backend and how many refinement iterations produce the islet-level clusters written for the Spatial tab; does not touch the pseudotime.  
Checked: Not UNKNOWN: scanpy 1.12.3 raises a FutureWarning naming exactly these three arguments (tools/_leiden.py:236-238, "To achieve the future defaults please pass: flavor='igraph' and n_iterations=2. directed must also be False"), so they are the library's announced future defaults copied verbatim, not a local tuning choice.

**`ALR_PSEUDOCOUNT` = `0.5`** — core/phenocycler/scripts/typing/trajectory_islet.py:55 — *ARBITRARY*  
Controls: Additive pseudocount on both the class count and the core-endocrine denominator at line 112-113, so it sets how hard zero- and low-count classes are pulled toward the denominator and therefore the geometry of every islet's ALR vector.  
Checked: `git log -S ALR_PSEUDOCOUNT --all` returns only 6b8c346, the single commit that created the file, whose message never mentions the pseudocount, and no sweep over it exists anywhere in the repo - the comment declares the value is recorded, not that 0.5 was arrived at.

**`METRIC` = `"euclidean"`** — core/phenocycler/scripts/typing/trajectory_islet.py:61 — *MEASURED_HERE*  
Controls: Distance sc.pp.neighbors uses to build the kNN graph the diffusion map and DPT run on (line 123).  
Checked: I re-ran the cited comparison on this cohort: cosine on the ALR at floor 0 gives rho(n_neighbors 30 vs 15) = 0.309, matching the comment's 0.31 to the digit, and the legacy lane it departs from does use metric="cosine" (compute_trajectory_senior.py:138,334) - only the "0.97 median pair" figure is loose, reproducing as 0.935 at floor 0 and 0.946 at floor 10.

**`N_NEIGHBORS` = `30`** — core/phenocycler/scripts/typing/trajectory_islet.py:62 — *INHERITED_DEFAULT*  
Controls: kNN size of the primary graph; sets the scale of local structure the diffusion map and pseudotime see.  
Checked: Identical uncommented constant at scripts/senior/compute_trajectory_senior.py:53 (and n_neighbors=30 at reaggregate_islets.py:148), with nothing in 6b8c346 justifying it - and it does not transfer per donor, where 6539 has only 67 islets above the floor, so k=30 is 45% of its data rather than 0.7% of the pooled 4,371.

**`N_NEIGHBORS_ALT` = `15`** — core/phenocycler/scripts/typing/trajectory_islet.py:63 — *ARBITRARY*  
Controls: Neighbourhood size of the stability re-run; the Spearman rho against the primary run is what stable_under_neighbourhood_size gates on (line 291), and the composition-floor sweep is defined as rho(30 vs 15).  
Checked: No measurement, comment or commit text chooses it; 15 is both scanpy 1.12.3's own n_neighbors default (neighbors/__init__.py:83) and the secondary size used in the two ancestor lanes (compute_trajectory_senior.py:353, reaggregate_islets.py:153,509), but the file cites neither.

**`COMPOSITION_FLOOR_CELLS` = `10`** — core/phenocycler/scripts/typing/trajectory_islet.py:80 — *MEASURED_HERE*  
Controls: Minimum core cells for an islet to enter the neighbourhood graph and receive a pseudotime (line 223); islets below it stay in the H5AD with an Unavailable reason.  
Checked: A real sweep - every islet count reproduces exactly off islet_celltype_proportions.csv (7,302 / 5,707 / 4,371 / 3,101 / 2,399 / 1,687, median 14 core cells) - but its rho column is a COSINE sweep, not the deployed metric: at floor 10 cosine gives exactly the table's 0.931 while the deployed euclidean run recorded rho_vs_smaller_neighbourhood = 0.9777 in the H5AD's uns, and under euclidean rho is >= 0.968 at every floor from 5 up, so nothing in the shipped configuration selects 10; separately, dropping the floor to 0 under euclidean leaves 380 islets non-finite and rho 0.829, failing both all_islets_finite and the 0.90 gate.

**`ENDOCRINE` = `("Beta_cell", "Alpha_cell", "Delta_cell")`** — core/phenocycler/scripts/typing/trajectory_islet.py:81 — *INHERITED_DEFAULT* **[claim overturned]**  
Controls: The classes summed into the ALR denominator (line 164-166), the endocrine_core > 0 keep rule, and endocrine_fraction, which is the root ranking key.  
Checked: Not UNKNOWN: the same three-member set is declared upstream in this very lane at core/phenocycler/scripts/typing/build_cells.py:47 (ENDOCRINE = ("Alpha", "Beta", "Delta"), display-renamed by the DISPLAY map two lines above), and the cohort CSV carries exactly eight classes with no fourth endocrine candidate, so the membership is inherited and constrained rather than unexplained.

**`STABILITY_MIN_RHO` = `0.90`** — core/phenocycler/scripts/typing/trajectory_islet.py:83 — *ARBITRARY*  
Controls: Spearman floor for both stability gates (lines 291-292); below it the build is written .rejected.h5ad and exits 1.  
Checked: A round number with a stated purpose and no cited or recoverable measurement, and it is load-bearing at close range - the deployed H5AD records rho_vs_alternate_root_median = 0.9243 (passing by 0.024) with rho_vs_alternate_root_min = 0.5295.

**`LEIDEN_RESOLUTIONS` = `(0.3, 0.5, 0.8, 1.0)`** — core/phenocycler/scripts/typing/trajectory_islet.py:86 — *INHERITED_DEFAULT*  
Controls: The four islet-level Leiden clusterings written to the H5AD and islet_leiden.parquet for the Spatial tab; does not enter the pseudotime.  
Checked: The identical four values are in both ancestors (compute_trajectory_senior.py:355, reaggregate_islets.py:519), though the stated reason overstates the constraint: mod_spatial_server.R:52 greps "^leiden_[0-9]" and builds its selector from whatever columns exist, so the app does not require these four.


### `core/phenocycler/scripts/typing/trajectory_singlecell.R`

**`set.seed` = `20260914`** — core/phenocycler/scripts/typing/trajectory_singlecell.R:40 — *ARBITRARY*  
Controls: Which cells the subsample draws; the only randomness in the script, so it fixes the matrix the diffusion map is built on.  
Checked: `git log` shows the file's only commit, 6b8c346, is dated Mon Sep 14 14:25:07 2026, so the seed is that day's date written as an integer; the comment states only that the subsample is reproducible, never why this value.

**`N_FIT (positional CLI arg 1)` = `100000L (default). The artifact on disk was built at 50000.`** — core/phenocycler/scripts/typing/trajectory_singlecell.R:42 — *ARBITRARY*  
Controls: How many cells are ordered; everything outside the subsample carries no pseudotime, and the destiny cost is roughly quadratic in it.  
Checked: Nothing states why 100,000 — the '25,000 cells take 94 s' timing appears only in this comment with no recorded run — and the in-file default has in fact never been exercised: singlecell_trajectory.json records n_ordered 49999 of n_population 800562, matching README.md:70's documented invocation with 50000, so the file's default and the lane's documented value disagree.

**`N_DCS` = `20L`** — core/phenocycler/scripts/typing/trajectory_singlecell.R:43 (used at :78 as n_eigs and at :81 as the slingshot input dimension) — *LIBRARY_DEFAULT*  
Controls: Eigenvectors computed by destiny (n_eigs) and the dimensionality of the space slingshot fits its curve in (rd = eigenvectors[, 1:20]).  
Checked: `args(destiny::DiffusionMap)` in the installed destiny 3.22.0 gives `n_eigs = min(20L, n - 2L)`, so 20 is the library default at this fit size — but it is written as a literal (it would not follow a change in destiny) and nothing states why all 20 eigenvectors are handed to slingshot, which is a choice destiny does not make.

**`proportional per-donor quota rule (round(N_FIT * table(donor_id) / nrow))` = `quota <- round(as.numeric(N_FIT) * table(cells$donor_id) / nrow(cells)); sample(rows, min(quota[[d]], length(rows)))`** — core/phenocycler/scripts/typing/trajectory_singlecell.R:63-67 — *ARBITRARY* **[claim overturned]**  
Controls: How the fit budget is split across donors. On the artifact this gave 1,576 cells to 6539 and 7,679 to 115, a 4.9x spread. In a single-donor run it collapses to a plain random subsample of size min(N_FIT, n).  
Checked: The 'coarse embedding' it is attributed to does not exist in the checkout — `grep -rni quota` over core/phenocycler/scripts finds proportional per-donor quotas only in this file and an unrelated test fixture (broadtype/tests/conftest.py), and 'coarse embedding' appears in no file outside commit 6b8c346, so there is no established lane precedent to inherit from.

**`DiffusionMap k (nearest neighbours), relied-on default` = `not passed -> find_dm_k(n - 1) = 100`** — core/phenocycler/scripts/typing/trajectory_singlecell.R:78 — *LIBRARY_DEFAULT*  
Controls: The neighbourhood size of the transition kernel: how many nearest neighbours each cell is connected to before the eigendecomposition.  
Checked: I read destiny::find_dm_k from the installed 3.22.0: `min_k = 100L, small = 1000L, big = 10000L`, returning min_k whenever n >= big, and find_dm_k(49998) evaluates to 100 for the 49,999-cell fit that is actually on disk.

**`DiffusionMap sigma / n_local, relied-on defaults` = `sigma = "local", n_local = seq(to = min(k, 7L), length.out = min(k, 3L)) = 5:7 at k = 100`** — core/phenocycler/scripts/typing/trajectory_singlecell.R:78 — *LIBRARY_DEFAULT*  
Controls: The kernel bandwidth per cell, taken from the distance to its 5th-7th nearest neighbours; the single largest lever on the shape of DC1.  
Checked: `args(destiny::DiffusionMap)` on the installed 3.22.0 shows `sigma = "local"` and `n_local = seq(to = min(k, 7L), length.out = min(k, 3L))`, which at k = 100 is 5, 6, 7; neither argument appears in the call at line 78 (only n_eigs and verbose are passed).

**`DiffusionMap density_norm, relied-on default` = `TRUE`** — core/phenocycler/scripts/typing/trajectory_singlecell.R:78 — *LIBRARY_DEFAULT*  
Controls: Whether the kernel is normalised for local sampling density, which decides how much dense regions dominate the diffusion components.  
Checked: `density_norm = TRUE` is in the installed destiny 3.22.0 signature and no density argument appears in the line-78 call.

**`DiffusionMap distance, relied-on default` = `"euclidean" (first of c("euclidean", "cosine", "rankcor", "l2"))`** — core/phenocycler/scripts/typing/trajectory_singlecell.R:78 — *LIBRARY_DEFAULT*  
Controls: The metric on the [.., 1]-scaled marker matrix, whose values run from -0.748 to 1.0 with a median of -0.004 once the lower clip is removed.  
Checked: The installed signature is `distance = c("euclidean", "cosine", "rankcor", "l2")`, so match.arg takes euclidean, and no distance argument appears in the call — worth flagging that the sibling islet lane (trajectory_islet.py) was explicitly moved off cosine because cosine on all-negative values collapsed similarity, per commit 6b8c346's message.

**`DiffusionMap n_pcs, relied-on default` = `NULL (no PCA before the kernel)`** — core/phenocycler/scripts/typing/trajectory_singlecell.R:78 — *LIBRARY_DEFAULT*  
Controls: Whether distances are computed on the raw 38-marker matrix or on a PCA of it; NULL means every marker enters at its scaled weight.  
Checked: `n_pcs = NULL` is in the installed destiny 3.22.0 signature and no n_pcs argument appears in the call at line 78.

**`slingshot() clusterLabels omitted` = `sl <- slingshot(rd) — no clusterLabels`** — core/phenocycler/scripts/typing/trajectory_singlecell.R:86 — *INHERITED_DEFAULT*  
Controls: With no clusters there is no MST and no branching: a single principal curve through the 20-DC cloud.  
Checked: The paper's Methods quote is in the file header at lines 8-16, and I confirmed the mechanism in the installed slingshot 2.16.0: getMethod("slingshot", c("matrix","ANY")) replaces a missing clusterLabels with rep("1", nrow(data)), and the run's manifest records n_lineages 1 with exactly one pseudotime_lineage_ column in the parquet.

**`slingshot() start.clus / end.clus omitted` = `not passed -> start.clus = NULL, end.clus = NULL`** — core/phenocycler/scripts/typing/trajectory_singlecell.R:86 — *LIBRARY_DEFAULT* **[claim overturned]**  
Controls: CORRECTED: where pseudotime zero sits is NOT inferred from lineage structure — with one cluster there is none. It is set by the sign of PC1 of the 20-DC matrix.  
Checked: getLineages' defaults are indeed start.clus = NULL / end.clus = NULL in slingshot 2.16.0, but I read getCurves' K == 1 branch: it builds the initial line as `rbind(ctr - 10*pca$sdev[1]*pca$rotation[,1], ctr, ctr + 10*pca$sdev[1]*pca$rotation[,1])` from prcomp(X), so pseudotime 0 lands at whichever end of PC1 the SVD sign convention put first — an arbitrary orientation, not an inferred root, and one that can flip between donors in a per-donor build.

**`getCurves approx_points, relied-on default` = `NULL -> 150`** — core/phenocycler/scripts/typing/trajectory_singlecell.R:86 — *LIBRARY_DEFAULT*  
Controls: The curve is represented by 150 sampled points for 49,999 cells, so pseudotime resolution and the curve's smoothness are set by this, not by the cell count.  
Checked: In the installed slingshot 2.16.0 the getCurves .local formals give `approx_points = NULL` and the body sets it to 150 when nrow(X) > 150 (else FALSE); the fit has 49,999 rows, so 150 it is.

**`getCurves smoother (and its df), relied-on defaults` = `smoother = "smooth.spline", df = 5, tol = 1e-04`** — core/phenocycler/scripts/typing/trajectory_singlecell.R:86 — *LIBRARY_DEFAULT*  
Controls: How flexible the principal curve is; df = 5 fixes how much of the 20-DC geometry the ordering can follow.  
Checked: `smoother = "smooth.spline"` is in the getCurves .local formals and the switch's smooth.spline branch hard-codes `df = 5, tol = 1e-04` inside the smootherFcn closure; it is genuinely live on the single-lineage path — the iteration loop calls smootherFcn per dimension regardless of L.

**`getCurves stretch and extend, relied-on defaults` = `stretch = 2 (live), extend = "y" (dead on this code path)`** — core/phenocycler/scripts/typing/trajectory_singlecell.R:86 — *LIBRARY_DEFAULT* **[claim overturned]**  
Controls: CORRECTED: stretch = 2 caps how far cells project past the curve's ends on each refit iteration. extend controls nothing here — the single-cluster branch never reads it.  
Checked: Both values are confirmed defaults in slingshot 2.16.0, but I traced their use: `stretch` is read inside the maxit loop (project_to_curve(X, s = s, stretch = stretch)), while `extend` is read only in the K > 1 branch — the K == 1 branch ends with `next` before those lines and hard-codes stretch = 9999 then 0 for its initial projection, so extend = "y" is unreachable for this run.

**`getCurves thresh and maxit, relied-on defaults` = `thresh = 0.001, maxit = 15`** — core/phenocycler/scripts/typing/trajectory_singlecell.R:86 — *LIBRARY_DEFAULT*  
Controls: The convergence tolerance and iteration cap of the principal-curve fit; hitting maxit returns an unconverged curve without an error.  
Checked: Confirmed in the getCurves .local formals of slingshot 2.16.0 (`thresh = 0.001, maxit = 15`), with the loop written `while (!hasConverged && it < maxit)` and no warning on exhaustion — and this script records neither the iteration count nor the convergence flag in singlecell_trajectory.json, so an unconverged fit would be indistinguishable from a converged one.


### `core/phenocycler/scripts/typing/trajectory_singlecell_features.py`

**`no lower clip on the scaled values (apply_scale's floor/span with the 0.0 clip removed)` = `scaled = (values - floor) / span, with only np.minimum(scaled, 1.0) applied — no lower bound`** — core/phenocycler/scripts/typing/trajectory_singlecell_features.py:106-107 (rationale in the comment at :101-102 and the docstring at :36-38) — *ARBITRARY* **[claim overturned]**  
Controls: Whether sub-floor values keep their negative magnitude or collapse to exactly zero; sets how much of the matrix is degenerate before the kernel is built.  
Checked: I measured the built matrix myself (800,562 cells x 38 markers = 30,421,356 values, the 10 CLEAN donors) and 57.233% of values fall below the floor, not the docstring's 55.1%; that figure appears nowhere else in the repo (`grep -rn "55.1"` hits only this docstring), names no cohort or procedure, and no measurement anywhere compares the clipped against the unclipped ordering — so the departure from subpanel_io.apply_scale rests on an unreproducible number, not a measurement.

**`upper clip at 1.0 (np.minimum(scaled, 1.0))` = `1.0`** — core/phenocycler/scripts/typing/trajectory_singlecell_features.py:107 — *INHERITED_DEFAULT*  
Controls: Caps every marker at the cohort top, so the brightest cells become indistinguishable on that marker.  
Checked: subpanel_io.py:73 is `np.clip((values - lo) / span, 0.0, 1.0)`, so the 1.0 is carried verbatim while the 0.0 is dropped; I measured 0.276% of finite values sitting at exactly 1.0 in the built matrix, so the cap is live but narrow.

**`feature_markers() source set (islet_subpanel + subpanels + hormones + bottom_tier_immune + SUBPANEL['exocrine'] + refinement own_markers)` = `38 markers: 36 read from data/development/typing/subpanel/subpanel_log.json, umap/refinement_log.json and subpanel_io.SUBPANEL['exocrine'], plus CD11b and IDO1 added in code at :78. Recorded in singlecell_features.json and singlecell_trajectory.json.`** — core/phenocycler/scripts/typing/trajectory_singlecell_features.py:64-80 — *INHERITED_DEFAULT*  
Controls: The entire feature space of the diffusion map: which markers define 'how similar are two cells'.  
Checked: I reproduced the function against the on-disk logs and got exactly the 38 markers the manifest records, 36 of them genuinely read from subpanel_log.json / refinement_log.json / subpanel_io.SUBPANEL — but the docstring's "rather than chosen here" is inexact, because two markers (CD11b, IDO1) are chosen here on line 78.

**`hard-coded marker append ["CD11b", "IDO1", "HLA_DR", "Granzyme_B"]` = `markers += ["CD11b", "IDO1", "HLA_DR", "Granzyme_B"] — net effect: adds CD11b and IDO1 (2 of the final 38); HLA_DR and Granzyme_B are already in bottom_tier_immune, so the append is redundant for those two.`** — core/phenocycler/scripts/typing/trajectory_singlecell_features.py:78 — *INHERITED_DEFAULT* **[claim overturned]**  
Controls: Adds two markers to the feature space that the classification logs do not supply (5.3% of the 38 features), weighting the embedding toward myeloid/activation signal.  
Checked: A source does exist: the four are exactly the immune channels in build_phenotype_rules.py:33-38 MARKERS that the refinement log's Immune own_markers omits — the log lists CD3e, CD68, Iba1, CD20, CD79a, MPO, CD66, CD11c, CD209, and MARKERS adds CD11b, IDO1, HLA_DR, Granzyme_B to exactly that block, so this line restores the classification's own immune marker list rather than inventing one.

**`population filter: cell_class != "Excluded" and cell_region_20um in ["core", "expansion"]` = `cells[(cells["cell_class"] != "Excluded") & cells["cell_region_20um"].isin(["core", "expansion"])]`** — core/phenocycler/scripts/typing/trajectory_singlecell_features.py:90-91 — *INHERITED_DEFAULT*  
Controls: Which cells are ordered at all: islet core plus the 20 um peri-islet expansion, everything else dropped.  
Checked: The identical two-clause filter is the lane idiom — build_app_data.py:85-86, build_crossed_fraction.py:66, build_galectin_layer.py:98, aggregate_islets.py:55 and :90-91 all use it — and cell_region_20um carries exactly three values, so on donor 115 it drops 729,312 'tissue' cells plus 8,600 Excluded, 8,360 of which sit inside islet cores (15.5% of that donor's core cells).


### `scripts/senior/compute_trajectory_senior.py`

**`nd.size < 10 (root-pool fallback in _robust_root)` = `10`** — scripts/senior/compute_trajectory_senior.py:113 — *ARBITRARY*  
Controls: If fewer than 10 ND islets clear ROOT_MIN_CELLS, the cell floor is abandoned and the root is taken from all ND islets — reopening the exact failure mode (a tiny saturated islet as root) the floor was added to prevent.  
Checked: Confirmed at line 113 inside _robust_root; the docstring says only 'Falls back to unfiltered ND if too few', the fallback has no counterpart in reaggregate_islets.py (which has no floor to fall back from), and 61f2006 states nothing about it.

**`g_root` = `(ds[root] == "ND") and (tc[root] >= ROOT_MIN_CELLS)`** — scripts/senior/compute_trajectory_senior.py:188 — *ARBITRARY* **[claim overturned]**  
Controls: Hard gate: the chosen root must be an ND islet with at least ROOT_MIN_CELLS cells; it can fail only via the line-113 fallback that drops the cell floor, and failure writes *.rejected.h5ad.  
Checked: The docstring attributes the medoid root to the Junior recipe, but reaggregate_islets.py has neither a cell floor nor any root gate (lines 170-200 are selection only, no validation), so both halves of this gate originate in 61f2006 uncommented and the threshold it leans on (ROOT_MIN_CELLS=20) is itself unsourced.

**`INS (ins_idx = list(a.var_names).index("INS"))` = `"INS"`** — scripts/senior/compute_trajectory_senior.py:375 — *INHERITED_DEFAULT*  
Controls: The marker that defines the root candidate pool (line 117) and the gate_INS hard gate (line 184); a run whose var_names lack INS raises ValueError before anything is computed, so the engine as written cannot run marker-blind.  
Checked: Location corrected — the lookup is at line 375, not 376 (line 376 is blank, 378 is the low_cell_qc assignment); the keying is inherited from reaggregate_islets.py:171 ('ins_idx = list(adata_islet.var_names).index("INS")') and 6f045dc shows it predates the 2026-07 rewrite, which only moved it.

**`K_ROOT` = `50`** — scripts/senior/compute_trajectory_senior.py:54 — *INHERITED_DEFAULT*  
Controls: Size of the ND-by-INS candidate pool whose medoid becomes the DPT root (line 118); larger K pulls the root toward the centre of mass of INS-high ND islets.  
Checked: Confirmed inherited from scripts/reaggregate_islets.py:170 ('K_ROOT = 50') whose rationale at lines 165-169 is explicitly qualitative — 'a single-islet argmax is brittle to staining outliers; the medoid of a size-50 candidate pool is anchored by the centre-of-mass' — a number with no cohort and no procedure, so it never rises above inherited.

**`ROOT_MIN_CELLS` = `20`** — scripts/senior/compute_trajectory_senior.py:55 — *ARBITRARY* **[claim overturned]**  
Controls: Cell floor an islet must clear to be eligible as the root (line 112) and the second half of the gate_root hard gate (line 188); below 10 eligible ND islets the filter is dropped entirely (lines 113-115).  
Checked: The docstring's '(Junior recipe, reaggregate_islets.py)' attribution does NOT cover this: reaggregate_islets.py:170-180 selects the root from a bare 'nd_mask = donor_status == ND' with no cell filter at all, and grep for total_cells/min_cells in that file returns only the unrelated aggregation arg — so the floor originates here at 61f2006 (2026-07-02) with no measurement, and trajectory_islet.py:64 copied 20 from here two months later, also bare.


## Dormant — die with Harmony, insulin-rooting or cohort pooling


### `apps/islet/R/mod_trajectory_server.R`

**`drill-down colour-by and palette defaults` = `colour by "phenotype"; palette choices Original / High Contrast / Colorblind Safe / Maximum Distinction, selected "High Contrast"`** — apps/islet/R/mod_trajectory_server.R:1002-1006 — *ARBITRARY*  
Controls: How drill-down cells are coloured; the default palette is not the colour-blind-safe one.  
Checked: Unrecorded in both files and in git; the observation that the default is not the accessible option is confirmed by reading the choice vector.

**`drill-down structure legend colours` = `Islet #0066CC, +20µm #00CCCC, Nerve #CC00CC, Capillary #CC0000, Lymphatic #00AA00`** — apps/islet/R/mod_trajectory_server.R:1013-1032 — *INHERITED_DEFAULT* **[claim overturned]**  
Controls: Hard-coded outline colours (and the hard-coded "+20µm" label) in the fallback legend when no per-cell data exists.  
Checked: Not unknown: the identical hexes are the renderer's structure palette at segmentation_helpers.R:561 and :605 (and the same block is copied into mod_plot_server.R:1285), so this legend is a hand-copy of that map -- one that omits the Non-endocrine class (#8C6BB1) the renderer also draws.

**`drill-down panel plot heights` = `segmentation view "550px"; composition summary "300px"`** — apps/islet/R/mod_trajectory_server.R:1042 and 1054 — *ARBITRARY*  
Controls: Size of the embedded drill-down plots.  
Checked: No comment; git -S shows the 550px arriving with layout commits (c012b10, 46c668e) that give no basis.

**`outlier-table threshold fallback` = `if (!is.finite(thr)) thr <- 3`** — apps/islet/R/mod_trajectory_server.R:1116-1117 — *ARBITRARY*  
Controls: The number printed in the table header and body when the shared threshold is unreadable; can disagree with the threshold actually used by flag_outliers.  
Checked: Unrecorded; it re-states the 3 independently of ieCore's own fallback, so a future change to one leaves the printed caption describing a rule that was not applied.

**`UMAP point style` = `geom_point(alpha = 0.6, size = 1.0)`** — apps/islet/R/mod_trajectory_server.R:1284 and 1323 — *ARBITRARY*  
Controls: Opacity and size of every islet in the two UMAP panels.  
Checked: No comment; git -S puts the tuning in 2f0e388 ("Refine UMAP visualizations: adjust point sizes..."), a polish commit with no stated basis.

**`UMAP axis expansion and aspect` = `expansion(mult = 0.02) on both axes; coord_fixed(); axis text/ticks/grid blanked`** — apps/islet/R/mod_trajectory_server.R:1286-1296 and 1324-1338 — *ARBITRARY* **[claim overturned]**  
Controls: Framing of the UMAP panels; coord_fixed is what forced the stacked layout.  
Checked: The quoted UI comment is an argument about card geometry ("~190px of width each against 500px of height") supporting the stacking, not a measurement of these values -- coord_fixed follows from a stated principle and the 0.02 expansion has no record at all.

**`UMAP legend key size override` = `guides(color = guide_legend(override.aes = list(size = 3)))`** — apps/islet/R/mod_trajectory_server.R:1292 — *ARBITRARY*  
Controls: Legend swatch size on the donor-status UMAP.  
Checked: Unrecorded; it exists only because the plotted points are size 1.0, itself an unrecorded value.

**`UMAP PNG download geometry` = `ISLET_PNG_DIMS = 10 x 7 in, 300 dpi, bg = "white"`** — apps/islet/R/mod_trajectory_server.R:1301-1306 and 1350-1355 — *INHERITED_DEFAULT* **[claim overturned]**  
Controls: Size, resolution and background of the saved UMAP PNGs; the file deliberately differs from the screen.  
Checked: The cited comment justifies bg = "white" by a rendering fact (a transparent PNG inverts on a dark slide) -- a reasoned convention, not a measurement -- and the dimensions are simply the app-wide constant at 00_globals.R:557 used by every download handler.

**`strip missing-bin colour` = `"#BBBBBB"`** — apps/islet/R/mod_trajectory_server.R:1432 — *ARBITRARY*  
Controls: Colour for a bin with no dominant status or no donor colour; visually similar to a faintly tinted bin.  
Checked: No record; last touched by the styling commit 61803f3 ("Port the Modernist design to Senior"), which changes skins and gives no rationale.

**`active_tab default` = `active_tab = reactive("Trajectory")`** — apps/islet/R/mod_trajectory_server.R:15 — *ARBITRARY* **[claim overturned]**  
Controls: Gates the embedded segmentation drill-down panel (line 953); the panel renders only when this equals "Trajectory".  
Checked: Not unknown: commit 92e8f6d ("add active tab reactive to prevent duplicate output IDs") introduced it and the drill-down uses NON-namespaced input/output ids shared with the Plot tab, so the default is a chosen open-gate fallback for standalone use -- app.R:409 always passes the real reactive.

**`donor_colors_reactive default` = `reactive(DONOR_COLORS_BRIGHT)`** — apps/islet/R/mod_trajectory_server.R:16 — *INHERITED_DEFAULT*  
Controls: Point, trend-line, legend and pseudotime-strip colours for ND/Aab+/T1D.  
Checked: Confirmed constant at 00_globals.R:35 (domain_group_colors("Bright")), but note it disagrees with the active domain: run20_default.yml:41 sets default_palette: Okabe-Ito, and app.R:409 overrides this default with active_donor_colors, so the Bright default only fires standalone/in tests.

**`multi-heatmap colour ramp` = `low "#2166ac", mid "#f7f7f7", high "#b2182b", midpoint = 0, na.value "#f7f7f7"`** — apps/islet/R/mod_trajectory_server.R:1608-1610 — *ARBITRARY* **[claim overturned]**  
Controls: Diverging ramp centred on a marker's own across-bin mean; an empty bin is painted the same colour as z = 0.  
Checked: Not a library default -- I checked the installed ggplot2 and scale_fill_gradient2 defaults to low = muted("red") (#832424), mid = "white", high = muted("blue") (#3A3A98), i.e. the opposite polarity -- these three hexes are hand-typed ColorBrewer RdBu entries (RdBu-11 positions 10, 6 and 2), and na.value repeating the mid colour makes an empty bin indistinguishable from a bin at the marker's mean.

**`multi-heatmap render height formula` = `max(220, 60 + n * 36)`** — apps/islet/R/mod_trajectory_server.R:1621-1624 — *ARBITRARY*  
Controls: On-screen pixel height per selected marker row.  
Checked: git -S traces it to 3f32595, a font/PNG standardization commit that records no basis for 36 px per row or the 220 px floor.

**`multi-heatmap download height formula` = `h <- max(6, 1 + n * 0.4) inches; width/dpi from ISLET_HEATMAP_PNG_DIMS (12 in, 300 dpi)`** — apps/islet/R/mod_trajectory_server.R:1628-1635 — *ARBITRARY* **[claim overturned]**  
Controls: Saved-PNG geometry for the heatmap; differs from the on-screen aspect.  
Checked: Only half is inherited: width and dpi do come from 00_globals.R:558 as claimed, but the height is a local formula from 3f32595 with no record, and it overrides the constant's own height = 8.

**`leiden_cols merged from prepared()$comp` = `c("leiden_0.3", "leiden_0.5", "leiden_0.8", "leiden_1.0")`** — apps/islet/R/mod_trajectory_server.R:390-391 — *INHERITED_DEFAULT* **[claim overturned]**  
Controls: Which precomputed Leiden resolutions can colour the scatter; other resolutions are not joined even if present.  
Checked: Not unknown: the producer explicitly mirrors the app, not the reverse -- trajectory_islet.py:84-85 says its LEIDEN_RESOLUTIONS match "the legacy column names (leiden_0.3 ... leiden_1.0) that mod_spatial_server.R greps for" -- and islet_leiden.parquet does carry exactly those four columns.

**`rank ties.method` = `ties.method = "average"`** — apps/islet/R/mod_trajectory_server.R:482 — *LIBRARY_DEFAULT*  
Controls: How tied pseudotime values share an x position.  
Checked: Confirmed in R 4.5.1: formals(rank)$ties.method is c("average","first","last","random","max","min"), so "average" is base R's own default written out.

**`traj_point_size (size-by) default mismatch` = `server fallback "uniform"; UI selected = "total_cells"`** — apps/islet/R/mod_trajectory_server.R:644 vs mod_trajectory_ui.R:39-43 — *ARBITRARY*  
Controls: Whether point area encodes cell count, islet diameter, or nothing; the two defaults disagree, so the first render before the input registers differs from the settled one.  
Checked: Confirmed both values by reading; neither has a comment and no commit message mentions the disagreement.

**`jitter RNG seed` = `set.seed(42)`** — apps/islet/R/mod_trajectory_server.R:651 — *ARBITRARY*  
Controls: Makes the jitter reproducible across renders -- and resets the session's global RNG stream on every scatter draw.  
Checked: No comment names it and it predates the module split (230af4d); 42 is the usual placeholder, and the side effect on the global stream is real since set.seed is called unconditionally inside renderPlotly.

**`point size transform for cell count` = `df$size_cells <- sqrt(df$total_cells)`** — apps/islet/R/mod_trajectory_server.R:701 — *ARBITRARY* **[claim overturned]**  
Controls: Compresses the cell-count range into point area.  
Checked: Not unknown: commit 4a1b5c9 introduced it in the same change as the LOESS weights, with a stated purpose ("well-measured islets drive the biological signal") but no measurement -- sqrt is the conventional area encoding, chosen not derived.

**`geom_point(stroke=)` = `stroke = 0`** — apps/islet/R/mod_trajectory_server.R:709, 711 — *ARBITRARY*  
Controls: Removes point outlines, which at high alpha changes the apparent density of overlapping islets.  
Checked: Unrecorded in both comments and git; it is a non-default (ggplot2's stroke default is 0.5).

**`ghost (excluded-outlier) marker style` = `shape = 4, size = point_size * 1.1, stroke = 0.9, colour = "grey45", alpha = 0.55`** — apps/islet/R/mod_trajectory_server.R:719-729 — *ARBITRARY* **[claim overturned]**  
Controls: How excluded points are drawn; they remain on the plot but not in the fits.  
Checked: The cited comment justifies the practice of drawing exclusions rather than deleting them -- a display principle, not a measurement -- and gives no basis for any of the five numbers.

**`donor_id palette rule` = `n_donors <= 12 -> brewer.pal(min(12, max(3, n)), "Paired"); else rainbow(n, s = 0.8, v = 0.8)`** — apps/islet/R/mod_trajectory_server.R:745-752 (legend copy at 1195-1200) — *ARBITRARY* **[claim overturned]**  
Controls: Per-donor colours in the scatter and its legend; above 12 donors the palette switches to rainbow.  
Checked: Only the 12 is a library fact (RColorBrewer::brewer.pal.info says Paired has maxcolors 12); the branch and the rainbow(s=0.8, v=0.8) fallback are hand-written with no record, so this is not a library default -- the deployed file has 10 donors, so the brewer branch is what runs.

**`Leiden palette rule` = `n_clust <= 12 -> brewer.pal(max(3, min(12, n)), "Set3"); else scales::hue_pal()`** — apps/islet/R/mod_trajectory_server.R:757-760 (legend copy at 1245-1248) — *ARBITRARY* **[claim overturned]**  
Controls: Cluster colours when colouring by a Leiden resolution.  
Checked: Same shape: brewer.pal.info confirms Set3 maxcolors 12 (and colorblind = FALSE), but the threshold rule and the hue_pal fallback are hand-written and unrecorded, so nothing here is a package default.

**`size scale range, islet diameter` = `c(0.5 * point_size, 2.5 * point_size)`** — apps/islet/R/mod_trajectory_server.R:769 — *ARBITRARY*  
Controls: How strongly islet diameter is expressed as point size.  
Checked: No comment and no commit rationale; the multipliers differ from the cell-count pair four lines later for no stated reason.

**`size scale range, cell count` = `c(0.3 * point_size, 3.0 * point_size)`** — apps/islet/R/mod_trajectory_server.R:773 — *ARBITRARY*  
Controls: How strongly sqrt(cell count) is expressed as point size.  
Checked: Unrecorded; introduced alongside the sqrt encoding in 4a1b5c9, which justifies the encoding but never the range.

**`trend line widths and overall-trend colour` = `size = 1 and color = "black" (overall); size = 0.8 (per status)`** — apps/islet/R/mod_trajectory_server.R:781-782 (overall), 800 and 807 (per status) — *ARBITRARY*  
Controls: Visual weight of the trend curves; the black overall curve is what leaks "colour: black" into the hover text (stripped at line 845).  
Checked: The only related comment describes the hover-text cleanup the black causes, not the choice itself; unrecorded in git.

**`scatter axis limits` = `coord_cartesian(xlim = c(0,1), ylim = range(value) +/- 0.05 * diff(range(value)))`** — apps/islet/R/mod_trajectory_server.R:840 — *ARBITRARY*  
Controls: Fixes the x-axis to the full rank range and pads y by 5% of the data range in each direction.  
Checked: Unrecorded (predates 230af4d); the xlim is only meaningful because of the rank normalization, and the 5% pad has no stated basis.

**`ggplotly tooltip fields` = `tooltip = c("x", "y", "colour", "text")`** — apps/islet/R/mod_trajectory_server.R:843 — *ARBITRARY*  
Controls: What the hover shows -- pseudotime, value, colour group and the cell-count string; islet id is not among them.  
Checked: No comment or commit explains the set; confirmed that the islet identity travels only in `key` (line 666), which is used for clicks, not hover.

**`plotly margins` = `list(l = 80, r = 20, t = 30, b = 70, pad = 5)`** — apps/islet/R/mod_trajectory_server.R:868 — *ARBITRARY*  
Controls: Plot area padding; the 80pt left margin is what keeps a five-figure marker-intensity axis label from clipping.  
Checked: git -S traces it to 3f32595 ("Standardize plot fonts and add publication-quality PNG downloads"), which records no measurement for the numbers.

**`plotly xaxis range` = `range = c(0, 1), fixedrange = FALSE`** — apps/islet/R/mod_trajectory_server.R:869 — *ARBITRARY*  
Controls: Restates the 0-1 pseudotime window in the plotly layout (the ggplot coord is not carried across by ggplotly).  
Checked: The CLAUDE.md rule the extractor cites explains why a layout() overlay is needed at all, not why the range is 0-1; that value simply tracks the rank transform and is written by hand in both places.

**`click coordinate-fallback tolerance` = `tolerance <- 0.01`** — apps/islet/R/mod_trajectory_server.R:917-924 — *ARBITRARY*  
Controls: How close a click must land to a stored point before the drill-down opens that islet; the first match within the window wins, and pt and value are compared with the same absolute tolerance despite living on different scales.  
Checked: No record anywhere; the extractor's scale observation is confirmed by measurement -- pt is 0-1 after ranking while value reaches 29,608 for marker means, so |value - click$y| < 0.01 is effectively an exact-match requirement on the y term.

**`drill-down marker column exclusions` = `setdiff(colnames(cells), c("X_centroid","Y_centroid","phenotype","cell_region","Cell.Area","Cell Area","Nucleus.Area","Nucleus Area", EXCLUDED_MARKERS))`** — apps/islet/R/mod_trajectory_server.R:966-970 — *ARBITRARY* **[claim overturned]**  
Controls: Which per-cell columns are offered as the drill-down colour-by; area columns and the four excluded markers are removed by name.  
Checked: Nothing is inherited -- the geometry names are a hand-written list carrying both the dot and space spellings defensively, with no comment and no commit message -- and the EXCLUDED_MARKERS half inherits a list that itself no longer matches its cited authority.

**`drill-down layer defaults` = `Single Cells TRUE, Peri Boundary TRUE, Structures TRUE, Show peri-islet cells TRUE`** — apps/islet/R/mod_trajectory_server.R:996-998 and 1007 — *ARBITRARY*  
Controls: What the embedded segmentation panel draws when an islet is clicked.  
Checked: No comment or commit rationale; note these are NON-namespaced input ids shared with the Plot tab's copy of the same panel, which is why the whole panel is gated on active_tab.


### `apps/islet/R/mod_trajectory_ui.R`

**`show_outlier_table default` = `value = FALSE`** — apps/islet/R/mod_trajectory_ui.R:10 — *ARBITRARY* **[claim overturned]**  
Controls: Whether the list of excluded islets is shown under the scatter.  
Checked: The quoted comment ("The excluded-outlier TABLE is a display choice and stays") argues for keeping the control when the mirrored outlier inputs were removed in c1a788e -- it says nothing about the FALSE default, which has no record.

**`scatter/heatmap card column split` = `column(7) scatter card, column(5) UMAP + drill-down card; control row column(3) x 4`** — apps/islet/R/mod_trajectory_ui.R:15 and 83 — *ARBITRARY* **[claim overturned]**  
Controls: Page layout, which is what forced the UMAP stacking decision above.  
Checked: "narrowed from col-8" is a change note recording what happened, not a measurement justifying 7/5; no commit supplies one.

**`traj_alpha slider` = `min = 0, max = 0.95, value = 0.4, step = 0.05`** — apps/islet/R/mod_trajectory_ui.R:46-47 (used as alpha_val <- 1 - x at line 646) — *INHERITED_DEFAULT* **[claim overturned]**  
Controls: Point opacity, hence how much overplotting is visible in dense pseudotime regions.  
Checked: Not unique to this tab: commit 019275f ("standardize point defaults ... Set default pt_size=3.0, pt_alpha=0.6") set it app-wide and mod_plot_ui.R:94 carries the identical slider, so the value is inherited from the app-wide convention -- which itself rests on no measurement.

**`traj_point_size_slider` = `min = 0.5, max = 5.0, value = 3.0, step = 0.1`** — apps/islet/R/mod_trajectory_ui.R:48-49 (server fallback 3.0, line 647) — *INHERITED_DEFAULT* **[claim overturned]**  
Controls: Base point size, and the multiplier base for both size-scale ranges and the ghost marker size.  
Checked: Same commit 019275f standardization; mod_plot_ui.R:92 has value = 3.0 too (with a different min), so the 3.0 is an app-wide convention rather than a per-tab choice -- still with no measurement behind it.

**`scatter plot height` = `height = 650`** — apps/islet/R/mod_trajectory_ui.R:67 (mounted) and mod_trajectory_server.R:575,577 (unmounted rdeck renderUI) — *ARBITRARY*  
Controls: Vertical extent of the pseudotime scatter.  
Checked: No comment, and the value predates the module split; the only mounted copy is the UI one, since output$traj_scatter_container is never referenced by any ns() call.

**`strip plot height and margins` = `height = 100; plot.margin = margin(5, 20, 10, 60, "pt")`** — apps/islet/R/mod_trajectory_ui.R:69 and mod_trajectory_server.R:1457 — *ARBITRARY*  
Controls: Vertical size and left offset of the strip; the 60pt left margin is hand-tuned to sit under the scatter's 80pt plotly margin and does not track it.  
Checked: Unrecorded in comment or commit, and the mismatch with the scatter's l = 80 (line 868) is visible in the code itself.

**`UMAP panel heights` = `height = 320 (both UMAP plotOutputs)`** — apps/islet/R/mod_trajectory_ui.R:95 and 104 — *ARBITRARY* **[claim overturned]**  
Controls: On-screen size of the two stacked UMAP cards.  
Checked: The comment measures the problem it replaced (two abreast at ~190 px wide against 500 px tall) but derives no number for the fix -- "the height comes down to match so each is about square" is an eyeballed target, and commit c883cad shows the same slot was 300 then 500 before landing on 320.


### `core/phenocycler/scripts/typing/summarise_singlecell_trajectory.py`

**`MIN_SAMPLED_CELLS` = `5`** — core/phenocycler/scripts/typing/summarise_singlecell_trajectory.py:38 (applied at :60-61) — *ARBITRARY*  
Controls: An islet region with fewer than five ordered cells gets NaN instead of a median. On the current artifact this is decisive, not marginal.  
Checked: The comment is a rationale with no number, cohort or procedure, and the same unmeasured 5 appears as MIN_CELLS in summarise_crossed_fraction.py:39 from the same commit (so neither inherits from the other); I measured its effect — it blanks 3,157 of 4,426 core islets (71%) and 4,821 of 6,966 peri islets, leaving only 1,020 of 7,179 islets with a peri-minus-core difference.

**`--no-merge = action="store_true" (default False, i.e. it does merge)` = `False — the script rewrites data/development/typing/app/parquet/islets.parquet in place by default`** — core/phenocycler/scripts/typing/summarise_singlecell_trajectory.py:43-45 (acted on at :77-87) — *INHERITED_DEFAULT* **[claim overturned]**  
Controls: Whether the sc_pseudotime_* / sc_n_* columns are written into the app's islets.parquet, dropping any existing ones first.  
Checked: A basis does exist: the identical flag with the identical default is in summarise_crossed_fraction.py:46, and README.md documents the ordering it serves ('build_parquet_views.py before the two summarise_* scripts ... Re-running build_parquet_views.py drops their columns'); it has already run — islets.parquet currently carries all seven sc_ columns of its 104.

**`non-finite pseudotime filter` = `cells = cells[np.isfinite(cells["pseudotime"])]`** — core/phenocycler/scripts/typing/summarise_singlecell_trajectory.py:48 — *UNKNOWN*  
Controls: Cells the curve fit left without a finite pseudotime are dropped before the per-islet median and before the sc_n_* counts. Currently drops nothing.  
Checked: No comment here and nothing in commit 6b8c346 about it; I measured the input — all 49,999 ordered cells have a finite pseudotime (range 0 to 0.2442) — so this is a silent guard that has never fired, and if it ever did it would shrink sc_n_* without printing anything, unlike the R script which at least reports the non-finite count.

**`islet statistic = median (agg ["median", "mean", "size"], median written to sc_pseudotime_{core,peri})` = `median`** — core/phenocycler/scripts/typing/summarise_singlecell_trajectory.py:55-61 — *ARBITRARY*  
Controls: How an islet's cells are collapsed to one position on the ordering; the mean is computed and carried into islets.parquet as sc_pseudotime_mean_* but never used.  
Checked: The docstring argues the case but records no comparison of the two statistics on this cohort, and the sampled counts make the argument load-bearing rather than cosmetic — the median core islet has only 2 sampled cells (p75 = 5), so for most islets the median and the mean are computed over a handful of cells.


### `core/phenocycler/scripts/typing/trajectory_islet.py`

**`mannwhitneyu(alternative=...)` = `"greater" (scipy default is "two-sided")`** — core/phenocycler/scripts/typing/trajectory_islet.py:319-321 — *INHERITED_DEFAULT*  
Controls: Makes the reported T1D-vs-ND p-value one-sided in the legacy lane's direction; reported, never gated here.  
Checked: Byte-identical to compute_trajectory_senior.py:172, where it fed a disease-ordering gate this file's docstring disavows; scipy 1.18.0 in .venv defaults to "two-sided", and the deployed H5AD records mwu_T1D_gt_ND_p = 0.0 - per donor the branch cannot fire, since one donor carries one status.

**`DONOR_ETA2_CEILING` = `0.15`** — core/phenocycler/scripts/typing/trajectory_islet.py:89 — *INHERITED_DEFAULT*  
Controls: Ceiling on within-status donor eta-squared; pseudotime_is_not_mostly_donor fails the build above it (line 293-294).  
Checked: compute_trajectory_senior.py:58 carries the same value for the same purpose (ETA2_GATE = 0.15) and neither file measures it; per donor it becomes inert, since a single-donor build leaves within_status_donor_eta2 NaN and gated_statuses empty.

**`DONOR_ETA2_MIN_DONORS` = `3`** — core/phenocycler/scripts/typing/trajectory_islet.py:95 — *ARBITRARY*  
Controls: How many donors a status needs before its eta-squared gates rather than is merely reported (line 283-285).  
Checked: Argued in prose, never measured, and raised from the legacy min_donors=2 (compute_trajectory_senior.py:72) - the raise is what admits this build: the deployed H5AD records within_status_donor_eta2 T1D = 0.4625 over 2 donors, three times the 0.15 ceiling, so at the inherited value of 2 the pseudotime_is_not_mostly_donor gate fails and the file is written .rejected.h5ad instead of deployed.


### `core/phenocycler/scripts/typing/trajectory_singlecell.R`

**`pseudotime column selection pt[, 1]` = `out$pseudotime <- pt[, 1]`** — core/phenocycler/scripts/typing/trajectory_singlecell.R:97 — *UNKNOWN*  
Controls: Which lineage becomes the single reported pseudotime if slingshot ever infers more than one. Currently inert: there is exactly one.  
Checked: No comment on the line, and the only commit touching the file (6b8c346) says nothing about lineage selection; the artifact confirms it is currently a no-op — singlecell_pseudotime.parquet has exactly one pseudotime_lineage_ column and the manifest records n_lineages 1 — but the selection is silent, so a two-lineage fit would discard the second without any record.


### `core/phenocycler/scripts/typing/trajectory_singlecell_features.py`

**`--donors = default=None -> COHORT.resolve(None) = the 10 CLEAN donors` = `None -> ["115", "6374", "6436", "6442", "6461", "6476", "6482", "6493", "6539", "6593"]`** — core/phenocycler/scripts/typing/trajectory_singlecell_features.py:118 (resolved at :121, cohort.py:23-24) — *INHERITED_DEFAULT* **[claim overturned]**  
Controls: Which donors' cells are concatenated into one feature matrix; the default pools the clean cohort into a single embedding.  
Checked: The value at this line is literally `default=None` — the measurement the extractor cites lives in cohort.py, whose docstring does state a procedure (hormone-high clusters outside every islet annotation, six settings: linking distance 15/25/40 um x percentile 99/99.5) and a table, but the sweep that produced it is not in the repo and the maintainer's own 2026-09-12 record says 'zero on nine of thirteen batch-1 donors' where CLEAN lists ten.

**`--include-marginal = action="store_true" (default False)` = `False`** — core/phenocycler/scripts/typing/trajectory_singlecell_features.py:119 — *INHERITED_DEFAULT* **[claim overturned]**  
Controls: Whether 6521, 6551 and 6380 join the build; off by default, so the pooled matrix has no Aab+ donor. Ignored entirely when --donors is given explicitly.  
Checked: The default False encodes cohort.py's COHORT_LIMITS policy sentence ('that is a maintainer decision, not a default') which rests on measured per-donor deficits, but the measurement is in cohort.py rather than here, and cohort.resolve() short-circuits on an explicit --donors list so this flag has no effect in a per-donor invocation.

**`PANEL_IDENTITY` = `("CD38", "b_Catenin1", "bCatenin1", "CD57")`** — core/phenocycler/scripts/typing/trajectory_singlecell_features.py:61 (applied at :79, written to the manifest at :164) — *ARBITRARY* **[claim overturned]**  
Controls: CORRECTED: nothing. It removes zero markers. None of the four is in the set feature_markers() builds, so `set(markers) - set(PANEL_IDENTITY)` is a no-op — yet singlecell_features.json still records "held_out": [CD38, b_Catenin1, bCatenin1, CD57], asserting a holdout that never happened.  
Checked: I re-ran feature_markers() against the real logs: the derived set is 36 markers and contains none of CD38/CD57/b_Catenin1 (and "bCatenin1" is not a marker name anywhere — the 59-row scale table spells it b_Catenin1), and the "coarse embedding" the comment attributes the holdout to appears nowhere in the checkout (`grep -rn "coarse embedding"` returns only this file and trajectory_singlecell.R, both from the single commit 6b8c346).


### `scripts/senior/compute_trajectory_senior.py`

**`metric="cosine" (sc.pp.neighbors, core graph)` = `"cosine"`** — scripts/senior/compute_trajectory_senior.py:138 — *INHERITED_DEFAULT*  
Controls: The distance the islet kNN graph — and hence the diffusion map and every DPT ordering — is built on; scanpy's default is euclidean, so this is a deliberate non-default.  
Checked: No basis in this file; identical at reaggregate_islets.py:148-149, also bare — and the only lane that actually measured the choice reached the opposite conclusion for its own representation (trajectory_islet.py:51-56: cosine 'put the median pair at 0.97 similarity ... the ordering did not survive a change of neighbourhood size (rho 0.31)'), so nothing here establishes cosine for the scVI latent.

**`sc.pp.neighbors — relied-on library defaults random_state / method / transformer` = `random_state=0, method="umap", transformer=None`** — scripts/senior/compute_trajectory_senior.py:138 (also :334, :353, :425) — *LIBRARY_DEFAULT*  
Controls: Fixes the graph seed and the UMAP fuzzy-simplicial-set kernel; with transformer=None and a non-euclidean metric scanpy switches from exact brute-force to approximate pynndescent kNN once n_obs >= 4096, so the same code builds an exact graph for a small per-donor run and an approximate one for the pooled cohort.  
Checked: No comment on any of the four calls; defaults confirmed in scanpy 1.11.5 (scvi-env, neighbors/__init__.py:72-86) and the 4096/8192 switch read directly at lines 647-650 there — the extractor's cited line 650 is right for 1.11.5, though in the repo .venv's scanpy 1.12.3 it has moved to line 737 and `metric` now defaults to None (resolved to 'euclidean' at line 209-210) rather than the literal 'euclidean'.

**`n_comps (sc.tl.diffmap)` = `10`** — scripts/senior/compute_trajectory_senior.py:141 (also :336) — *INHERITED_DEFAULT*  
Controls: Number of diffusion components computed; DC1/DC2 are taken from it and the bare sc.tl.dpt then consumes n_dcs=10, so this exactly saturates the DPT input and lowering it would silently truncate DPT.  
Checked: Bare call with no comment; identical at reaggregate_islets.py:202 ('sc.tl.diffmap(adata_islet, n_comps=10)') with no basis, and confirmed non-default — scanpy 1.11.5 and 1.12.3 both default diffmap to n_comps=15.

**`sc.tl.dpt(a) — relied-on library defaults n_dcs / n_branchings / min_group_size` = `n_dcs=10, n_branchings=0, min_group_size=0.01`** — scripts/senior/compute_trajectory_senior.py:142 (also :337) — *LIBRARY_DEFAULT*  
Controls: How many diffusion components the pseudotime is computed over (10, matching diffmap) and that no branching/segment detection runs — the output is a single monotone pseudotime, never a branch assignment.  
Checked: Call is bare at both sites; confirmed against the real runtime env (scanpy 1.11.5 in scvi-env, tools/_dpt.py:41-50: n_dcs=10, n_branchings=0, min_group_size=0.01, allow_kendall_tau_shift=True) and identical in the repo .venv's scanpy 1.12.3, so the default is stable across both.

**`g_ins threshold (r_ins < -0.3)` = `-0.3`** — scripts/senior/compute_trajectory_senior.py:184 — *ARBITRARY*  
Controls: Hard gate: Spearman(pseudotime, INS) must be below -0.3 or the run writes *.rejected.h5ad and exits 1; it also filters the theta-sweep fallback branch (line 242).  
Checked: git log -L 184,189 shows all five gate expressions arriving together, uncommented, in 61f2006, whose message says only 'Introduced validation gates to ensure data integrity' — no number, cohort or procedure; the later lane explicitly refuses this gate (trajectory_islet.py:17-22: 'A check that passes by reproducing a previous finding measures the expectation, not the build').

**`g_gcg threshold (r_gcg > 0.1)` = `0.1`** — scripts/senior/compute_trajectory_senior.py:185 — *ARBITRARY*  
Controls: Hard gate: Spearman(pseudotime, GCG) must exceed +0.1 when GCG is in var_names; absent GCG passes vacuously, and failure sends the run to *.rejected.h5ad.  
Checked: Same uncommented block in 61f2006 (git log -L 184,189) with no stated basis; confirmed vacuous-when-absent by reading lines 165-167 and 185.

**`g_order thresholds (mwu_p < 0.05, alternative="greater")` = `0.05, alternative="greater"`** — scripts/senior/compute_trajectory_senior.py:186 (statistic at :172) — *ARBITRARY*  
Controls: Hard gate: median pseudotime T1D must exceed ND and the one-sided Mann-Whitney p must be under 0.05, else the run writes *.rejected.h5ad; the one-sided direction encodes the expected ND<T1D result as a pass condition.  
Checked: Conventional alpha with no stated basis, introduced uncommented in 61f2006; note this gate is not dead per-donor but hard-fails there — with a single donor one of nd_pt/t1d_pt is empty, mannwhitneyu raises ValueError (caught at line 173-174 to NaN) and med_t1d > med_nd is False, so a per-donor run always exits 1 unless --force.

**`z-score epsilon in _struct_block` = `1e-8`** — scripts/senior/compute_trajectory_senior.py:281 — *ARBITRARY*  
Controls: Denominator floor when z-scoring structural proportions and counts, so a constant column becomes ~0 rather than NaN; reached only in combined mode.  
Checked: Confirmed as an inline literal inside the lambda at line 281 with no comment and no counterpart in reaggregate_islets.py's struct block; a conventional numerical guard, but nothing in the file or 61f2006 states it.

**`--force` = `store_true, default False`** — scripts/senior/compute_trajectory_senior.py:366 — *ARBITRARY* **[claim overturned]**  
Controls: The only CLI parameter: off, a hard-gate failure writes islets_core_senior.rejected.h5ad and exits 1 (lines 440-445); on, the deployed file is overwritten anyway and a gate-failing combined column is written too (line 406).  
Checked: Downgraded from UNKNOWN: the provenance is legible, not unknown — the flag and the fail-before-write design both arrive in 61f2006 with the docstring stating the intent ('hard gates that RAISE ... BEFORE overwriting the deployed file') and the help text scoping it ('debug only'), which is a design decision taken here with no measurement and no upstream source, i.e. arbitrary in the sense used by the other entries.

**`n_comps (sc.pp.pca, display UMAP)` = `15`** — scripts/senior/compute_trajectory_senior.py:424 — *INHERITED_DEFAULT*  
Controls: PCs of the raw marker matrix feeding the display-only UMAP; it never enters DPT, but the call errors if a run has fewer than 15 islets or 15 markers, aborting before the write.  
Checked: The in-file comment at lines 421-422 attributes the block ('Junior recipe.') but not the number; confirmed inherited from reaggregate_islets.py:448 ('sc.pp.pca(adata_islet, n_comps=15)'), and confirmed non-default — scanpy signs n_comps: int | None = None in both 1.11.5 and 1.12.3.

**`F = SD / "islets_core_senior.h5ad" (SD = parents[2] / "data")` = `<repo>/data/islets_core_senior.h5ad`** — scripts/senior/compute_trajectory_senior.py:48 — *ARBITRARY* **[claim overturned]**  
Controls: The single hardcoded input, also overwritten in place on a pass (line 447) or written as F.with_suffix('.rejected.h5ad') on a gate failure (line 441); there is no output path argument.  
Checked: Downgraded from UNKNOWN: git log shows the path was set deliberately in 86aa55d (2026-06-17, 'Relocate data into the repo (was ~/IO60panc2nd/senior_data)'), so its origin is documented in the commit record even though the line itself carries no comment — a project convention, not an unknown.

**`PERI = SD / "islets_peri_senior.h5ad"` = `<repo>/data/islets_peri_senior.h5ad`** — scripts/senior/compute_trajectory_senior.py:49 — *ARBITRARY* **[claim overturned]**  
Controls: Its presence decides whether combined mode runs at all (line 289); if missing, the combined pseudotime is silently skipped and the app falls back to the core axis.  
Checked: Downgraded from UNKNOWN for the same reason as F — the relocation to <repo>/data is recorded in 86aa55d (2026-06-17) and the constant was added alongside the Step-5b port in 61f2006; the only in-file text is the runtime message at line 290, but the provenance is traceable, so 'unknown' overstates it.

**`N_NEIGHBORS` = `30`** — scripts/senior/compute_trajectory_senior.py:53 — *INHERITED_DEFAULT*  
Controls: k for the kNN graph under both the core DPT (line 138) and the combined DPT (line 334); sets the locality of the diffusion operator and therefore the whole ordering.  
Checked: Bare assignment with no comment; confirmed identical at scripts/reaggregate_islets.py:148 ('sc.pp.neighbors(adata_islet, n_neighbors=30, use_rep="X_scVI_harmony_pt", metric="cosine")', 2026-04-25, no basis stated) and re-copied into core/phenocycler/scripts/typing/trajectory_islet.py:62 on 2026-09-14 — scanpy's default is 15 in both 1.11.5 and 1.12.3, so this is a carried non-default with no measurement behind it.

**`ALPHA_STRUCT` = `0.15`** — scripts/senior/compute_trajectory_senior.py:59 — *INHERITED_DEFAULT*  
Controls: Weight on the z-scored structural-phenotype block (proportions + counts, lines 313 and 317) when concatenated to core+peri scVI for the combined-mode pseudotime.  
Checked: The sweep it inherits is a genuine measurement — reaggregate_islets.py:223-234 gives an alpha/span/order/INS/CD8a/PDPN table dated 2026-04-25 (matching that file's commit 54aa678) with 'Disease ordering (ND<Aab+<T1D) breaks at alpha >= 0.20' — but it was measured on the Junior legacy-CODEX lane, is not re-measured here, and the script named as its producer (alpha_sweep.py) is not in this repo, so it cannot be reproduced from this checkout.

**`STRUCT_PHENOTYPES` = `["Endothelial", "Fibroblast", "Muscle", "Blood vessel", "Blood Vessel", "Neural", "Lymphatic"]`** — scripts/senior/compute_trajectory_senior.py:62-63 — *ARBITRARY* **[claim overturned]**  
Controls: Which phenotype columns enter the combined-mode structural block; names absent from uns['phenotype_names'] are silently dropped (line 304) and case variants de-duplicated (lines 306-310), so the block's width and the combined latent shift with the taxonomy.  
Checked: Only four of the seven entries are inherited — reaggregate_islets.py:235 is exactly ['Neural', 'Blood Vessel', 'Endothelial', 'Lymphatic'] — while Fibroblast, Muscle and the lower-case 'Blood vessel' variant were added here in 61f2006 with no stated basis, and the alpha=0.15 sweep that justified the block was run against the four-name list, not this seven-name one.

**`min_donors (within_status_batch_eta2 keyword default)` = `2`** — scripts/senior/compute_trajectory_senior.py:72 — *ARBITRARY*  
Controls: A status with fewer donors returns NaN and is excluded from gated_statuses (line 178); when gated_statuses is empty g_eta is vacuously True (line 187), so in a single-donor run the batch gate silently stops constraining anything.  
Checked: Added in 6e96cef (2026-07-03, git log -L 72,72) as a stated rule with no measurement, and actively superseded since: trajectory_islet.py:80-87 (2026-09-14) sets DONOR_ETA2_MIN_DONORS = 3 arguing 'eta-squared over two donors is one pairwise contrast, not a variance decomposition ... A status needs at least three donors before the number is a check rather than a description.'


## Dead in this path


### `apps/islet/R/mod_trajectory_server.R`

**`islet number regex` = `stringr::str_extract(combined_islet_id, "(?<=_Islet_)[0-9]{1,3}")`** — apps/islet/R/mod_trajectory_server.R:155-157 — *ARBITRARY* **[claim overturned]**  
Controls: Derives obs2$islet_key -- but nothing reads that column; result_df$islet_key comes from base_islet_id at lines 357-358.  
Checked: Unrecorded (predates 230af4d), and the extractor's failure mode is understated: str_extract is unanchored, so 575 of 7,302 deployed islets with numbers above 999 (max Islet_1304) are silently TRUNCATED to a colliding key ("Islet_130"), not set to NA -- harmless only because obs2$islet_key is never read anywhere in the module.

**`rdeck backend switch` = `if (is_rdeck_available()) rdeckOutput else plotlyOutput`** — apps/islet/R/mod_trajectory_server.R:573-598 — *ARBITRARY* **[claim overturned]**  
Controls: Would swap the ggplotly scatter for a WebGL deck whose x/y are aliased to pt/value; the deck path draws no trend lines, no jitter and no ghost outliers.  
Checked: Dead code, not an inherited default: mod_trajectory_ui.R:67 mounts plotlyOutput(ns("traj_scatter")) directly and nothing mounts traj_scatter_container, which tests/manual/mount_audit.R:16-23 records as one of exactly two accepted unmounted outputs ("the WebGL path is unreachable").


### `core/phenocycler/scripts/typing/trajectory_islet.py`

**`root ranking key and pool restriction` = `pool = (status == "ND") & (cells >= ROOT_MIN_CELLS); top = argsort(-endocrine_fraction)[:K_ROOT]; root = medoid of top in the ALR representation`** — core/phenocycler/scripts/typing/trajectory_islet.py:139-147 — *INHERITED_DEFAULT*  
Controls: The whole root-selection rule, and with it the root_is_ND and root_above_cell_floor gates and the re-root stability check.  
Checked: The recipe and the three-cell anecdote both come from the legacy lane (compute_trajectory_senior.py:16-18 and 109-118); the only substitution is the ranking key, endocrine fraction here in place of argsort(-ins) there, and the deployed H5AD shows it picked 6476_Islet_551.

**`root pool fallback threshold (pool.size < 10)` = `10`** — core/phenocycler/scripts/typing/trajectory_islet.py:140 — *INHERITED_DEFAULT*  
Controls: If fewer than 10 ND islets clear ROOT_MIN_CELLS the cell floor is dropped from root selection entirely.  
Checked: The same branch with the same constant is at compute_trajectory_senior.py:113 ("if nd.size < 10:") under a docstring at :109 saying "Falls back to unfiltered ND if too few"; neither file measures the 10.

**`re-root Spearman minimum overlap (ok.sum() > 10)` = `10`** — core/phenocycler/scripts/typing/trajectory_islet.py:255 — *ARBITRARY*  
Controls: How many islets must have a finite pseudotime under both roots before that re-root comparison contributes a rho to the median that gates stable_under_reroot.  
Checked: No comment, no commit text and no ancestor to inherit from - the legacy lane has no re-root loop at all, its only sensitivity check being the >=5-cell subset (compute_trajectory_senior.py:255-274) - so this threshold was introduced bare in 6b8c346.

**`K_ROOT` = `50`** — core/phenocycler/scripts/typing/trajectory_islet.py:64 — *INHERITED_DEFAULT*  
Controls: Size of the top-ranked ND pool the root medoid is taken over (line 144) and the number of alternative roots the re-root stability check iterates (line 250).  
Checked: Verbatim in both ancestors - compute_trajectory_senior.py:54 and reaggregate_islets.py:170 - with no measurement in any of the three, and the deployed H5AD records root_pool_size = 50.

**`ROOT_MIN_CELLS` = `20`** — core/phenocycler/scripts/typing/trajectory_islet.py:65 — *INHERITED_DEFAULT*  
Controls: Minimum core cells to enter the root candidate pool (line 139) and the threshold of the root_above_cell_floor gate (line 290).  
Checked: Carries no comment here; compute_trajectory_senior.py:55 holds the same value and its docstring credits the Junior recipe (reaggregate_islets.py), so the 20 is inherited twice over and measured nowhere.


### `core/phenocycler/scripts/typing/trajectory_singlecell_features.py`

**`span fallback 1.0 for non-finite or non-positive span` = `span = np.where(np.isfinite(span) & (span > 0), span, 1.0)`** — core/phenocycler/scripts/typing/trajectory_singlecell_features.py:105 — *INHERITED_DEFAULT*  
Controls: A marker whose cohort span is missing or zero would be left on raw arcsinh units instead of rescaled. It never fires on the current data.  
Checked: Line 105 is character-identical to subpanel_io.py:72 inside apply_scale, and I checked the loaded scale table: all 59 markers have finite, strictly positive span, so the fallback branch is never taken on this build.

**`markers absent from a donor's panel filled with np.nan` = `cells[missing] = np.nan`** — core/phenocycler/scripts/typing/trajectory_singlecell_features.py:111-112 — *INHERITED_DEFAULT*  
Controls: CORRECTED: nothing on this cohort. The loop body never executes — the `missing` list is empty for every donor.  
Checked: The comment matches the repo's standing rule, but I read the aligned schema for all 20 donors under data/development/celltype/uniform and every one carries all 38 features, so `missing` is always empty and the built matrix has 0 non-finite values — the extractor's claim that 'for a single-donor run every such column is entirely NaN' is false for every donor in this checkout, including the three marginal and seven batch-2 ones.


### `scripts/senior/compute_trajectory_senior.py`

**`non-finite eta^2 substitution in choose_theta fallback` = `1.0`** — scripts/senior/compute_trajectory_senior.py:245 and :246 — *ARBITRARY*  
Controls: When no theta passes the hard gates, a NaN within-status eta^2 (a single-donor status) is scored as 1.0 — the worst possible value — so a status that cannot be measured is penalised as if it were maximally batch-driven when picking the fallback theta.  
Checked: Confirmed as an inline literal with no comment, and the extractor's location needs correcting — the substitution appears twice, once for within_status_eta2_ND at line 245 and once for within_status_eta2_T1D at line 246, both inside the max(...) lambda.

**`worst_drop > 0.05 (sensitivity NOTE threshold)` = `0.05`** — scripts/senior/compute_trajectory_senior.py:271 — *ARBITRARY*  
Controls: Above this drop in within-status eta^2 the script prints a recommendation to switch to a graph-drop fallback; it is print-only and gates nothing.  
Checked: Confirmed print-only by reading lines 271-273 (the branch body is a single print) and the literal carries no comment; introduced with the rest of sensitivity_check in 61f2006 with no stated basis.

**`theta (hm.run_harmony in compute_combined_dpt)` = `0.1`** — scripts/senior/compute_trajectory_senior.py:326 — *INHERITED_DEFAULT*  
Controls: Harmony strength for the combined core+peri+structural latent — hardcoded, so combined mode is NOT subject to the theta sweep that chose the core latent's theta.  
Checked: Bare 'theta=0.1' here; confirmed inherited from reaggregate_islets.py:326-333 ('Light Harmony (theta=0.1) on the combined latent — same theta as core'), whose core sweep table at lines 125-128 does record measured values ('theta=0.1 : span=0.084 INS=-0.57 CD8a=+0.19 Aab+ eta^2=0.09 <- chosen') — a Junior-lane measurement carried over, not re-measured here, and no longer 'same theta as core' since the core theta is now swept.

**`theta (hm.run_harmony in compute_leiden)` = `2.0`** — scripts/senior/compute_trajectory_senior.py:348 — *INHERITED_DEFAULT*  
Controls: Harmony strength for the separate islet-Leiden latent feeding the app's Spatial tab; stronger donor correction than the trajectory graph uses, and it does not touch dpt_pseudotime.  
Checked: Confirmed inherited from reaggregate_islets.py:498 (theta=2.0) with the rationale at lines 132-134 ('Step 8 still uses theta=2 because cluster-level donor invariance is the goal there') — a stated goal, not a measurement — and it coincides with harmonypy's own effective default, confirmed in both installed versions ('theta ... Default is 2 for each batch').

**`n_neighbors / metric (sc.pp.neighbors, key_added="leiden_nn")` = `n_neighbors=15, metric="cosine"`** — scripts/senior/compute_trajectory_senior.py:353 — *INHERITED_DEFAULT*  
Controls: The separate graph the islet Leiden partitions are computed on (15, not the trajectory's 30); it does not feed DPT.  
Checked: No basis in this file; confirmed identical at reaggregate_islets.py:504-505 ('sc.pp.neighbors(adata_clust, use_rep="X_scVI_harmony", n_neighbors=15, metric="cosine")'), and 15 is also scanpy's own n_neighbors default in both 1.11.5 and 1.12.3, so only the cosine metric is a carried non-default.

**`leiden resolutions` = `[0.3, 0.5, 0.8, 1.0]`** — scripts/senior/compute_trajectory_senior.py:355 — *INHERITED_DEFAULT*  
Controls: The four resolutions written as leiden_0.3 … leiden_1.0 and offered in the app's Spatial tab selector.  
Checked: No basis in this file; confirmed identical at reaggregate_islets.py:519 ('resolutions = [0.3, 0.5, 0.8, 1.0]'), recorded downstream as the app's options at apps/islet/docs/spatial.md:34 ('Leiden resolution (0.3/0.5/0.8/1.0)'), and re-copied again into trajectory_islet.py:82 explicitly to match the legacy column names.

**`sc.tl.leiden — relied-on library defaults flavor / n_iterations / random_state` = `flavor="leidenalg", n_iterations=-1, random_state=0`** — scripts/senior/compute_trajectory_senior.py:356 — *LIBRARY_DEFAULT*  
Controls: Which Leiden implementation and seed produce the cluster labels; the labels are pinned only by the installed scanpy version.  
Checked: Call passes only resolution, key_added and neighbors_key; confirmed in the real runtime env (scanpy 1.11.5, tools/_leiden.py:46 flavor: Literal[...] = 'leidenalg', n_iterations=-1, random_state=0) — but the claim is version-fragile, because in the repo .venv's scanpy 1.12.3 the signature default is already flavor=None, resolved to 'leidenalg' only via a FutureWarning at _validate_flavor (lines 234-248) telling callers to pass flavor='igraph' and n_iterations=2.

**`n_neighbors / metric (sc.pp.neighbors, key_added="neighbors_viz")` = `n_neighbors=30, metric="euclidean"`** — scripts/senior/compute_trajectory_senior.py:425 — *INHERITED_DEFAULT*  
Controls: The display-only graph the visualization UMAP is laid out from; separate from the trajectory graph and not used by DPT.  
Checked: Nothing in this file beyond the block-level 'Junior recipe.' at line 422; confirmed identical at reaggregate_islets.py:449-450, and metric='euclidean' merely restates scanpy 1.11.5's own default while n_neighbors=30 is the carried non-default.

**`min_dist / spread (sc.tl.umap)` = `min_dist=0.5, spread=2.0`** — scripts/senior/compute_trajectory_senior.py:426 — *INHERITED_DEFAULT*  
Controls: Display-only UMAP layout spacing; spread=2.0 is non-default, min_dist=0.5 equals scanpy's default.  
Checked: Nothing in this file beyond 'Junior recipe.' at line 422; confirmed identical at reaggregate_islets.py:451 and echoed in that file's own log line at :516, and confirmed against scanpy 1.11.5 tools/_umap.py:44-45 (min_dist=0.5, spread=1.0) that only spread departs from stock.

**`THETA_SWEEP` = `[0.0, 0.1, 0.5, 1.0]`** — scripts/senior/compute_trajectory_senior.py:52 — *ARBITRARY*  
Controls: Candidate Harmony strengths swept for the core trajectory graph; the smallest theta passing the hard gates becomes the deployed latent (0.0 = raw X_scVI_mean, no correction). The chosen rep is also what _robust_root takes the medoid in.  
Checked: Confirmed at line 52 verbatim; the whole tunables block arrived in one commit 61f2006 (2026-07-02) whose message lists features and cites no number, cohort or procedure, and the only measured theta table in the repo (reaggregate_islets.py:125-128) covers raw/0.1/2.0 — so 0.5 and 1.0 appear nowhere as measured points and 2.0, the one measured non-light value, was dropped.

**`LOW_CELL_QC` = `10`** — scripts/senior/compute_trajectory_senior.py:56 — *ARBITRARY*  
Controls: Sets the boolean obs column low_cell_qc (line 378); nothing in this script filters, weights or reweights on it, so the written pseudotime is unchanged by it — it only ships a label into the app's H5AD.  
Checked: Comment gives only a description ('obs flag threshold (noisy small-islet marker)'); the one genuinely measured 10-cell floor in the repo (trajectory_islet.py:66-79, a floor-vs-rho(n_neighbors 30 vs 15) table) is a 2026-09-14 file postdating this one by two months and floors an ALR composition, not this flag.

**`SENSITIVITY_MIN_CELLS` = `5`** — scripts/senior/compute_trajectory_senior.py:57 — *ARBITRARY*  
Controls: The subset sensitivity_check (lines 254-275) uses to ask whether tiny islets pollute the manifold; print-only, it never gates or changes the written pseudotime.  
Checked: Confirmed print-only by reading lines 254-275 and 390 (the return value goes only into uns['sensitivity']); comment states the use, not a basis, and commit 61f2006 gives none.

**`ETA2_GATE` = `0.15`** — scripts/senior/compute_trajectory_senior.py:58 — *ARBITRARY*  
Controls: Ceiling on within-status donor eta^2; a chosen theta whose per-status donor eta^2 exceeds it fails gate_within_status_eta2 (line 187), which with the other hard gates writes islets_core_senior.rejected.h5ad and exits 1.  
Checked: No measurement anywhere: this file states only the ceiling, and the later lane confirms the direction of inheritance — trajectory_islet.py:85-87 (2026-09-14) says 'The legacy lane used the same value for the same reason: above it, the axis is a donor axis', i.e. it took 0.15 from here, and here it was introduced bare in 61f2006.

**`min_islets (within_status_batch_eta2 keyword default)` = `3`** — scripts/senior/compute_trajectory_senior.py:72 — *ARBITRARY*  
Controls: A donor contributes to its status's eta^2 only if it has at least this many islets (line 84); raising it drops sparse donors out of the batch metric.  
Checked: Confirmed by git log -L 72,72: min_islets=3 entered bare in 61f2006 and survived the 6e96cef signature change untouched; the docstring states the rule only, and no commit message, comment or sibling lane gives a number behind 3.

**`max_iter_harmony (hm.run_harmony)` = `20`** — scripts/senior/compute_trajectory_senior.py:98 (also :326, :348) — *INHERITED_DEFAULT*  
Controls: Iteration cap for every Harmony correction here — the sweep latent, the combined latent and the Leiden latent; convergence behaviour therefore differs from stock.  
Checked: Confirmed non-default against both runtime envs (harmonypy 0.2.0 in scvi-env and 2.0.0 in scimap both sign max_iter_harmony=10), and confirmed inherited from reaggregate_islets.py:144 and :498-499, where 20 likewise appears with no basis; note the repo .venv has no harmonypy at all, so this call cannot be exercised there.

