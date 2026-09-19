# Statistics Tab — Detailed Reference

Detailed reference for the Statistics tab. Linked from `apps/islet/CLAUDE.md`.

> **SCOPE, added 2026-09-18.** The *machinery* below — donor-level aggregation, the mixed-effects
> model, the shared per-bin coverage rule, the fixed bin width and the fixed alpha — is
> domain-independent and current. **Every NUMBER in it is the LEGACY CODEX DOMAIN's**
> (`config/islet_default.yml`, 15 donors, 5 per group, 5,214 islets), because that was the only
> dataset when this was written. The PhenoCycler run cohort (`config/run20_default.yml`) has 20
> donors — 8 ND / 3 Aab+ / 9 T1D, 18,384 islets, measured 2026-09-18 — so read `N=15` and
> `N=5,214` as examples of the shape, never as the app's cohort. The tab itself derives both from
> the loaded data; no cohort size is written into the UI.

## Statistics Tab (Phase 6+14+16, Feb-Mar 2026)

### Pseudoreplication Fix (Phase 16)

The app has 5,214 islets from only **15 donors** (5 per group). Islets within a donor are correlated. Phase 16 fixes this textbook pseudoreplication by:

1. **Donor-level aggregation**: All tests use `aggregate_to_donor()` to compute per-donor means (N=15), not islet-level (N=5,214). Produces realistic p-values (0.01-0.05 range, not 1e-50).
2. **Mixed-effects sensitivity**: `lmer_test_donor()` runs `lmerTest::lmer(value ~ donor_status + (1 | case_id))` on islet-level data with donor as random intercept. ICC quantifies donor-level clustering.
3. **Min-cells filter**: Available on all 4 analysis tabs — Plot, Trajectory, Statistics, Spatial. All use `total_cells_core + total_cells_peri >= threshold` (Trajectory uses `total_cells` directly). Defaults to 1.
4. **Normality testing**: Shapiro-Wilk on donor-level means per group, with auto-suggestion for test type.
5. **No-binning option**: `stats_no_binning` checkbox hides Section 3 (size-dependent analysis).
6. **Bin slider fix**: Max increased from 75 to 150, step=5. Heatmap x-axis labels angled (-45°).
7. **Per-bin donor-level**: `per_bin_donor_anova()` and `per_bin_donor_kendall()` aggregate to donor means within each bin. Bins with <2 donors per group are skipped.
8. **Demographics redesign**: Age scatter shows all islets (not donor means) with Pearson r. Sex analysis uses box plot faceted by sex. Autoantibody analysis (Aab+ only): per-donor AAb profile table + feature by AAb count box plot. Use 'Sex' not 'Gender' in user-facing labels.

### Key Utilities (`utils_stats.R`)
- `aggregate_to_donor(rdf, agg_fn)` -- Groups by `Case ID` + `donor_status`, computes mean/median of `value`, preserves age/gender
- `lmer_test_donor(rdf)` -- Mixed-effects model, returns `list(p_value, icc, fit, error_msg)`. Uses `lmerTest::anova()` (Satterthwaite F-test). Returns diagnostic `error_msg` on failure (singular fit, convergence, missing packages).
- `per_bin_donor_anova(df, ...)` -- Donor-level per-bin ANOVA. Skips bins unless ≥2 donor-status groups each have ≥2 donors (shared coverage rule).
- `per_bin_donor_kendall(df, ...)` -- Donor-level per-bin Kendall tau. Same coverage rule as `per_bin_donor_anova()` so the heatmap and trend plot include the same bins (May 2026 harmonization).
- `code_donor_status(x)` -- Maps `ND/Aab+/T1D` → `0/1/2` for Kendall trend coding. Single source of truth, used by both `per_bin_kendall()` and `per_bin_donor_kendall()`.
- `normality_tests(ddf)` -- Shapiro-Wilk per group, returns data.frame with `W`, `p_value`, `is_normal`
- `cohens_d(x, y)` -- returns `list(d, ci_lo, ci_hi)` (NOT `ci_lower`/`ci_upper`)
- `eta_squared(fit)` -- from `anova()` output
- `pairwise_wilcox(df, group_col, value_col)` -- returns data.frame with columns `group1`, `group2`, `p_value` (NOT `p.adj`/`statistic`)
- `flag_outliers(df, value_col, group_col, method, threshold, transform, exclude_zeros, enabled)` -- Single source of truth for outlier flagging. Appends `is_outlier` (logical) and `outlier_score` (z or fence-distance) columns. Attributes carry method/threshold/enabled for downstream summaries.
- `summarize_outlier_filter(df)` -- Compact summary list for badges: `n_total`, `n_excluded`, `n_kept`, `method`, `threshold`, `enabled`, `label`.
- `outlier_title_suffix(summary)` -- Formats the `" — N outliers excluded (label)"` string for plot titles. Returns `""` when no exclusion.
- `GLOBAL_OUTLIER_TOOLTIP` -- String constant used as `title=` attribute on every outlier control across tabs.
- `strip_plotly_color_leak(gg)` -- Scrubs each trace's `text`/`hovertemplate`/`hovertext` to remove `colour:\s*black` (and `color:\s*black`) lines, plus collapses adjacent duplicate `donor_status:\s*<value>` lines. Apply after every `ggplotly()` call where any geom uses `color = "black"` as a parameter (currently: Plot scatter, Plot distribution panel, Trajectory scatter). Defensive — no-op when no leak is present.

### Unified Outlier Handling (May 2026 — Statistics decoupled May 14; defaults OFF + dynamic z May 15)
- One non-namespaced checkbox + slider pair per tab (`sidebar_remove_outliers` + `sidebar_outlier_z`, `traj_remove_outliers` + `traj_outlier_z`, `spatial_remove_outliers` + `spatial_outlier_z`) — synced bidirectionally via six `observeEvent`s in `app.R` (three for the checkboxes, three for the sliders). Reactives `remove_outliers_global()` and `outlier_z_global()` are injected into `plot_server`, `trajectory_server`, and `spatial_server` as the `remove_outliers` and `outlier_threshold` arguments.
- **Statistics tab is independent**: namespaced controls `stats_remove_outliers` and `stats_outlier_z` live inside `mod_statistics_ui.R` Section 1 and are read via `input$…` inside `stats_run()`. They do NOT sync with the Plot/Trajectory/Spatial toggles. The `statistics_server()` function still accepts the `remove_outliers` argument for backward compatibility but ignores it.
- Defaults: **OFF** for all four checkboxes (raw distribution renders out of the box; users opt in to outlier exclusion). z-threshold slider default 3.0, range 0.5–10, step 0.5, on every tab. Toggling the Plot sidebar propagates to Trajectory and Spatial but leaves Statistics alone (and vice versa); same for the slider.
- Plot, Trajectory, Statistics, Spatial Card C: `flag_outliers(method = "zscore", threshold = outlier_threshold())` per donor status — threshold pulled from the live slider. Card C transforms via `log1p` to handle the long-tailed distance distribution.
- **Robust two-sided z (May 15)**: the `method = "zscore"` branch uses the Iglewicz–Hoaglin modified z-score, `0.6745·(x − median)/MAD`, instead of the classical mean/SD z. Because every metric in this app is non-negative and right-skewed, mean/SD pushed the achievable lower z toward 0 and the rule only flagged the upper tail in practice. Median/MAD is symmetric and genuinely two-sided. Fallback: when `MAD = 0` (zero-inflated groups), the function reverts per-group to classical mean/SD. The threshold slider scale is unchanged; the 0.6745 rescale keeps `|z| > 3` roughly equivalent to the old classical cutoff. The active rule is exposed on `attr(df, "outlier_statistic")` (`"modified_z"`, `"classical_z"`, or `"modified_z (classical_z fallback used for some groups)"`).
- Spatial Card A: `flag_outliers(method = "iqr", threshold = 1.5, exclude_zeros = TRUE)` per donor status. Keeps the Tukey logic because phenotype enrichment is zero-inflated; intentionally NOT routed through the slider (different scale).
- Visual surfacing: every affected plot shows a **count badge** in its title or subtitle and renders excluded points as **shape=4 grey × ghost markers** (Plot, Trajectory, Spatial Card C). Stats forest plot title shows the badge; the Statistics outlier checkbox carries a small caption beneath it reporting `N excluded last Run` (or a yellow "Filter OFF" notice). Statistics Methods Reference paragraph names the local rule and cites the exact `n_excluded` for the current Run. Download CSV header records `Outlier rule: …` and the count.
- Spatial Card B's "Clip z > 5" checkbox is a **display saturation only** (caps the colour scale) — left alone, NOT an outlier control.
- Trajectory bug fix: previously discarded outliers silently before plotting and logged only to console; now drawn as ghosts and counted in a badge above the scatter.
- **Spatial per-card outlier tables** (May 13 follow-up): Cards A and C now expose `infiltration_show_outlier_table` and `distance_show_outlier_table` checkboxes that surface the list of excluded islets in a yellow warning panel below each plot. Reactive state lives in module-scope `infiltration_outliers` / `distance_outliers` `reactiveVal`s, populated inside each renderPlotly after the `flag_outliers()` call. Renders three states: table, "no outliers detected" blue note, or filter-OFF amber label. Matches the Plot/Trajectory `show_outlier_table` UX exactly.
- **Plot distribution panel independence** (May 13 follow-up): `output$dist` now reads `raw_df_base()` (not `raw_df()`) and applies its own `exclude_zero_dist` filter, so left and right panels' "Exclude zero values" controls no longer cascade. The dist panel also re-runs `flag_outliers()` on its own subset and renders ghost × markers + a title-suffix badge, fixing the gap where outliers only appeared on the left panel. Both panels' `position_jitter()` calls now use explicit `seed = 17` (normal points) / `seed = 31` (ghost ×) so visuals stay stable across reactive re-renders.
- **`strip_plotly_color_leak` post-process** (May 13 follow-up): added because hard-coded `color = "black"` parameters on `geom_point(shape=21, …, color="black")` (Plot summary points) and `geom_smooth(…, color="black")` (Trajectory overall trend) were surfacing as a literal `colour: black` line in plotly hover tooltips, alongside a duplicate `donor_status` entry from ggplotly emitting both `fill` and inherited `colour` aesthetics. Applied at `mod_plot_server.R:684`, `mod_plot_server.R:891`, and `mod_trajectory_server.R:662`.

### R Package Dependencies for Phase 16
- `lmerTest` (3.2-1) -- Mixed-effects model with Satterthwaite df
- `lme4` (1.1.37) -- `VarCorr()` for ICC computation

### Shared Sidebar Architecture
The Plot sidebar (mode, feature, region, donor status, AAb, age, gender filters) is visible on both the Plot and Statistics tabs. Achieved via `conditionalPanel` condition `"input.tabs == 'Plot' || input.tabs == 'Statistics'"` and matching JS `adjustLayout()` logic in `app.R`. The Statistics module consumes `plot_returns$raw_df` and `plot_returns$summary_df` directly -- no data duplication.

### Zero-Value Exclusion (May 2026)
- UI: `stats_exclude_zero` checkbox in Section 1 controls row (`mod_statistics_ui.R`). Default OFF.
- Server: `stats_run()` filters `value != 0` immediately after outlier flagging and before min-cells / diameter / aggregation. State and drop count (`zero_dropped`) propagate through the result list.
- Scope: applies to mixed-effects model, donor-level ANOVA / Kruskal, pairwise Cohen's d, per-bin tests, demographics regressions, and the covariate-adjusted model. Does NOT alter AUC (read from the Plot tab's `summary_df`) or the normality table (computed before filtering would change power conclusions for a comparison user's reference).
- Surfacing: overview banner shows an orange "Zero-valued excluded (N dropped)" badge; forest plot title appends `— N zero-valued dropped`; Methods Reference inserts a dedicated paragraph; CSV header records `Zero-value filter: ON/OFF; N dropped = …`.
- Use case: zero-inflated features (rare phenotypes, peri-islet enrichment z-scores) where many islets register exactly 0. Lets users sanity-check whether group differences survive on the non-zero subset.

### Sex-Stratified Tests — BH Correction (May 2026)
- Pairwise tests within each sex are run UNADJUSTED (`p.adjust.method = "none"`) and then BH-corrected across all 6 contrasts (3 pairs × 2 sexes) in `gender_strat$p_adj`. Replaces the prior "minimum raw p, no correction" reporting.
- `gender_plot` subtitle reports minimum BH-adjusted q per sex with a `(BH-corrected across 6 contrasts)` annotation. The `low_power` flag (n=2-3 donors per sex × group) remains as a separate badge.
- CSV export includes a `# Sex-Stratified Pairwise (BH-corrected across 6 contrasts: 3 pairs x 2 sexes)` section.

### 5-Section Narrative Layout

A "About this tab" banner above Section 1 describes the tab's purpose, primary test (mixed-effects with donor random intercept), and how to use it. Each section heading carries a `(?)` icon with a native-tooltip one-line description (helper: `section_heading(step, title, subtitle, tip = NULL)` in both `mod_statistics_ui.R` and `mod_statistics_server.R`).

1. **Configure Analysis** -- Overview banner, Run/CSV buttons, test type, alpha, outliers, min-cells filter, no-binning checkbox, **zero-exclusion checkbox**, bin width, diameter range, normality test results.
2. **Primary Results** -- Pseudoreplication info banner (ICC, mixed-effects p) + hypothesis table (col-4) + forest plot (col-4) + AUC (col-4). Cards equal-height via `display: flex` row.
3. **Size-Dependent Patterns** -- Wrapped in `conditionalPanel` for no-binning toggle. Donor-level stratified tests heatmap + trend analysis. Heatmap and trend plot share a single coverage rule (≥2 groups × ≥2 donors).
4. **Confounders & Deeper Analysis** -- Demographics: donor summary table, age scatter (islet-level, all points), sex box plot (faceted by sex, BH-corrected across 6 contrasts), autoantibody profile (Aab+ only), covariate model (donor-level)
5. **Methods Reference** -- Pseudoreplication, donor-level aggregation, ICC, mixed-effects, normality testing, zero-value exclusion state, harmonized per-bin rule.

### Multiple Testing Correction
- Per-bin ANOVA, Kendall τ, and pairwise tests are BH-corrected across bins via `p.adjust(method = "BH")`. Heatmap shows `-log10(q)` (corrected) with star annotations using corrected values.
- Donor-level pairwise contrasts are BH-corrected across the 3 pairs.
- Sex-stratified pairwise contrasts are BH-corrected across all 6 contrasts (3 pairs × 2 sexes).

### Plotly Legend Positioning
Forest plot uses `theme(legend.position = "none")` in ggplot (prevents ggplotly auto-placement), then `layout(margin = list(b = 60), legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.25))` in plotly. The `margin(b=60)` creates room below the x-axis title to prevent overlap. The trend plot uses `margin(b=100)` with legend at `y = -0.3`.

