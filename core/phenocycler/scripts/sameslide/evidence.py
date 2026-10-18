"""Frozen, protein-independent RNA evidence for the four-slide development lane."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import poisson, chi2

ALPHA = 0.003  # nominal per-gene test, not a panel-wide FDR guarantee
Y_CLEAN = ('TBL1Y', 'USP9Y', 'EIF1AY')
Y_CROSSHYB = ('DDX3Y', 'KDM5D')
SEX_GENES = frozenset((*Y_CLEAN, *Y_CROSSHYB, 'TSPY1', 'XIST', 'TSIX'))
# Recovered verbatim from the established pln_control_sets.py curation.
CLEAN_NEG = (
    'ABCA4','ARR3','CNGB1','CNGB3','GUCA1A','GUCY2D','IMPG1','NR2E3','OPN1LW','PDE6H',
    'PRPH2','RPE65','VSX2','DSG1','DSG3','TGM3','COL2A1','PHEX','TNNI3','TNNT2','NPPA',
    'NPPB','TBX5','AKAP4','BOLL','DDX25','DPY19L2','MAGEA1','MAGEA10','MAGEA12','MAGEA3',
    'MAGEA4','MEIOB','PIWIL1','PRM1','SPATA16','SPATA22','SPATA4','STK31','SYCP1','TEX14',
    'DAZL','GDF9','KHDC3L','NLRP5','NLRP7','TLE6',
)


def validate_counts(counts):
    x = sparse.csr_matrix(counts, copy=True)
    x.sum_duplicates(); x.eliminate_zeros(); x.sort_indices()
    if not np.isfinite(x.data).all() or (x.data < 0).any() or (x.data != np.floor(x.data)).any():
        raise ValueError('raw counts must be nonnegative finite integers')
    return x


def detected(counts, rates, alpha=ALPHA):
    """Stable P(Poisson(lambda) >= observed), evaluated in bounded row chunks."""
    x = validate_counts(counts)
    rates = np.asarray(rates, float)
    if rates.shape != (x.shape[0],) or not np.isfinite(rates).all() or (rates < 0).any():
        raise ValueError('background rates must match cells and be finite nonnegative')
    if not 0 < alpha < 1:
        raise ValueError('invalid nominal threshold')
    out = x.copy().astype(bool)
    for start in range(0, x.shape[0], 10000):
        end = min(start + 10000, x.shape[0]); lo, hi = x.indptr[start], x.indptr[end]
        row = np.repeat(rates[start:end], np.diff(x.indptr[start:end + 1]))
        out.data[lo:hi] = poisson.sf(x.data[lo:hi] - 1, row) <= alpha
    out.eliminate_zeros()
    return out


def qc_table(counts, cells, rates, n_panel_genes):
    x = validate_counts(counts)
    rates = np.asarray(rates, float)
    if len(cells) != x.shape[0] or rates.shape != (len(cells),):
        raise ValueError('QC axes differ')
    if not np.isfinite(rates).all() or (rates < 0).any() or n_panel_genes < 1:
        raise ValueError('invalid background evidence')
    area = cells.cell_area.to_numpy(float); nuc = cells.nucleus_area.to_numpy(float)
    depth = np.asarray(x.sum(axis=1)).ravel(); genes = np.diff(x.indptr)
    valid_area = np.isfinite(area) & (area > 0)
    has_nucleus = np.isfinite(nuc) & (nuc > 0)
    biological = valid_area & has_nucleus
    if not biological.any():
        raise ValueError('no cells with a measured nucleus and valid area')
    density = np.divide(depth, area, out=np.full(len(area), np.nan), where=valid_area)
    cap = float(np.quantile(density[biological], .98))
    expected = rates * n_panel_genes
    signed = depth - expected
    reasons = np.full(len(cells), '', object)
    rules = {'no_nucleus': ~has_nucleus, 'invalid_cell_area': ~valid_area,
             'fewer_than_20_genes': genes < 20, 'adjusted_transcripts_below_50': signed < 50,
             'density_above_p98': density > cap}
    for name, mask in rules.items():
        reasons[mask] = np.char.add(reasons[mask].astype(str), name + ';')
    table = cells.copy()
    table['vendor_total_counts'] = cells.total_counts.to_numpy()
    table['raw_gene_counts'] = depth; table['raw_detected_genes'] = genes
    table['gdna_lambda'] = rates; table['expected_background'] = expected
    table['adjusted_qc_total'] = signed; table['count_density'] = density
    table['background_dominated'] = expected >= .5 * depth
    table['background_exceeds_observed'] = expected >= depth
    table['biological_eligible'] = biological
    table['exclusion_reasons'] = reasons
    table['qc_eligible'] = reasons == ''
    return table, {'min_genes': 20, 'min_adjusted_counts': 50, 'density_quantile': .98,
                   'density_cap': cap, 'nucleus_required': True,
                   'n_cells': len(cells), 'n_eligible': int((reasons == '').sum()),
                   'exclusion_counts': {k: int(v.sum()) for k, v in rules.items()}}


def wilson(k, n):
    if not n:
        return (np.nan, np.nan)
    z = 1.959963984540054; p = k / n; den = 1 + z*z/n
    center = (p + z*z/(2*n))/den
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
    return max(0., center-half), min(1., center+half)


def group_summary(counts, support, rates, genes):
    """Poisson count CI shifted by expected background; conditional on fitted rates."""
    x = sparse.csr_matrix(counts); n = x.shape[0]
    observed = np.asarray(x.sum(axis=0)).ravel().astype(np.int64)
    expected = float(np.asarray(rates).sum())
    low = np.where(observed > 0, .5*chi2.ppf(.025, 2*observed), 0.)
    high = .5*chi2.ppf(.975, 2*(observed + 1))
    raw = np.asarray((x > 0).sum(axis=0)).ravel()
    corr = np.asarray(support.sum(axis=0)).ravel()
    bg_detection = float(np.mean(-np.expm1(-np.asarray(rates)))) if n else np.nan
    return pd.DataFrame({'gene': genes, 'n_cells': n, 'observed_counts': observed,
        'expected_background': expected, 'signed_excess': observed - expected,
        'excess_ci95_low_conditional': low - expected, 'excess_ci95_high_conditional': high - expected,
        'raw_detected_cells': raw, 'corrected_detected_cells': corr,
        'raw_detection_rate': raw/max(n,1), 'corrected_detection_rate': corr/max(n,1),
        'below_section_background': raw/max(n,1) <= bg_detection})


def control_report(counts, support, genes, cells, sex, tissue):
    if sex not in ('F', 'M'):
        raise ValueError('recorded sex must be verified before control assignment')
    raw = sparse.csr_matrix(counts) > 0; support = sparse.csr_matrix(support)
    pos = {g: i for i, g in enumerate(genes)}
    negatives = [g for g in CLEAN_NEG if tissue == 'pancreas' or g != 'NLRP7']
    sets = {'Y_clean': Y_CLEAN, 'Y_crosshyb': Y_CROSSHYB,
            'X_diagnostic': ('XIST','TSIX'), 'curated_negative': negatives}
    strata = {'section': np.ones(len(cells), bool)}
    for column in ('raw_gene_counts','gdna_lambda'):
        # Quantile boundaries, never split equal values by row order.
        bins = pd.qcut(cells[column], 5, labels=False, duplicates='drop')
        for b in sorted(bins.dropna().unique()):
            strata[f'{column}_q{int(b)+1}'] = bins.to_numpy() == b
    rows = []
    for group, names in sets.items():
        available = [g for g in names if g in pos]
        for stratum, mask in strata.items():
            n = int(mask.sum())
            for name, members in [(g, [g] if g in pos else []) for g in names] + [('ANY', available)]:
                ids = [pos[g] for g in members]
                k = int(np.asarray(support[mask][:,ids].sum(axis=1)).ravel().astype(bool).sum()) if ids else 0
                kr = int(np.asarray(raw[mask][:,ids].sum(axis=1)).ravel().astype(bool).sum()) if ids else 0
                low, high = wilson(k,n)
                female_test = sex == 'F' and group.startswith('Y_') and name != 'ANY'
                status = ('missing_probe' if not ids else
                          ('pass' if k/n <= ALPHA else 'fail') if female_test and n else
                          'female_Y_negative_not_applicable' if sex == 'M' and group.startswith('Y_') else 'diagnostic')
                rows.append({'set':group,'probe':name,'stratum':stratum,'n_cells':n,
                    'available_probes':';'.join(members),'n_detected':k if ids else np.nan,
                    'raw_detection_rate':kr/n if n and ids else np.nan,
                    'corrected_detection_rate':k/n if n and ids else np.nan,
                    'ci95_low':low if ids else np.nan,'ci95_high':high if ids else np.nan,
                    'target':ALPHA if female_test else np.nan,'status':status})
    return pd.DataFrame(rows)
