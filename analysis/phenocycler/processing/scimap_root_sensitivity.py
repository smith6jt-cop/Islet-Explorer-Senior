"""Compare five broad pLN root calls with SCIMAP 2.3.5 on donor- and status-stratified cells.

Use the original scimap environment. Its continuous scores receive the active
dedicated-section gates, with the floor, cutoff and top mapped to 0, .5 and 1.
This sensitivity does not apply the archived manual 6551 gates to other slides.
"""
import argparse
from pathlib import Path
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from multiplex5000.tissues.lymph_node import load_rule_source
from scimap_all_cells import DRAW, SEED, root_name, scorer


def section(donor, roots, markers, phenotype_cells, *, calls_dir, dedicated_root, gate_table, disease, source_roots):
    entry = donor + "__pln_1"
    mask = pd.read_parquet(calls_dir / f"{donor}.parquet",
                           columns=["analysis_eligible", "candidate_state", "candidate_phenotype"])
    mask = mask.loc[mask.analysis_eligible].copy()
    states = mask.groupby("candidate_state", observed=True).size()
    chosen = [v.sample(min(len(v), DRAW), random_state=SEED + int(donor)) for _, v in
              mask.groupby("candidate_state", observed=True, sort=True)]
    sample = pd.concat(chosen)
    branch = sample.candidate_phenotype.astype("string").map(source_roots)
    sample["source_root"] = branch.fillna(sample.candidate_state.astype("string"))
    thresholds = gate_table[gate_table.protein_entry.eq(entry)].set_index("marker").loc[markers]
    raw = pq.ParquetFile(dedicated_root / "03_protein_measurement" / "cells" / f"{entry}.parquet").read(
        columns=["object_id"] + ["marker:" + m for m in markers]).to_pandas().set_index("object_id")
    values = raw.loc[sample.index].to_numpy(dtype=float) * thresholds.shift_scale.to_numpy()
    valid = np.isfinite(values) & (values >= 0)
    scaled = np.column_stack([np.interp(values[:, index], [t.floor, t.cutoff, t.top], [0, .5, 1])
                              for index,t in enumerate(thresholds.itertuples())])
    if not np.array_equal((scaled >= .5)[valid], (values >= thresholds.cutoff.to_numpy())[valid]):
        raise ValueError(f"{donor}: alternate scale changed gate classification")
    complete = valid.all(axis=1)
    sample["scimap_root"] = "Unassessable root marker"
    obj = ad.AnnData(scaled[complete], obs=pd.DataFrame(index=sample.index[complete]),
                     var=pd.DataFrame(index=markers))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=FutureWarning)
        warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
        phenotype_cells(obj, roots.copy(), gate=.5, label="sensitivity", verbose=False)
    result = obj.obs.sensitivity.astype("string")
    if not result.isin(set(roots.Phenotype) | {"Unknown"}).all():
        raise ValueError(f"{donor}: unexpected SCIMAP broad-root output")
    sample.loc[result.index, "scimap_root"] = result
    out = sample.groupby(["candidate_state", "source_root", "scimap_root"], observed=True,
                         sort=True).size().rename("sampled_cells").reset_index()
    drawn = sample.groupby("candidate_state", observed=True).size()
    out["stratum_population"] = out.candidate_state.astype("string").map(states)
    out["stratum_sampled"] = out.candidate_state.astype("string").map(drawn)
    out["estimated_cells"] = out.sampled_cells * out.stratum_population / out.stratum_sampled
    return out.assign(donor=donor, disease=disease)


def main(data_dir, calls_dir, dedicated_root, manifest):
    phenotype_cells = scorer()
    source = load_rule_source(manifest, dedicated=True)
    rules = [rule for rule in source.rules if rule.phenotype in source.candidate_lineage_roots]
    sheet = pd.DataFrame([{"CellType":rule.parent, "Phenotype":rule.phenotype, **dict(rule.constraints)}
                          for rule in rules])
    proteins = list(sheet.columns[2:])
    donors = sorted(f.stem for f in calls_dir.glob("*.parquet"))
    gate_table = pd.read_csv(dedicated_root / "06_typing/protein_thresholds.csv")
    source_roots = {r.phenotype: root_name(source,r.phenotype) for r in source.rules}
    clinical = pd.read_parquet(data_dir/"donor_counts.parquet", columns=["tissue","donor","disease"])
    disease = clinical[clinical.tissue.eq("pLN")].drop_duplicates().set_index("donor").disease
    all_sections = pd.concat([section(d,sheet,proteins,phenotype_cells,calls_dir=calls_dir,dedicated_root=dedicated_root,
                                      gate_table=gate_table,disease=disease.at[d],source_roots=source_roots)
                              for d in donors],ignore_index=True)
    target = data_dir / "pln_no_match_scimap_sensitivity.parquet"
    tmp = target.with_suffix(".parquet.tmp")
    all_sections.to_parquet(tmp, index=False)
    tmp.replace(target)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["data-dir","calls-dir","dedicated-root","manifest"]:
        parser.add_argument("--"+name,type=Path,required=True)
    main(**vars(parser.parse_args()))
