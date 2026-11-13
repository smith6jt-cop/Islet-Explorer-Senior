"""Count measured positive defining markers in formerly unmatched cells receiving each score-based root."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from multiplex5000.tissues.lymph_node import load_rule_source


def build(data_dir, dedicated_root, manifest):
    source = load_rule_source(manifest,dedicated=True)
    thresholds = pd.read_csv(dedicated_root/"06_typing/protein_thresholds.csv")
    counts = pd.read_parquet(data_dir/"pln_scimap_candidate_counts.parquet")
    counts = counts[counts.candidate_state.eq("no_match")]
    records = []
    for donor, cohort in counts.groupby("donor",observed=True):
        cells = pd.read_parquet(data_dir/"pln_scimap_candidate_cells.parquet",
                                filters=[("donor","=",donor),("candidate_state","=","no_match")],
                                columns=["object_id","scimap_candidate_root"]).set_index("object_id")
        if not cells.index.is_unique or len(cells) != cohort.cells.sum():
            raise ValueError(f"{donor}: incomplete strict-to-SCIMAP membership")
        entry = donor+"__pln_1"
        gates = thresholds[thresholds.protein_entry.eq(entry)].set_index("marker")
        panels = {root:[(m,op) for m,op in source.by_phenotype[root].constraints if op!="allneg"]
                  for root in source.candidate_lineage_roots}
        markers = sorted({m for rules in panels.values() for m,_ in rules})
        raw = pq.ParquetFile(dedicated_root/"03_protein_measurement"/"cells"/f"{entry}.parquet").read(
            columns=["object_id",*["marker:"+m for m in markers]],use_threads=False).to_pandas().set_index("object_id")
        if not raw.index.is_unique:
            raise ValueError(f"{donor}: repeated native ID")
        values = raw.loc[cells.index].rename(columns=lambda c:c.removeprefix("marker:"))
        for root, group in panels.items():
            pick = cells.scimap_candidate_root.eq(root)
            if pick.sum() != cohort.loc[cohort.scimap_candidate_root.eq(root),"cells"].sum():
                raise ValueError(f"{donor}: inconsistent selected root membership")
            for marker, role in group:
                readings = values.loc[pick,marker]
                valid = np.isfinite(readings) & readings.ge(0)
                threshold = gates.loc[marker]
                positive = valid & (readings*threshold.shift_scale).ge(threshold.cutoff)
                records.append(dict(donor=donor,scimap_candidate_root=root,marker=marker,rule_term=role,
                                    cells=int(pick.sum()),measured=int(valid.sum()),positive=int(positive.sum())))
    output = data_dir/"pln_scimap_no_match_root_markers.parquet"
    stage = output.with_suffix(".building")
    pd.DataFrame(records).to_parquet(stage,index=False)
    stage.replace(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["data-dir","dedicated-root","manifest"]:
        parser.add_argument("--"+name,type=Path,required=True)
    build(**vars(parser.parse_args()))
