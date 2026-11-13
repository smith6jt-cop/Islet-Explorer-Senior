"""Fit a separate low-inflammation-rooted PhenoCycler composition trajectory for each donor."""
from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from threadpoolctl import threadpool_limits


def digest(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def fit_donor(data):
    """Use native composition inputs. Density selects the root only after graph construction."""
    data = data[np.argsort(data.obs_names)].copy()
    o = data.obs
    if not o.index.is_unique or not np.array_equal(o.index.astype(str), o.islet_id.astype(str)):
        raise ValueError('Native islet UUID mismatch')
    donor = o.donor.astype(str).unique().item()
    disease = o.disease.astype(str).unique().item()
    values = np.asarray(data.X, dtype=float)
    counts = np.asarray(data.layers['counts'], dtype=float)
    denominator = o.n_core_endocrine.to_numpy(float)
    with np.errstate(invalid='ignore', divide='ignore'):
        expected = np.log1p(counts / np.where(denominator > 0, denominator, np.nan)[:, None])
    np.testing.assert_allclose(values, expected, equal_nan=True)
    density = o.combined_immune_per_mm2
    direct = (1e6 * (o.n_core_immune + o.n_peri_immune) / (o.core_area_um2 + o.peri_area_um2))
    supported = o.n_core.gt(0) & o.n_peri.gt(0) & o.core_area_um2.gt(0) & o.peri_area_um2.gt(0)
    np.testing.assert_allclose(density, direct.where(supported), equal_nan=True)
    out = o[['islet_id', 'donor', 'disease', 'combined_immune_per_mm2', 'composition_dpt']].copy()
    out['graph_exclusion'] = ''
    conditions = [
        (~(o.geometry_valid & o.n_peri.gt(0) & o.core_pancreatic_tissue_fraction.ge(1 - 1e-6)), 'FEATURE_INPUT_UNAVAILABLE'),
        (~np.isfinite(values).all(axis=1), 'NONFINITE_COMPOSITION'),
        (o.n_core.lt(10), 'BELOW_10_ELIGIBLE_CORE_CELLS'),
        (o.n_core_endocrine.le(0), 'NO_ELIGIBLE_CORE_ENDOCRINE_CELL'),
        (o.original_object_class.ne('Islet'), 'NOT_NATIVE_ISLET'),
    ]
    for mask, label in conditions:
        out.loc[mask, 'graph_exclusion'] = out.loc[mask, 'graph_exclusion'].map(
            lambda reason: reason + ';' + label if reason else label)
    out['graph_cluster'] = pd.Series(index=out.index, dtype='string')
    for column in ('graph_x', 'graph_y', 'graph_dpt'):
        out[column] = np.nan
    out['graph_root'] = False
    record = dict(donor=donor, disease=disease, native_islets=len(o), n_fitted=0, n_components=None,
                  n_outside_component=0, root_islet_id=None, root_density=None, status=None)
    def unavailable(reason):
        out.loc[out.graph_exclusion.eq(''), 'graph_exclusion'] = reason
        record['status'] = reason
        return out.reset_index(drop=True), [], [], record

    keep = np.flatnonzero(out.graph_exclusion.eq(''))
    if len(keep) < 50:
        return unavailable('FEWER_THAN_50_ELIGIBLE_ISLETS')
    columns = np.ptp(values[keep], axis=0) > 1e-12
    if columns.sum() < 3:
        return unavailable('FEWER_THAN_3_VARIABLE_FEATURES')
    features = values[keep][:, columns]
    zero = np.abs(features).sum(axis=1) <= 1e-12
    out.loc[o.index[keep[zero]], 'graph_exclusion'] = 'ZERO_VARIABLE_FEATURE_VECTOR'
    keep = keep[~zero]
    if len(keep) < 50:
        return unavailable('FEWER_THAN_50_NONZERO_ELIGIBLE_ISLETS')
    unique, inverse = np.unique(np.round(features[~zero], 9), axis=0, return_inverse=True)
    inverse = inverse.ravel()
    if len(unique) < 12:
        return unavailable('FEWER_THAN_12_DISTINCT_PROFILES')
    graph = ad.AnnData(unique)
    sc.pp.neighbors(graph, n_neighbors=min(15, len(unique)-1), metric='cosine', use_rep='X', random_state=0)
    ncomponents, components = connected_components(graph.obsp['connectivities'], directed=False)
    primary = np.flatnonzero(components == np.argmax(np.bincount(components)))
    inside = np.isin(inverse, primary)
    record.update(n_components=int(ncomponents), n_outside_component=int((~inside).sum()),
                  n_distinct_profiles=len(unique), features_used=data.var_names[columns].tolist(),
                  component_profile_counts=np.bincount(components).tolist())
    out.loc[o.index[keep[~inside]], 'graph_exclusion'] = 'OUTSIDE_LARGEST_CONNECTED_COMPONENT'
    if len(primary) < 12 or inside.sum() < 50:
        return unavailable('PRIMARY_COMPONENT_BELOW_50_ISLETS_OR_12_PROFILES')
    graph = graph[primary].copy()
    islets = o.iloc[keep[inside]]
    mapping = np.searchsorted(primary, inverse[inside])
    candidates = islets.dropna(subset=['combined_immune_per_mm2', 'root_reference_distance'])
    if candidates.empty:
        return unavailable('NO_ASSESSED_ROOT')
    root = candidates.sort_values(['combined_immune_per_mm2', 'root_reference_distance', 'islet_id']).iloc[0]
    profile = int(mapping[islets.index.get_loc(root.name)])
    graph.uns['iroot'] = profile
    sc.tl.diffmap(graph, n_comps=10, random_state=0)
    sc.tl.dpt(graph, n_dcs=10, n_branchings=0)
    sc.tl.leiden(graph, resolution=1, flavor='igraph', n_iterations=2, directed=False, random_state=0)
    sc.tl.paga(graph, groups='leiden')
    sc.pl.paga(graph, plot=False, random_state=0)
    sc.tl.draw_graph(graph, layout='fr', init_pos='paga', random_state=0)
    dpt = graph.obs.dpt_pseudotime.to_numpy(float)[mapping]
    finite = np.isfinite(dpt)
    index = islets.index[finite]
    out.loc[islets.index[~finite], 'graph_exclusion'] = 'NONFINITE_DIFFUSION_DISTANCE'
    out.loc[index, 'graph_dpt'] = dpt[finite]
    out.loc[index, ['graph_x','graph_y']] = graph.obsm['X_draw_graph_fr'][mapping[finite]]
    out.loc[index, 'graph_cluster'] = graph.obs.leiden.astype(str).to_numpy()[mapping[finite]]
    out.loc[root.name, 'graph_root'] = True
    if out.loc[root.name, 'graph_dpt'] != 0 or not len(index):
        raise ValueError(f'{donor}: invalid diffusion root')
    nodes = out.loc[index].groupby('graph_cluster', observed=True).agg(
        graph_x=('graph_x', 'median'), graph_y=('graph_y', 'median'),
        n_islets=('islet_id', 'size'), median_dpt=('graph_dpt', 'median'),
        median_density=('combined_immune_per_mm2', 'median')).reset_index()
    nodes['donor'] = donor
    nodes['is_root_cluster'] = nodes.graph_cluster.eq(out.loc[root.name, 'graph_cluster'])
    groups = graph.obs.leiden.cat.categories.astype(str)
    connectivity = sparse.triu(graph.uns['paga']['connectivities'], k=1, format='coo')
    tree = graph.uns['paga']['connectivities_tree'].tocsr()
    tree = tree.maximum(tree.T)
    edges = [dict(donor=donor, source=groups[i], target=groups[j], connectivity=float(value),
                  in_tree=bool(tree[i,j] > 0)) for i,j,value in zip(connectivity.row, connectivity.col, connectivity.data)
             if value > 0]
    g = graph.obsp['connectivities'].tocsr()
    h = hashlib.sha256()
    for array in (primary.astype('<i8'), g.indices.astype('<i8'), g.indptr.astype('<i8'), g.data.astype('<f8')):
        h.update(array.tobytes())
    record.update(status='fitted', n_fitted=len(index), n_clusters=len(nodes), n_paga_edges=len(edges),
                  root_islet_id=str(root.name), root_density=float(root.combined_immune_per_mm2),
                  root_reference_distance=float(root.root_reference_distance),
                  root_profile_islets=int((mapping == profile).sum()), n_minimum_density_roots=int(
                      candidates.combined_immune_per_mm2.eq(root.combined_immune_per_mm2).sum()),
                  graph_sha256=h.hexdigest())
    return out.reset_index(drop=True), nodes.to_dict('records'), edges, record


def build(path, output):
    data = ad.read_h5ad(path)
    if not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ValueError('Duplicated UUID or composition feature')
    if not set(data.var.region.astype(str)).issubset({'core','peri'}):
        raise ValueError('Unexpected region in composition inputs')
    for name in data.var_names:
        if not name.startswith(('core_','peri_')):
            raise ValueError(f'Unexpected composition feature {name}')
    points, nodes, edges, records = [], [], [], []
    with threadpool_limits(limits=1):
        for donor, indices in data.obs.groupby('donor', observed=True).indices.items():
            rows, vertices, links, record = fit_donor(data[indices])
            points.append(rows); nodes.extend(vertices); edges.extend(links); records.append(record)
            print(donor, record['status'], record['n_fitted'], flush=True)
    output.mkdir(parents=True, exist_ok=True)
    dfs = {'points':pd.concat(points, ignore_index=True), 'nodes':pd.DataFrame(nodes),
           'edges':pd.DataFrame(edges, columns=['donor','source','target','connectivity','in_tree'])}
    with TemporaryDirectory(dir=output) as temp:
        base = Path(temp)
        for name, frame in dfs.items():
            frame.to_parquet(base/f'donor_islet_trajectory_{name}.parquet', index=False)
        method = dict(version='within-donor-low-density-root-dpt-v1', source_file=Path(path).name,
                      input_sha256=digest(path), producer_sha256=digest(__file__),
                      software={p:version(p) for p in ['scanpy','anndata','igraph','numpy','scipy']},
                      graph='Each donor independently: log1p(region phenotype count / core endocrine count), variable columns; distinct profiles rounded to 9 decimals; cosine 15-NN, largest connected component; 10 diffusion components, DPT with no branch assignment.',
                      root='Minimum measured combined core + 20 µm peri-islet immune density among fitted islets; break ties by distance to the original nonimmune non-diabetic composition reference (omits the current donor), then native UUID.',
                      coordinates='Within-donor PAGA-initialized Fruchterman–Reingold layout of distinct profiles, Leiden resolution 1. Node positions are median member coordinates. Edges are actual Scanpy PAGA connectivities; in_tree identifies the undirected Scanpy connectivity tree. Cluster labels and layout axes are local to the donor.',
                      interpretation='Cross-sectional distance through islet compositions from a measured low-density root. Immune counts enter composition features, so association with density is descriptive. Diffusion distance does not imply progressively increasing density, branching in time, or observed disease progression.',
                      output_sha256={p.name:digest(p) for p in sorted(base.glob('*.parquet'))}, donors=records)
        receipt = base/'donor_islet_trajectory_method.json'
        receipt.write_text(json.dumps(method, indent=2, allow_nan=False)+'\n')
        for p in sorted(base.glob('*.parquet')):
            p.replace(output/p.name)
        receipt.replace(output/receipt.name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--islets', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    build(args.islets, args.output)
