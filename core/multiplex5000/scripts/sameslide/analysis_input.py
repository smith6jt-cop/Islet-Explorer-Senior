"""Common eligible-cell contract for four-slide downstream consumers."""
from pathlib import Path
import numpy as np
from scipy import sparse


def load_base(path, *, donor, roi, serial):
    import anndata as ad
    a=ad.read_h5ad(Path(path))
    if not a.obs_names.is_unique or not a.var_names.is_unique:raise ValueError('nonunique RNA axes')
    for col,value in [('donor_id',donor),('roi',roi),('sample',serial)]:
        if set(a.obs[col].astype(str))!={str(value)}:raise ValueError(f'section mismatch: {col}')
    for col in ('nucleus_area','cell_area','gdna_lambda'):
        if not np.isfinite(a.obs[col]).all():raise ValueError(f'nonfinite {col}')
    if (a.obs.nucleus_area<=0).any() or (a.obs.cell_area<=0).any() or not a.obs.qc_eligible.all():
        raise ValueError('RNA population is not nucleus-aware QC eligible')
    x=sparse.csr_matrix(a.layers['counts'])
    if not np.issubdtype(x.dtype,np.integer) or (x.data<0).any():raise ValueError('raw integer counts required')
    if not np.isfinite(a.obsm['spatial']).all():raise ValueError('invalid spatial coordinates')
    from evidence import SEX_GENES
    if any(bool(a.var.loc[g,'identity_eligible']) or bool(a.var.loc[g,'null_reference_eligible']) for g in SEX_GENES if g in a.var_names):
        raise ValueError('sex controls cannot enter identity/null-reference fits')
    return a
