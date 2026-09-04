"""Replace the area-filter diagnostic cell in H18 so the x-axis is in mm²
(log scale), making the cutoff line and the legend "5.0 mm²" align 1:1.

Idempotent: identifies the target cell by the marker substring
`HE_area_filter_diagnostic` and replaces its source with the mm²-axis version.
"""
from pathlib import Path
import nbformat as nbf

PROJECT = Path(__file__).resolve().parent.parent
NB_PATH = PROJECT / "analysis" / "H18_HE_follicles.ipynb"
MARKER = "HE_area_filter_diagnostic"

NEW_SOURCE = """follicles_all = raw[raw['Classification'] == 'Follicle'].copy()
tissue_all = raw[raw['Classification'] == 'Tissue'].copy()

# Plot in mm² (log scale) so the cutoff label and the line position read 1:1.
area_mm2 = follicles_all['Area µm^2'] / 1e6
cutoff_mm2 = AREA_FILTER_UM2 / 1e6

fig, ax = plt.subplots(figsize=(8, 4))
upper = max(area_mm2.max(), cutoff_mm2) * 1.1
bins = np.logspace(np.log10(area_mm2.min()), np.log10(upper), 80)
ax.hist(area_mm2, bins=bins, color='#888', edgecolor='white')
ax.set_xscale('log')
ax.axvline(cutoff_mm2, color='crimson', ls='--', lw=2,
           label=f'Filter cutoff: {cutoff_mm2:.1f} mm²')
ax.set_xlabel('Follicle Area (mm², log scale)')
ax.set_ylabel('Count')
ax.set_title('Follicle area distribution (raw, before filter)')
ax.legend()
save_figure(fig, 'HE_area_filter_diagnostic')
plt.show()

n_drop = (follicles_all['Area µm^2'] > AREA_FILTER_UM2).sum()
drop_per_donor = follicles_all[follicles_all['Area µm^2'] > AREA_FILTER_UM2].groupby('Sample').size()
print(f'Dropping {n_drop} follicles > {AREA_FILTER_UM2/1e6:.1f} mm² (out of {len(follicles_all):,})')
if n_drop > 0:
    print('Per-donor drop counts:')
    print(drop_per_donor)
"""


def main():
    nb = nbf.read(NB_PATH, as_version=4)
    hits = 0
    for cell in nb.cells:
        if cell.cell_type != "code":
            continue
        src = cell.source if isinstance(cell.source, str) else "".join(cell.source)
        if MARKER in src:
            cell.source = NEW_SOURCE
            hits += 1
    if hits != 1:
        raise SystemExit(f"Expected exactly 1 cell matching marker '{MARKER}', found {hits}")
    nbf.write(nb, NB_PATH)
    print(f"Replaced source of 1 cell in {NB_PATH.relative_to(PROJECT)}")


if __name__ == "__main__":
    main()
