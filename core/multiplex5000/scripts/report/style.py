"""Drawing rules for the report, enforced by reuse rather than restated per figure.

* **Status colours** are pinned to the app's display contract
  (``apps/islet/config/islet_default.yml`` -> ``grouping.default_palette``) and the order to
  ``phenocycler.config.STATUS_ORDER``.  ``AAB`` is displayed as ``Aab+`` in exactly one place,
  :func:`status_display`; ``tests/test_style.py`` compares the hex strings to the YAML text.
* **Image panels** go through :func:`image_panel`, which accepts only a
  ``scripts/restain/display.DisplayRange`` -- ``vmin`` is 0 by construction, the white point and
  ``frac_at_zero`` are printed in the panel caption, ``interpolation='nearest'``, gamma 1.  Nothing
  in this lane calls ``imshow`` directly; the tests grep for it.
* **Axes start at 0**, the measured floor.  Where the range demands it the axis is ``asinh``
  (``scripts/celltype/figures._asinh_axis``), and the axis label says which.
* **Unavailable is grey, with the reason in the legend**; Ambiguous is hatched.  Neither is a
  bucket: every grey patch names why.
* **Every figure** carries ``pins.DEV_LINE`` as its second title line and a provenance footer.
"""
from __future__ import annotations

import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

# `import pins` used to be relied on for the side effect of putting scripts/restain and
# scripts/celltype on sys.path (pins.py does that itself). That only happens the FIRST time
# 'pins' is imported in the process; when this lane is collected alongside another that has
# already cached its own same-named 'pins' module, the side effect never runs and `import
# display` fails with "No module named 'display'" even though scripts/restain/display.py
# exists. Insert both directories explicitly so this module does not depend on import order.
_SCRIPTS = pathlib.Path(__file__).resolve().parent.parent
for _sub in ("restain", "celltype"):
    _p = str(_SCRIPTS / _sub)
    if _p not in sys.path:
        sys.path.insert(0, _p)

import pins  # noqa: E402
import display as DISP  # noqa: E402  (scripts/restain/display.py)
from figures import _asinh_axis  # noqa: E402  (scripts/celltype/figures.py)
from phenocycler.config import STATUS_ORDER  # noqa: E402
from phenocycler.integration.figures import COMMON_COLORS  # noqa: E402

DisplayRange = DISP.DisplayRange
unsigned_display_range = DISP.unsigned_display_range
UNAVAILABLE = DISP.UNAVAILABLE                     # "#BBBBBB"
UNAVAILABLE_ALT = "#7F7F7F"                        # a second Unavailable reason on the same panel
AMBIGUOUS_HATCH = "///"

#: ``islet_default.yml`` -> grouping.palettes[grouping.default_palette]; keyed by the STORED token.
STATUS_PALETTE_NAME = "Bright"
STATUS_COLORS = {"ND": "#3366CC", "AAB": "#FF6600", "T1D": "#109618"}
_STATUS_DISPLAY = {"ND": "ND", "AAB": "Aab+", "T1D": "T1D"}

#: Eligible / excluded pair (pheno_qc.py).
ELIGIBLE = "#4c78a8"
EXCLUDED = "#e45756"
ROLE_COLORS = {"analysis": "#4c78a8", "estimation": "#f2cf5b", "spillover_context": "#59a14f"}
THRESHOLD = "#222222"
MODE = "#d62728"

DPI = 200


def status_display(status: str) -> str:
    """The ONE place the stored ``AAB`` token becomes the display string ``Aab+``."""
    s = str(status)
    if s in _STATUS_DISPLAY:
        return _STATUS_DISPLAY[s]
    if s in _STATUS_DISPLAY.values():          # already a display string
        return s
    raise KeyError(f"unknown donor status {status!r}; STATUS_ORDER is {STATUS_ORDER}")


def status_color(status: str) -> str:
    s = str(status)
    if s not in STATUS_COLORS:
        inv = {v: k for k, v in _STATUS_DISPLAY.items()}
        s = inv[s]
    return STATUS_COLORS[s]


def status_axis() -> list:
    """The ND / Aab+ / T1D axis in STATUS_ORDER, as display strings."""
    return [status_display(s) for s in STATUS_ORDER]


# --------------------------------------------------------------------------- figure frame
def figure(nrows: int, ncols: int, *, width: float = 5.2, height: float = 4.2, **kw):
    fig, axes = plt.subplots(nrows, ncols, figsize=(width * ncols, height * nrows), squeeze=False, **kw)
    return fig, axes


def title(fig, line1: str) -> None:
    fig.suptitle(f"{line1}\n{pins.DEV_LINE}", fontsize=11, fontweight="bold", y=0.995)


def footer_text(prov: dict, donor: str | None = None) -> str:
    ids = prov.get("stage_content_ids", {})
    bits = [f"run {prov['run_id']}", f"run.json {prov['run_content_id'][:12]}",
            f"expression manifest {prov['expression_content_id'][:12]} / output {prov['expression_output_sha256'][:12]}"]
    for s in ("ingest", "geometry", "redsea"):
        if ids.get(s):
            bits.append(f"{s} {ids[s][:12]}")
    if prov.get("uniform_commit"):
        bits.append(f"UniFORM {prov['uniform_commit'][:8]}")
    line1 = " · ".join(bits)
    d_bits = []
    donors = [donor] if donor else list(prov.get("donors", {}))
    for d in donors:
        row = prov.get("donors", {}).get(str(d), {})
        if row.get("base_sha256"):
            d_bits.append(f"Xenium base {row.get('serial', '')} sha {row['base_sha256'][:12]}")
    code = prov.get("code", {})
    d_bits.append(f"code {code.get('commit', '')[:10]}{' (dirty)' if code.get('dirty') else ''}")
    rel = prov.get("identity_releases") or {}
    rel_txt = " / ".join(f"{m} {r}" for m, r in rel.items() if r) or "none (QC figures need none)"
    d_bits.append(f"identity release: {rel_txt}")
    d_bits.append(prov.get("written_utc", ""))
    return line1 + "\n" + " · ".join(d_bits)


def footer(fig, prov: dict, donor: str | None = None) -> None:
    fig.text(0.005, 0.003, footer_text(prov, donor), fontsize=6, family="monospace", va="bottom", ha="left",
             color="#333333")


def finish(fig, *, top: float = 0.90, bottom: float = 0.06) -> None:
    fig.tight_layout(rect=(0, bottom, 1, top))


def savefig(fig, out_dir: pathlib.Path, stem: str) -> list:
    """PNG at 200 dpi plus a PDF, same stem.  Returns the paths written."""
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    png, pdf = out_dir / f"{stem}.png", out_dir / f"{stem}.pdf"
    fig.savefig(png, dpi=DPI)
    fig.savefig(pdf)
    plt.close(fig)
    return [png, pdf]


def denominator(ax, text: str, *, loc: str = "upper right") -> None:
    """Every proportion or distribution panel prints what it is a proportion OF."""
    x, ha = (0.99, "right") if loc.endswith("right") else (0.01, "left")
    y, va = (0.98, "top") if loc.startswith("upper") else (0.02, "bottom")
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=7,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#999999", lw=0.5, alpha=0.9))


def wrap(text: str, width: int = 88) -> str:
    """Hard-wrap a panel title or note so it never runs past its axes (matplotlib's wrap=True is unreliable)."""
    import textwrap
    return "\n".join(textwrap.fill(line, width) for line in str(text).split("\n"))


def note(ax, text: str, *, y: float = 0.5, fontsize: int = 8, width: int = 60) -> None:
    ax.text(0.5, y, wrap(text, width), transform=ax.transAxes, ha="center", va="center", fontsize=fontsize)


def unavailable_patch(reason: str, color: str = UNAVAILABLE) -> Patch:
    return Patch(facecolor=color, edgecolor="#555555", lw=0.5, label=f"Unavailable: {reason}")


def ambiguous_patch(reason: str) -> Patch:
    return Patch(facecolor="white", edgecolor="#555555", hatch=AMBIGUOUS_HATCH, label=f"Ambiguous: {reason}")


# --------------------------------------------------------------------------- axes from zero
def range_demands_asinh(values) -> bool:
    """A stated rule, not a per-figure choice: asinh when the max is far beyond the bulk."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 10:
        return False
    p99 = float(np.percentile(v, 99))
    return bool(p99 > 0 and float(v.max()) > 25.0 * p99)


def asinh_axis(ax, hi: float, axis: str = "x") -> None:
    """From 0 (the measured floor) to ``hi`` on an arcsinh scale; reuses the celltype helper for x."""
    if axis == "x":
        _asinh_axis(ax, hi)
    else:
        ax.set_yscale("asinh", linear_width=max(hi / 1e3, 1.0))
        ax.set_ylim(0, hi)


def hist_from_zero(ax, values, *, hi: float | None = None, bins: int = 120, asinh: bool | None = None,
                   color: str = ELIGIBLE, label: str | None = None, alpha: float = 0.85, density: bool = False):
    """Histogram whose axis starts at 0.  Values above ``hi`` are counted and folded into the last bin
    (the count is printed by the caller via the returned dict); NaN is never drawn as 0."""
    v = np.asarray(values, dtype=float)
    n_nan = int((~np.isfinite(v)).sum())
    v = v[np.isfinite(v)]
    if hi is None:
        hi = float(np.percentile(v, 99.9)) * 1.2 if v.size else 1.0
    hi = max(float(hi), 1e-9)
    if asinh is None:
        asinh = range_demands_asinh(v)
    n_over = int((v > hi).sum())
    vv = np.minimum(v, hi)
    if asinh:
        edges = np.sinh(np.linspace(0.0, np.arcsinh(hi), bins + 1))
    else:
        edges = np.linspace(0.0, hi, bins + 1)
    ax.hist(vv, bins=edges, color=color, alpha=alpha, label=label, density=density)
    if asinh:
        asinh_axis(ax, hi, "x")
    else:
        ax.set_xlim(0, hi)
    return {"n": int(v.size), "n_nan": n_nan, "n_over_axis": n_over, "hi": hi, "asinh": bool(asinh)}


def axis_label_from_zero(name: str, unit: str, asinh: bool) -> str:
    scale = "asinh" if asinh else "linear"
    unit_part = f"{unit}; " if unit else ""
    return f"{name} ({unit_part}{scale} axis from 0)"


def threshold_line(ax, x: float, text: str, *, color: str = THRESHOLD, ls: str = "--") -> None:
    ax.axvline(x, color=color, ls=ls, lw=1.0)
    ax.text(x, 0.97, f" {text}", transform=ax.get_xaxis_transform(), fontsize=7, va="top", ha="left",
            color=color, rotation=90)


# --------------------------------------------------------------------------- image panels
def image_panel(ax, img, rng, title_text: str, *, extent=None, cmap: str = "gray", caption_extra: str = ""):
    """The only path to ``imshow`` in this lane.  Refuses anything but a ``DisplayRange`` (vmin 0)."""
    if not isinstance(rng, DisplayRange):
        raise DISP.DisplayError("image panels take a display.DisplayRange; no ad-hoc vmin/vmax")
    if rng.vmin != 0.0:
        raise DISP.DisplayError("black point must be the measured floor 0")
    im = ax.imshow(np.asarray(img), cmap=cmap, vmin=rng.vmin, vmax=rng.vmax, interpolation="nearest",
                   extent=extent)
    ax.set_xticks([]); ax.set_yticks([])
    cap = rng.caption() + (f"; {caption_extra}" if caption_extra else "")
    if rng.fallback_reason:
        cap += f"; {rng.fallback_reason}"
    ax.set_title(f"{title_text}\n{cap}", fontsize=8)
    return im


def filled_panel(ax, label_img, lut, rng, title_text: str, *, unavailable_mask=None, cmap: str = "magma",
                 extent=None, caption_extra: str = ""):
    """Cells filled by a per-cell value: ``lut[label]`` per pixel, NaN drawn as Unavailable grey.

    The white point is a ``DisplayRange`` over the SELECTION-INDEPENDENT per-cell population (the
    caller decides which; the caption states it), so two panels sharing a range are comparable and
    zero stays at 0.  ``unavailable_mask`` marks pixels of a second Unavailable reason (drawn
    ``UNAVAILABLE_ALT``); the legend is the caller's, via :func:`unavailable_patch`.
    """
    if not isinstance(rng, DisplayRange):
        raise DISP.DisplayError("filled panels take a display.DisplayRange; no ad-hoc vmin/vmax")
    lab = np.asarray(label_img)
    vals = np.asarray(lut, dtype=float)[lab]
    vals = np.where(lab > 0, vals, np.nan)
    cm = DISP.mapped_colormap(cmap)                    # bad = UNAVAILABLE grey
    im = ax.imshow(vals, cmap=cm, vmin=rng.vmin, vmax=rng.vmax, interpolation="nearest", extent=extent)
    if unavailable_mask is not None:
        alt = np.zeros((*lab.shape, 4), dtype=float)
        rgba = matplotlib.colors.to_rgba(UNAVAILABLE_ALT)
        alt[np.asarray(unavailable_mask, dtype=bool)] = rgba
        ax.imshow(alt, interpolation="nearest", extent=extent)
    bg = np.zeros((*lab.shape, 4), dtype=float)
    bg[lab == 0] = matplotlib.colors.to_rgba("#f8f9fa")
    ax.imshow(bg, interpolation="nearest", extent=extent)
    ax.set_xticks([]); ax.set_yticks([])
    cap = rng.caption() + (f"; {caption_extra}" if caption_extra else "")
    ax.set_title(f"{title_text}\n{cap}", fontsize=8)
    return im


def scale_bar(ax, px_per_um: float, length_um: float = 200.0, *, x_frac: float = 0.04, y_frac: float = 0.95,
              color: str = "white") -> None:
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    w, h = x1 - x0, y1 - y0
    bx = x0 + x_frac * w
    by = y0 + y_frac * h
    ax.plot([bx, bx + length_um * px_per_um], [by, by], color=color, lw=3, solid_capstyle="butt")
    ax.text(bx, by - 0.02 * h, f"{length_um:.0f} µm", color=color, fontsize=7, va="bottom", ha="left")


def core_outline(ax, pad: int, size: int, color: str = "#f76707") -> None:
    """Outline the scored core inside a padded window (frame idiom: 1024 core in a 2048 window)."""
    from matplotlib.patches import Rectangle
    ax.add_patch(Rectangle((pad - 0.5, pad - 0.5), size, size, fill=False, ec=color, lw=0.8, ls="--"))
