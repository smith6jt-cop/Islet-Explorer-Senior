"""Export one PDF page per pancreas donor using the spatial notebook's plot function."""

import argparse
import csv
from pathlib import Path
from tempfile import NamedTemporaryFile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from islet_trajectory_plots import load_trajectories, plot_donor_trajectory


def export_trajectories(data: Path, output: Path) -> None:
    with (data / "pancreas_manifest.csv").open(newline="") as manifest:
        donors = [row["donor"].strip() for row in csv.DictReader(manifest)]
    if not donors or not all(donors) or len(donors) != len(set(donors)):
        raise ValueError("The pancreas manifest must list each donor exactly once.")

    inputs = load_trajectories(data)
    metadata = {
        "Title": "Per-donor PhenoCycler islet trajectories",
        "Subject": "Islet trajectories and core and peri-islet immune-cell density by donor",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(dir=output.parent, suffix=".pdf", delete=False) as tmp:
        staged = Path(tmp.name)
    try:
        with matplotlib.rc_context({"pdf.fonttype": 42}), PdfPages(staged, metadata=metadata) as pages:
            for donor in donors:
                fig = plot_donor_trajectory(inputs, donor)
                try:
                    pages.savefig(fig, dpi=200)
                finally:
                    plt.close(fig)
            if pages.get_pagecount() != len(donors):
                raise ValueError("The PDF page count does not match the donor manifest.")
        staged.replace(output)
    finally:
        staged.unlink(missing_ok=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True, help="Private Scanpy bundle directory")
    parser.add_argument("--output", type=Path, required=True, help="Destination PDF")
    args = parser.parse_args()
    export_trajectories(args.data, args.output)
