"""Create a publication-ready flowchart of the completed GHOST F1 workflow."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from bliss_alma.provenance import write_script_provenance


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data" / "ghost_f1_5min" / "f1_snr_pb_grid" / "figures"


def box(ax, x, y, width, height, text, facecolor, edgecolor="#333333"):
    patch = FancyBboxPatch(
        (x - width / 2, y - height / 2), width, height,
        boxstyle="round,pad=0.02,rounding_size=0.025",
        linewidth=1.25, edgecolor=edgecolor, facecolor=facecolor,
    )
    ax.add_patch(patch)
    ax.text(x, y, text, ha="center", va="center", fontsize=9)


def arrow(ax, start, end):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=12,
                                linewidth=1.2, color="#555555",
                                connectionstyle="arc3,rad=0"))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(13.2, 7.2), constrained_layout=True)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    blue, green, amber, rose, grey = "#DCEAF7", "#DDEFD9", "#F7E9BF", "#F5DCDD", "#EEEEEE"
    w, h = 0.15, 0.105

    box(ax, .09, .83, w, h, "Zero-valued\n2D FITS", grey)
    box(ax, .28, .83, w, h, "GHOST F1\n5 min, 230 GHz", blue)
    box(ax, .47, .83, w, h, "Blank MS +\nblank dirty image", blue)
    arrow(ax, (.165, .83), (.205, .83))
    arrow(ax, (.355, .83), (.395, .83))

    box(ax, .66, .72, w, h, "Inject 44 sources\n11 S/N × 4 PB bins", green)
    box(ax, .84, .72, w, h, "CASA dirty imaging\n10 batches", green)
    arrow(ax, (.53, .80), (.595, .745))
    arrow(ax, (.735, .72), (.765, .72))

    box(ax, .84, .49, w, h, "Positive-clump extraction\n4σ seed / 3σ flood", green)
    box(ax, .66, .49, w, h, "Match detections\nto input positions", green)
    arrow(ax, (.84, .665), (.84, .545))
    arrow(ax, (.765, .49), (.735, .49))

    box(ax, .48, .62, w, h, "Completeness\nby input S/N and PB", amber)
    box(ax, .48, .40, w, h, "Input S/N → output S/N\nselection-bias check", amber)
    arrow(ax, (.595, .515), (.555, .59))
    arrow(ax, (.595, .465), (.555, .425))

    box(ax, .28, .28, w, h, "Invert blank image", rose)
    box(ax, .47, .19, w, h, "Negative-clump catalogue\ndiagnostic floor: 2σ", rose)
    box(ax, .68, .19, w, h, "Nmock / Npositive\ndifferential + cumulative Q", amber)
    arrow(ax, (.44, .79), (.31, .335))
    arrow(ax, (.355, .255), (.405, .215))
    arrow(ax, (.545, .19), (.605, .19))
    arrow(ax, (.82, .435), (.715, .245))

    ax.text(.84, .60, "PRODUCTION DETECTION", ha="center", va="center",
            fontsize=9, color="#2A6F3A")
    ax.text(.39, .07, "CONTAMINATION DIAGNOSTIC — does not change the adopted 4σ threshold",
            ha="center", va="center", fontsize=9, color="#8B3A3A")
    ax.set_title("DICE GHOST F1 continuum injection–recovery workflow", fontsize=16, pad=12)

    fig.savefig(OUT / "f1_workflow.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT / "f1_workflow.pdf", bbox_inches="tight")
    plt.close(fig)
    write_script_provenance(OUT/"f1_workflow.provenance.json",root=ROOT,script=__file__,
        outputs=[OUT/"f1_workflow.png",OUT/"f1_workflow.pdf"],
        parameters={"workflow":"GHOST F1 continuum injection-recovery"})
    print(OUT / "f1_workflow.png")


if __name__ == "__main__":
    main()
