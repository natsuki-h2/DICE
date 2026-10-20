"""Plot the controlled target + bright-neighbour sidelobe pilot."""

from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/ghost_matrix/controlled_sidelobe_pair_test"
OUT = INPUT / "figures"
CONFIGS = ("F1", "F5", "F6")
BRIGHT = (10.0, 30.0, 100.0)
SEPARATIONS = (2.0, 5.0, 10.0, 20.0)


def load(config):
    with (INPUT/config.lower()/"pair_results.csv").open(newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    OUT.mkdir(parents=True, exist_ok=True); metric_rows=[]
    fig, axes = plt.subplots(2, 3, figsize=(14.2, 8.2), constrained_layout=True)
    for ci, config in enumerate(CONFIGS):
        rows=load(config); dc=np.zeros((3,4)); dp=np.zeros_like(dc)
        for bi, bright in enumerate(BRIGHT):
            for si, separation in enumerate(SEPARATIONS):
                cell=[r for r in rows if float(r["target_snr"])==4 and
                      float(r["bright_snr"])==bright and float(r["separation_beams"])==separation]
                isolated=np.mean([r["isolated_recovered"].lower()=="true" for r in cell])
                combined=np.mean([r["combined_recovered"].lower()=="true" for r in cell])
                dc[bi,si]=combined-isolated
                dp[bi,si]=np.median([float(r["delta_peak_snr"]) for r in cell])
                metric_rows.append({"configuration":config,"target_snr":4,
                    "bright_snr":bright,"separation_beams":separation,
                    "isolated_completeness":isolated,"combined_completeness":combined,
                    "delta_completeness":dc[bi,si],"median_delta_peak_snr":dp[bi,si]})
        for row_index,(array,title,vmin,vmax,cmap,fmt) in enumerate([
            (dc,"Change in completeness",-.4,.4,"RdBu_r","{:+.2f}"),
            (dp,"Median change in target peak S/N",-3,3,"RdBu_r","{:+.2f}")]):
            ax=axes[row_index,ci]
            im=ax.imshow(array,origin="lower",aspect="auto",cmap=cmap,vmin=vmin,vmax=vmax)
            ax.set(xticks=range(4),xticklabels=[str(int(v)) for v in SEPARATIONS],
                   yticks=range(3),yticklabels=[str(int(v)) for v in BRIGHT],
                   xlabel="Separation (beams)",ylabel="Bright-source S/N",
                   title=f"{config}: {title}")
            for y in range(3):
                for x in range(4):
                    ax.text(x,y,fmt.format(array[y,x]),ha="center",va="center",fontsize=8,
                            color="white" if abs(array[y,x])>.55*max(abs(vmin),abs(vmax)) else "black")
            fig.colorbar(im,ax=ax,shrink=.82,label="ΔC" if row_index==0 else "Δ peak S/N")
    fig.suptitle("Controlled bright-neighbour pilot: target S/N = 4, paired noise realizations")
    fig.savefig(OUT/"controlled_sidelobe_heatmaps.png",dpi=200)
    fig.savefig(OUT/"controlled_sidelobe_heatmaps.pdf"); plt.close(fig)

    fig,axes=plt.subplots(1,3,figsize=(13.8,4.3),constrained_layout=True)
    for ax,config in zip(axes,CONFIGS):
        payload=np.load(INPUT/config.lower()/"psf_display.npz"); psf=payload["psf"]
        cell=float(payload["cell_arcsec"]); beam=np.sqrt(float(payload["beam_major_arcsec"])*float(payload["beam_minor_arcsec"]))
        center=np.asarray(psf.shape)//2; radius=int(np.ceil(22*beam/cell))
        crop=psf[center[0]-radius:center[0]+radius+1,center[1]-radius:center[1]+radius+1]
        extent=np.array([-radius,radius,-radius,radius])*cell/beam
        im=ax.imshow(crop.T,origin="lower",extent=extent,cmap="RdBu_r",vmin=-.25,vmax=.25)
        for separation in SEPARATIONS:
            ax.add_patch(plt.Circle((0,0),separation,fill=False,color="black",lw=.5,alpha=.5))
        ax.set(title=f"{config} dirty PSF",xlabel="Offset (beams)",ylabel="Offset (beams)")
        fig.colorbar(im,ax=ax,shrink=.8,label="Normalized response (core clipped)")
    fig.suptitle("Native dirty PSFs and tested pair separations")
    fig.savefig(OUT/"controlled_sidelobe_psfs.png",dpi=200); plt.close(fig)

    with (OUT/"controlled_sidelobe_metrics.csv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(metric_rows[0])); writer.writeheader(); writer.writerows(metric_rows)
    print(OUT/"controlled_sidelobe_heatmaps.png")
    print(OUT/"controlled_sidelobe_psfs.png")


if __name__=="__main__": main()
