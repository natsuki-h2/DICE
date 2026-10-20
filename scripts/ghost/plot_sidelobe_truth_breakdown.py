"""Plot injection-matched recovery and contamination for bright-neighbour tests."""

from pathlib import Path
import argparse
import csv

import matplotlib.pyplot as plt
import numpy as np
from bliss_alma.provenance import write_script_provenance

ROOT=Path(__file__).resolve().parents[2]
INPUT=ROOT/"data/ghost_matrix/controlled_sidelobe_pair_test"
OUT=INPUT/"figures"
CONFIGS=("F1","F5","F6"); BRIGHT=(10.,30.,100.); SEPARATIONS=(2.,5.,10.,20.)

def load(config):
    with (INPUT/config.lower()/"pair_results.csv").open(newline="") as handle:
        return list(csv.DictReader(handle))

def main():
    global INPUT,OUT
    parser=argparse.ArgumentParser()
    parser.add_argument("--input-root",type=Path,default=INPUT)
    args=parser.parse_args()
    INPUT=args.input_root.resolve(); OUT=INPUT/"figures"
    OUT.mkdir(parents=True,exist_ok=True); summary=[]
    fig,axes=plt.subplots(3,3,figsize=(14.2,11.0),constrained_layout=True)
    for ci,config in enumerate(CONFIGS):
        rows=load(config); false=np.zeros((3,4)); boosted=np.zeros((3,4)); unmatched=np.zeros((3,4))
        for bi,bright in enumerate(BRIGHT):
            for si,separation in enumerate(SEPARATIONS):
                cell=[r for r in rows if float(r["target_snr"])==4 and
                      float(r["bright_snr"])==bright and float(r["separation_beams"])==separation]
                false[bi,si]=np.mean([r["recovery_class"]=="sidelobe_only_false_association" for r in cell])
                boosted[bi,si]=np.mean([r["recovery_class"]=="sidelobe_boosted_recovery" for r in cell])
                unmatched[bi,si]=np.median([int(r["combined_n_unmatched"]) for r in cell])
                classes={name:np.mean([r["recovery_class"]==name for r in cell]) for name in
                    ("robust_target_recovery","sidelobe_only_false_association",
                     "sidelobe_boosted_recovery","sidelobe_suppressed_recovery","missed")}
                summary.append({"configuration":config,"target_snr":4,"bright_snr":bright,
                    "separation_beams":separation,**classes,
                    "median_unmatched_detections":unmatched[bi,si],
                    "median_truth_contamination":np.median([float(r["combined_truth_contamination"]) for r in cell])})
        for ri,(array,title,vmin,vmax,cmap,fmt,label) in enumerate([
            (false,"Bright-only false association fraction",0,1,"magma","{:.2f}","Fraction"),
            (boosted,"Sidelobe-boosted recovery fraction",0,1,"viridis","{:.2f}","Fraction"),
            (np.log10(unmatched+1),"Unmatched detections per 512-pixel cutout",0,3,"cividis","{:.0f}","log10(N+1)")]):
            ax=axes[ri,ci]; im=ax.imshow(array,origin="lower",aspect="auto",cmap=cmap,vmin=vmin,vmax=vmax)
            ax.set(xticks=range(4),xticklabels=[str(int(v)) for v in SEPARATIONS],
                   yticks=range(3),yticklabels=[str(int(v)) for v in BRIGHT],
                   xlabel="Separation (beams)",ylabel="Bright-source S/N",title=f"{config}: {title}")
            annotation=unmatched if ri==2 else array
            for y in range(3):
                for x in range(4): ax.text(x,y,fmt.format(annotation[y,x]),ha="center",va="center",fontsize=8,
                                             color="white" if array[y,x]>.58*vmax else "black")
            fig.colorbar(im,ax=ax,shrink=.82,label=label)
    n_realizations=len({int(r["realization"]) for r in load(CONFIGS[0])})
    fig.suptitle(f"Injection-matched sidelobe outcomes: target S/N = 4; {n_realizations} realizations")
    fig.savefig(OUT/"controlled_sidelobe_truth_breakdown.png",dpi=200)
    fig.savefig(OUT/"controlled_sidelobe_truth_breakdown.pdf"); plt.close(fig)
    with (OUT/"controlled_sidelobe_truth_breakdown.csv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(summary[0])); writer.writeheader(); writer.writerows(summary)
    write_script_provenance(OUT/"controlled_sidelobe_truth_breakdown.provenance.json",root=ROOT,script=__file__,
        inputs=[INPUT],outputs=[OUT/"controlled_sidelobe_truth_breakdown.png",OUT/"controlled_sidelobe_truth_breakdown.pdf",OUT/"controlled_sidelobe_truth_breakdown.csv"],
        parameters={"target_snr":4,"classification":"injection-matched"})
    print(OUT/"controlled_sidelobe_truth_breakdown.png")

if __name__=="__main__": main()
