"""Validate one image-domain pair model against direct visibility injection."""

from pathlib import Path
import argparse
import json
import math
import shutil

import numpy as np
from casatasks import ft, imstat, tclean, uvsub
from casatools import componentlist, image

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"data/ghost_matrix/controlled_sidelobe_pair_test"
BEAM_MINOR={"F5":0.04867600277066231,"F6":0.01248726062476635}
BEAM_MAJOR={"F5":0.06545752286911011,"F6":0.013991966843605042}
BASE={"F5":ROOT/"data/ghost_matrix/f5_5min/extracted/output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms",
      "F6":ROOT/"data/ghost_matrix/f6_5min/extracted/output_custom_skymodel_SIM_MIXED_5m_-27_480000_ALMA.ms"}
IMSIZE=512

def remove(path):
    if path.is_dir(): shutil.rmtree(path)
    elif path.exists(): path.unlink()

def remove_images(prefix):
    for path in prefix.parent.glob(prefix.name+".*"): remove(path)

def read(path):
    ia=image()
    try: ia.open(str(path)); return np.squeeze(np.asarray(ia.getchunk(),dtype=np.float32))
    finally: ia.done()

def shifted(array,shift):
    result=np.zeros_like(array); source=[]; destination=[]
    for size,delta in zip(array.shape,shift):
        if delta>=0: source.append(slice(0,size-delta)); destination.append(slice(delta,size))
        else: source.append(slice(-delta,size)); destination.append(slice(0,size+delta))
    result[tuple(destination)]=array[tuple(source)]; return result

def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--config",choices=("F5","F6"),required=True)
    args=parser.parse_args(); config=args.config; cell=BEAM_MINOR[config]/4
    out=OUT/config.lower(); work=out/"validation_work"; remove(work); work.mkdir()
    ms=work/"pair.ms"; shutil.copytree(BASE[config],ms)
    blank_prefix=work/"blank"; pair_prefix=work/"pair"; components=work/"pair.cl"
    common=dict(vis=str(ms),specmode="mfs",imsize=[IMSIZE,IMSIZE],cell=f"{cell}arcsec",
                weighting="natural",robust=.5,niter=0,pblimit=.2,gridder="standard",interactive=False)
    tclean(imagename=str(blank_prefix),**common)
    blank=read(Path(str(blank_prefix)+".image")); psf=read(Path(str(blank_prefix)+".psf"))
    rms=float(imstat(str(blank_prefix)+".image",listit=False)["rms"][0]); center=(IMSIZE//2,IMSIZE//2)
    peak=np.unravel_index(np.argmax(psf),psf.shape); psf=np.roll(psf,(center[0]-peak[0],center[1]-peak[1]),axis=(0,1))
    beam=math.sqrt(BEAM_MAJOR[config]*BEAM_MINOR[config]); dx=int(round(10*beam/cell))
    east_arcsec=-dx*cell; ra=180+east_arcsec/(3600*math.cos(math.radians(-27.8)))
    cl=componentlist()
    try:
        cl.addcomponent(dir="J2000 180deg -27.8deg",flux=4*rms,fluxunit="Jy",shape="point")
        cl.addcomponent(dir=f"J2000 {ra}deg -27.8deg",flux=100*rms,fluxunit="Jy",shape="point")
        cl.rename(str(components))
    finally: cl.done()
    ft(vis=str(ms),spw="0:0",complist=str(components),usescratch=True); uvsub(vis=str(ms),reverse=True)
    tclean(imagename=str(pair_prefix),**common); actual=read(Path(str(pair_prefix)+".image"))
    predicted_plus=blank+4*rms*psf+100*rms*shifted(psf,(dx,0))
    predicted_minus=blank+4*rms*psf+100*rms*shifted(psf,(-dx,0))
    valid=np.isfinite(actual)&np.isfinite(blank); margin=abs(dx)+16
    region=np.zeros_like(valid); region[margin:-margin,margin:-margin]=True; valid&=region
    def metrics(predicted):
        diff=(actual-predicted)[valid]/rms
        return {"difference_rms_sigma":float(np.sqrt(np.mean(diff**2))),
                "maximum_absolute_difference_sigma":float(np.max(np.abs(diff)))}
    result={"configuration":config,"target_snr":4,"bright_snr":100,"separation_beams":10,
            "dx_pixels":dx,"same_shift":metrics(predicted_plus),"opposite_shift":metrics(predicted_minus)}
    (out/"visibility_validation.json").write_text(json.dumps(result,indent=2)+"\n")
    remove(work); print(json.dumps(result,indent=2))

if __name__=="__main__": main()
