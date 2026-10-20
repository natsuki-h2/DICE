# ALMA, ALMA2/WSU, and ALMA2040 comparison for DICE

## Terminology

`ALMA2` is the Japanese name for the enhancement phase aligned with the
international ALMA2030 Wideband Sensitivity Upgrade (WSU). It is not a new
array with a separately published, final antenna configuration. The committed
WSU program primarily replaces the signal chain and receivers while using the
present antennas and infrastructure. The earlier Japanese ALMA2 science goals
also state approximately two times sensitivity, more than two times angular
resolution, and more than two times instantaneous bandwidth; the angular-
resolution goal should not be treated as a committed WSU baseline extension.

## Facility-level comparison

| Property | Current ALMA | ALMA2 / ALMA2030 WSU | ALMA2040 assumptions used by DICE |
|---|---|---|---|
| Status | Operating facility; Cycle 12 reference | Funded/staged upgrade under development | Community technical vision and simulation assumptions, not a finalized facility |
| Antennas relevant to 12-m interferometry | Up to 50 12-m antennas; sensitivity calculations commonly assume 43. ALMA overall has 54 x 12-m plus 12 x 7-m antennas. | Present antennas and infrastructure retained; initial release after at least 36 antennas are upgraded, whose improved digital efficiency is expected to match current 43-antenna point-source sensitivity. | Six supplied DICE/GHOST configuration files each contain 255 x 12-m antennas. |
| Approx. 12-m collecting area | 5,655 m2 for 50 dishes; 4,863 m2 for a 43-antenna observing assumption | Nominal geometric area unchanged; effective sensitivity improves through receivers, digitization, bandwidth, and efficiency | 28,840 m2 geometric area; 5.1x the 50-antenna ALMA area and 5.9x the 43-antenna assumption |
| Configuration / maximum baseline | Ten representative 12-m configurations; physical design reaches about 16.2 km, although Cycle 12 offers up to 8.5 km | No separately defined WSU configuration set in the current international baseline; WSU is principally a signal-chain upgrade | F1-F6 maximum baselines approximately 0.39, 1.09, 3.49, 12.51, 20.78, and 28.37 km |
| Band-6 angular resolution reference | About 0.018 arcsec for the full-design longest configuration; Cycle-dependent offered resolution is coarser when C-9/C-10 are unavailable | Same array geometry under the currently defined WSU; older Japanese ALMA2 goal says more than 2x angular resolution, but implementation/configuration is not fixed in current WSU documentation | At 230 GHz: compact 0.68 arcsec, intermediate 0.097 arcsec, extended 0.013 arcsec in the White Paper guidelines; measured native F1-F6 beams span roughly 0.72 to 0.014 arcsec |
| Continuum bandwidth convention | 7.5 GHz usable continuum bandwidth | Initial 2x system bandwidth; final goal 4x. In current-continuum convention this is approximately 15 GHz initially and 30 GHz at goal, receiver dependent. | At least 128 GHz instantaneous bandwidth anticipated; goal 512 GHz through simultaneous multi-band observing. GHOST simulations currently still use 7.5 GHz. |
| Spectral capability | Bandwidth and spectral resolution trade against one another in the current correlator | ATAC/TPGS intended to process the full instantaneous bandwidth at high spectral resolution; at least 0.1 km/s across full bandwidth and access to at least 15 m/s modes | Guideline target is 0.1 km/s over 35-1000 GHz; implementation and correlator architecture are not fixed |
| Published performance shorthand | Cycle-12 Band-6 reference: about 0.12 mJy/beam in 60 s with 43 antennas and 7.5 GHz, or about 15.5 microJy/beam in one hour by square-root time scaling | Initial WSU: at least 3x continuum imaging speed and at least 2x spectral-line imaging speed; Japanese ALMA2 goal: about 2x sensitivity | Guideline Table 3: 1.56 microJy/beam in one hour at Band 6 and 7.5 GHz for the proposed configurations |
| Rollout | Current | Initial 2x-bandwidth science observations targeted for the early 2030s; staged completion during the 2030s | 2040s vision; specifications remain subject to the White Paper process |

## Directly comparable scaling at 230 GHz

The following are planning calculations, not guaranteed performance:

- Increasing from the 43-antenna sensitivity assumption to 255 identical
  antennas improves point-source RMS by approximately
  `sqrt[255*254/(43*42)] = 5.99` before receiver or bandwidth improvements.
- Increasing usable bandwidth from 7.5 to 128 GHz gives another ideal
  radiometer improvement of `sqrt(128/7.5) = 4.13`.
- Applying only that bandwidth scaling to the ALMA2040 guideline value of
  1.56 microJy/beam gives approximately 0.38 microJy/beam in one hour. This is
  an inference and assumes unchanged system temperature, efficiency, usable
  bandwidth, weighting, and calibration performance.
- The 512-GHz goal is simultaneous multi-band coverage and should not be
  inserted as one contiguous same-noise continuum band without a defined
  receiver and combination model.

## Recommended DICE comparison

For a defensible recovery comparison, hold sky frequency, source population,
declination, integration time, weighting, detection algorithm, and PB cut
fixed, while using three separate measurement-set families:

1. current ALMA configuration files and current 7.5-GHz sensitivity;
2. the same ALMA geometry with WSU/ALMA2 bandwidth and receiver assumptions;
3. the supplied 255-antenna ALMA2040 F1-F6 configurations, first at 7.5 GHz
   for geometry-only comparison and then at 128 GHz for the facility forecast.

This separates collecting area/configuration, bandwidth, and receiver
sensitivity instead of changing all three simultaneously.

## Sources

- [ALMA Cycle 12 Technical Handbook](https://almascience.eso.org/documents-and-tools/cycle12/alma-technical-handbook/view)
- [ALMA Cycle 12 Science Primer](https://almascience.nrao.edu/documents-and-tools/cycle12/alma-science-primer)
- [ALMA Observatory factsheet](https://www.almaobservatory.org/en/factsheet/)
- [ALMA2030 WSU overview](https://www.almaobservatory.org/en/scientists/alma-2030-wsu/)
- [ALMA2030 WSU program and rollout](https://www.almaobservatory.org/en/scientists/alma-2030-wsu/wsu-program/)
- [ALMA2030 WSU technology](https://www.almaobservatory.org/en/scientists/alma-2030-wsu/wsu-technology/)
- [NAOJ description of ALMA2](https://alma-telescope.jp/en/news/alma2-project-202202.html)
- Local ALMA2040 guideline: `docs/ALMA2040____White_paper_Document.pdf`
- Local proposed configuration files: `configs/alma2040/alma.F1.cfg` through
  `alma.F6.cfg`
