# DICE terminology

Use direct, measurable terms in documentation, command output, and figures.
The following conventions keep the software documentation consistent with the
manuscript and reduce ambiguity for new users.

| Use | Avoid or restrict | Reason |
|---|---|---|
| injection--recovery framework | source finder, when referring to DICE as a whole | Source extraction is one measurement stage. DICE uses completeness, contamination, flux recovery, and astrometric accuracy to connect a measured catalogue to corrected source counts and population constraints. |
| untargeted search | blind search | State that no individual source position or line frequency was preselected. Use the historical term only when quoting a paper or legacy filename. |
| injected-source catalogue or injected position | truth catalogue or truth position | Name the information that is actually available in a simulation. |
| flux recovery ratio | flux fidelity | Define the quantity explicitly as recovered flux divided by injected flux. The function name `flux_fidelity` remains only for backward compatibility. |
| astrometric accuracy or positional offset | position fidelity | Report the angular offset or the offset normalized by the synthesized beam. |
| catalogue fidelity | — | Retain this established literature term only for a catalogue reliability estimator such as `1 - N_negative / N_positive`. Do not use it for flux or position measurements. |
| completeness--contamination relation | trajectory | State which two measured quantities are plotted. |
| forward modelling of source populations | population-forward modelling | Put the scientific object before the operation. Existing module and directory names are retained for compatibility. |
| diagnostic comparison with specified fixed inputs | controlled stress test | State what is varied, what is fixed, and whether the cases are empirical or synthetic. |
| uncertainty estimated by resampling complete realizations | whole-realization uncertainty | Describe the resampling unit instead of introducing an unexplained label. |

Terms such as *robust RMS*, CASA `robust` weighting, and seed/flood thresholds
have specific technical meanings and should not be simplified. Likewise,
legacy filenames and configuration keys should be changed only with a
documented compatibility path.
