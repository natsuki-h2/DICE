# Fixed-width line-injection demo

This deterministic example places five artificial point-source lines in one
cube. Their FWHM is fixed at 10 channels while peak S/N, sky position, and
centre channel vary. Other FWHM values are evaluated in separate cubes.

The left panel follows the simple position-map presentation of Hayatsu et al.
(2025, Fig. 4): colour gives input peak S/N over a PB background. The right
panel shows the corresponding Gaussian profiles and four-FWHM protected
intervals. The mixed-width multi-helix planner remains an advanced option.

From the repository root, run:

```bash
python scripts/tutorials/plot_line_injection_helix.py
```

The script writes the PNG/PDF figure, placement CSV, and provenance JSON in
this directory. The committed products provide a reference result for student
testers.
