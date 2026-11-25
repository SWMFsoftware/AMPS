# Input data of `make_qfit_figures.py`

Copied 2026-11-25 from `validation_runs/` (not tracked by git), so that the figure script runs from git.
Byte-identical copies, SHA-256 below. (The folder is not called `data`: the repository's
`.gitignore` excludes every directory of that name.)

| File | Content | Copied from |
|---|---|---|
| `iuvs_o1304_case5_digitized.csv` | MAVEN/IUVS O I 130.4-nm radial profile of 3I/ATLAS, 2025-10-29 17:32 (case 5, perihelion): 7 bins, bin bounds, radiance (R, kR) | `validation_runs/i3_master_campaign3_oxygen_20261119/report_oxygen/observations/` |
| `iuvs_o1304_case5_slide_model.csv` | the IUVS team model curve of the same slide panel (Q(H2O) = 1.89e29 s^-1), 153 points | same |
| `ohfit_numbers_output4.json` | OH fit of the cancelled SA1982_VH runs at output 4: fit factors, fitted Q(H2O) (the input of `validation_runs/i3_qfit_campaign_20261123`), predicted H and case-5 O model/obs | `validation_runs/i3_fixed_campaign_20261122/report_interim_20261123/data/ohfit_numbers.json` |

| File | SHA-256 |
|---|---|
| `iuvs_o1304_case5_digitized.csv` | `8cd7ec66d5465aa74f449fffddafd84cea043a28a6c76f167c3bba29d0fc0d78` |
| `iuvs_o1304_case5_slide_model.csv` | `e11d6cd982c98f7cf7bb644d6a5e45ad9d425d1b60b526f31ffb47c070d66947` |
| `ohfit_numbers_output4.json` | `76eba93e934f2d1833bc76430ce75675e6b3feebd7d0c7a1bfbd90c1152860eb` |

## The IUVS O I 130.4-nm profile

- Source: the right-hand panel of the presentation slide "MAVEN/IUVS FUV Observations at Perihelion:
  O 130.4 nm" (photographed 2026-11-18), digitized on 2026-11-21. Axes 0.03-20 x 10^5 km by 1-3000 R,
  log-log; labels "Q(H2O)=1.89e+29 s^-1", "Model", "Oct 29 17:32", "Delta=1.8e+05 km".
- Case 5 only: 2025-10-29 17:32 lies in the case-5 window of `srcHartley2/i3/reference/i3_cases.csv`
  (2025-10-28T04:59:25Z to 2025-10-31T06:05:51Z), and Delta = 1.8e5 km is case 5's resolution.
- Bins are spatial elements folded about the nucleus: bin 1 = [0, w], bin k = [(2k-3)w, (2k-1)w], with
  w = 0.918e5 km (bin width 1.836e5 km); the distance column is the bin midpoint.
- Accuracy: radiances about +-1.5 %, distances better than 1 %. The vertical error bars are not resolved in
  the photograph; no uncertainty column is given.
- The full digitization record (perspective correction, tick calibration, residuals, check images, the
  digitizing script) stays in `validation_runs/i3_master_campaign3_oxygen_20261119/report_oxygen/observations/`.
  It is not copied: it holds the photograph's local path and camera metadata, and the photograph itself
  is not kept (it shows the participants of a video call).
