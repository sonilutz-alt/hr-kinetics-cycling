# hr-kinetics-cycling

Analysis code for the manuscript *"Enhancing Endurance Performance: Time-Series Monitoring of Autonomic HR Regulation in Professional Cycling"* (European Journal of Applied Physiology, manuscript EJAP-D-26-00784, revision 1).

The repository contains the full pipeline used in the revised manuscript, from race files to statistical output. **Athlete data are not included** (see *Data availability*).

## Pipeline

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `scripts/01_extract.py` | one `.xlsx` per race file (1 Hz timestamp, power, heart rate, exported from the device `.fit` file) in `data/raw/` | `output/BBDD_revision1_*.xlsx` (one row per race file × quartile) |
| 2 | `scripts/02_results_and_figure2.py` | output of step 1 | `output/Resultados_finales_*.xlsx` (Table 1, Table 2, Supplementary Tables S1–S4) and Figure 2 |
| 3 | `scripts/03_figure1.py` | one race file | Figure 1 (method illustration) |

Run from the repository root:

```bash
pip install -r requirements.txt
python scripts/01_extract.py
python scripts/02_results_and_figure2.py
python scripts/03_figure1.py
```

`data/rider_map.csv` (columns `archivo,corredor_id`) links each race file to an anonymous rider code. It is kept locally and is not published.

## What step 1 computes

- **Pre-processing.** A true 1-s time grid is built (duplicated timestamps averaged; missing seconds kept as missing). Power > 1800 W or < 0 W and HR < 30 or > 220 bpm are set to missing. Gaps of up to 5 s are linearly interpolated; longer gaps stay missing (no imputation). Signals are smoothed with a 2-s moving average.
- **Race quartiles.** Four equal time-based quartiles (Q1–Q4) of each race file.
- **Power-change segments.** The first derivative of smoothed power is computed. Consecutive seconds in which it exceeds +0.3 SD (power-increase segments) or falls below −0.3 SD (power-decrease segments) form a segment.
- **HR derivatives (bpm·s⁻¹).** For each segment, the mean first derivative of smoothed HR is computed. Increase segments with a negative mean and decrease segments with a positive mean are discarded, following the original definition. The per-quartile value is the mean over the retained segments. Versions without this sign filter are also exported for transparency.
- **Relative power (%CP).** Mean power over the power-increase segments of each quartile divided by critical power. CP is estimated per race file from the best 1-, 5- and 12-min mean powers of that file. No default CP value is substituted.
- **Maximal HR.** Highest smoothed HR within the power-increase segments of each quartile.
- **Rate of HR decrease** is reported as an absolute value (higher = faster).
- **Inclusion.** Files with physiologically impossible power values (any sample > 1800 W) are excluded. Quartiles with < 80% valid power or HR are flagged for sensitivity analysis but not excluded.
- **Validation.** Step 1 also re-runs a line-by-line Python port of the earlier R extraction script, so that the two can be compared.

## Statistics (step 2)

Linear mixed models with quartile, group (top-10 vs 11th–30th) and their interaction as fixed effects and a random intercept per race file (`statsmodels` MixedLM, REML), fitted as cell-mean models. Main effects and the interaction are tested with F tests whose denominator degrees of freedom are assigned by level, as in `nlme::anova` (estimates and F statistics were checked against `nlme`). Partial eta-squared is derived from F. Pairwise quartile comparisons are Tukey-adjusted. Group differences within quartiles are Holm-adjusted, and Hedges' g is computed from the observed data. The rate of HR decrease is analysed as an absolute value. Sensitivity analyses: race as covariate, exclusion of Milan–San Remo, exclusion of quartiles with < 80% valid data, and finishing position (1–30) as a continuous predictor. A sensitivity power analysis is also run.

## Data availability

The race files belong to professional riders and were provided under confidentiality. They cannot be shared publicly. Aggregated, anonymised per-quartile data may be requested from the corresponding author.

## Citation

Please cite the article (reference to be added upon publication) and this repository (Zenodo DOI to be added upon release).

## License

MIT (see `LICENSE`).
