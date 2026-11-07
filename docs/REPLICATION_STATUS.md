# Replication status: Martiny et al. 2026, Sci. Adv. 12, eaed8089

**State: the main quantitative results reproduce closely from the authors' posted data. Fig. 1A, 2B, 1D, S1, S3 and S6
are not reproduced, Fig. 4D only by proxy, Fig. 1C only partly. Raw reprocessing was not attempted.**

Correction: an earlier version of this file said the published Fig. S4 did not come from the posted code. That was wrong. See D2.

## What was read
- Main text and supplement (full text, PMC13240183). Equations 1-4 are images and did not extract, so they
  are taken from calc_thetaprime_v2.m, not read off the paper.
- Authors' code: github.com/ljustick/genomic_to_space_nut_stress (cloned), same project as Zenodo 15168175.
- File lists for Zenodo 15168108 and Dryad 10.5061/dryad.f4qrfj71p. Contents not seen.

## Environment
Python port of the MATLAB steps (MATLAB is not available). R 4.2.3 (paper: 4.2.1), randomForest 4.7.1.1 (same).

## Data used
Zenodo 15168108 (all files), ONI (PSL, now ERSST v6; paper used v5), PDO (NCEI ERSST v5), ETOPO5 from CDT.
Blocked: WOA 2023 (NOAA server 503/timeouts). Dryad .mat files need a login (401).
Not public: BathyLL.mat, Dataset S1_WOD.csv, gene abundance table, Fig 1C code.

## Targets from the paper
| Result | Paper value | Script |
|---|---|---|
| Fig 1A nutricline vs theta (all grid boxes) | r = -0.52 | needs WOA |
| Fig S4 R2 at 2x2 deg, +-20 d | 59.94 (read from the figure's embedded numbers) | fig1_rf.R |
| Fig S4 full grid | 15x6 values, in the supplement | fig1_rf.R |
| Fig 1B | genes nearly double R2 over nutricline alone | fig1_rf.R |
| Seasonal variance | 59% | run_all.py |
| EOF PC1-4 / PC5-8 / PC1-8 | 59.4 / 4.2 / 63.6 % | run_all.py |
| EOF SST PC1-8 | 89.7 % | run_all.py |
| SST warming area | 89 % | run_all.py |
| Median SST trend | 2.3e-2 C/yr (0.43 C total) | run_all.py |
| Median theta trend | -1.3e-3 /yr | run_all.py |
| Fig S16 seasonal / annual r | -0.51 / -0.23 | to add once data is in |
| Sensor agreement | 0.98 for means | to add, needs Dryad |
| Temperature effect on theta | <1% of variance | to add |

## Discrepancies found (paper vs code vs supplement)
- D1 calc_thetaprime_v2.m has `CChlDM.*chlCSM`, an undefined name. It must mean CChlSM, so the script cannot run as posted.
- D2 random_forest_analysis_v2.R loads the 1p file but lists sat_data_1d in slot 1, so the 1P column is really the 1-degree file run
  a second time. The published 1P column matches my 1-degree column (r=0.999) and not the true 1-pixel file (r=0.32), mean gap 0.73 vs 29.8.
  So published Fig S4's 1P column is mislabelled. True 1-pixel R2 stays below about 40% with a large drop from +-14 to +-24 days.
  Published 1P and 1x1D differ by RF noise, so the figure used another seed or order than the posted set.seed(123).
- D3 Sample size: n = 990 in Methods; n = 1137 in Fig S2.
- D4 Time scale: ~40 days in the main text, ~3 weeks in the supplement text, ~1 month in Fig S4 caption.
- D5 Lowess window: 20 in Methods, 30 in climate_index_comparison.m. run_all.py does both.
- D6 EOF top eight: 60% in supplement text, 63.6% in Fig S12 caption (59.4 + 4.2 = 63.6).
- D7 Sensor correlation: Methods 0.74 / 0.97; Fig S6 0.75 / 0.98; main text 0.98.
- D8 Methods cite fig. S15 for the RF grid; that is the SST figure. It should be fig. S4.
- D9 Grid points: 20,822 in Methods, 22,468 in a code comment. Fig S16 N values imply 29,520 total cells (82x360).
- D10 Fig 1C has six bars incl. medium P, but the stated RF has five biomarkers (no medium P). Fig 1C code is not posted.
- D11 Repo README cites older Zenodo DOIs (8064615, 8065231); the paper cites 15168108 and 15168175.

- D12 CSV lat/lon label columns are recycled vectors, not true pairs per row. Reshaping 180x360 column-first gives a correct world map; labels ignored.
- D13 SI says the ratio of long-term to seasonal variance is 1:4 for both theta and SST. Paper's own 4.2/59.4 is 1:14; I get 1:14 and 1:47.
- D14 "59% of variance" is not defined. ANOVA median gives 54.7%, EOF1-4 gives 59.8% (paper 59.4), so it likely means the EOF number.
- D15 Fig S16 N (354,240; 590,400) equal 29,520 cells x 12 or 20, so they count NaN cells too.
- D16 Time axis ends with one composite dated 2022 day 1, and 2002 day 217 is missing; 897 steps. Paper says 20 years.

## Assumptions in the port (check before trusting close results)
- A1 BathyLL.mat missing. Substitute: ETOPO5 mean cell depth >300 m, threshold picked to get close to 20,822 (gives 20,828).
- A2 Latitude rows 50:131 taken as |lat| < 41.
- A3 bbp443 to bbp532 power law form is assumed; only matters if theta is recomputed.
- A4 CDT eof checked against CDT source: unweighted, matches. Sign is arbitrary.
- A5 inpaint_nans(.,1) on a vector treated as linear fill.
- A6 MATLAB lowess: even spans reduced to odd, end handling from memory of smooth.m; not checked against MATLAB.
- A7 bootci in MATLAB defaults to BCa; port uses percentile intervals.
- A8 Cell area uses cos(lat), CDT also varies radius with latitude. Small.
- A9 DataLL_SST script not posted. SST got the same steps as theta without the 0-4 clip. Likely why SST EOF is 91.9 vs 89.7.
- A10 The single 2022 composite was joined to 2021 so there are 20 year levels.
- A11 Fig 1C bar definition is not stated; results depend on my choice.
- A12 Fig 4D proxy nutricline: mean of 8 nearest stations, cells within 10 degrees only (13,586 of 20,828). Not the paper's method.

