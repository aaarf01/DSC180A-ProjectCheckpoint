# Replicating "Genomic-to-space measurements reveal large-scale ocean nutrient stress"

DSC 180A Quarter 1 project. This repository replicates Martiny et al. (2026), *Science Advances* 12, eaed8089
([paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13240183/)), starting from the authors' posted data
(Zenodo [15168108](https://zenodo.org/records/15168108)) and their code
([ljustick/genomic_to_space_nut_stress](https://github.com/ljustick/genomic_to_space_nut_stress)).

The authors' MATLAB steps are ported to Python. The random forest steps use R, as in the paper.
We start from the authors' computed Θ′ values and DNA stress index (Ω). We did not reprocess raw satellite data or metagenomes.

## Main results

| Result | Paper | This replication |
|---|---|---|
| Random forest R², nutricline + genes (Fig. 1B) | 59.94% | 60.2% |
| Random forest R², nutricline only (Fig. 1B) | about 0.36 | 0.34 |
| All 90 Fig. S4 cells | printed values | r = 0.999 |
| Top 4 EOF modes | 59.4% | 59.8% |
| EOF modes 5-8 | 4.2% | 4.2% |
| Median Θ′ trend per year | -1.3e-3 | -1.39e-3 |
| Median SST trend (°C per year) | 2.3e-2 | 2.29e-2 |
| Share of area warming | 89% | 88% |
| Seasonal ΔSST vs ΔΘ′ (Fig. S16) | -0.51 | -0.51 |
| Annual ΔSST vs ΔΘ′ (Fig. S16) | -0.23 | -0.22 |
| PC5 vs ONI, PC6 vs PDO | no r given | \|r\| = 0.84, 0.60 |

Not reproduced yet: Fig. 1A global r = -0.52 and Fig. 2B (World Ocean Atlas server was down), Fig. 1D (only partly; some genes are not public),
Fig. S3, S6 and S1. Fig. 4D uses a station-based stand-in for the nutricline. Full list in `docs/REPLICATION_STATUS.md`.

## Problems found in the paper and code
- `calc_thetaprime_v2.m` has a typo (`chlCSM`) and cannot run as posted.
- In `random_forest_analysis_v2.R` the 1-degree file is used in the slot meant for the 1-pixel file, so the "1P" column in Fig. S4 is mislabelled.
- The supplement gives a 1:4 ratio of long-term to seasonal variance; the paper's own numbers give about 1:14.
- Several numbers differ between sections (sample size 990 vs 1,137; matching timescale; smoothing window 20 vs 30).

## Layout
```
src/        analysis code (Python and R)
tests/      unit tests for the ported MATLAB functions
results/    output tables
figures/    rebuilt figures
docs/       replication status, discrepancy log, plan for recomputing Ω from raw reads
data/       input data, not committed (see data/README.md)
run.sh      runs everything in order
```

## How to run
1. Download the files listed in `data/README.md` into `data/`.
2. Python 3 with the packages in `requirements.txt`. R 4.2 with `randomForest` 4.7-1.1 and `pdp`
   (we used R 4.2.3; the paper used R 4.2.1).
3. `bash run.sh` from the repository root. Outputs go to `out/`.
4. `python3 tests/test_lib.py` runs the unit tests.

## Choices we made that the paper does not specify
- The authors' ocean mask (`BathyLL.mat`) is not public. We use ETOPO5 mean cell depth below 300 m, which gives 20,828 grid points (paper: 20,822).
- The SST preprocessing script is not public. SST gets the same steps as Θ′, without the 0-4 limit.
- ONI is NOAA PSL's current version (ERSST v6); the paper used v5.
- MATLAB `bootci` defaults to BCa intervals; we use percentile intervals.

## AI and Human Contributions

- AI: Claude (Opus 5.5) generated scripts to process the original study’s data and reproduced its analysis and figures, and assisted with debugging and revising the replication code.

- Humans: Team members selected replication scope/procedure, reviewed and ran the code, compared outputs from Claude code with the published results, and investigated and interpreted discrepancies; responsible for topics, data sources, and content of the written report.

- Testing: Validation included comparing reproduced figures and numerical results with those reported in the paper. Most reproduced results generally aligned with the published findings, while remaining discrepancies and data availability limitations are documented in the report.
