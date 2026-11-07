#!/bin/bash
# Full replication, in order. Run from the repository root after filling data/ (see data/README.md).
set -e
mkdir -p out
cd src
python3 convert.py ../data/theta_v2_25-04.csv ../data/theta
python3 convert.py ../data/SST_v2_25-04.csv ../data/SST
python3 prep_sat.py ../data
python3 build_mask.py ../data/global_topography.mat ../data/mask_mean300.npy
python3 get_indices.py ../data/oni.data ../data/pdo.dat ../data
Rscript fig1_rf.R ../data ../out
Rscript fig1_pdp.R ../data ../out
python3 sat_analysis.py ../data ../out ../data/mask_mean300.npy
python3 figures.py
python3 figures_si.py
# Optional, needs files not redistributed here:
# python3 fig1D_partial.py ../data/abe6301-ustick-sm-data-s1.csv ../data/data_insitu_v2_25-04.csv ../out
# python3 woa_steps.py ../data ../out        # needs WOA 2023 nitrate; not yet run (NOAA server down)
