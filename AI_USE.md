# AI Usage Documentation

## Tool and Scope

- AI tool: Claude
- Model: Opus 5.5

Claude assisted with paper replication and coding, including generating scripts to process the original study’s data and reproduce its analysis and figures. Claude also helped with result validation.

## Summary of AI-assisted tasks:

We used the original paper, the datasets used by the authors (listed at the end of the replicated paper), and their publicly available code as the basis for replication. Claude helped generate and revise scripts to process the data and reproduce the analyses and figures. After Claude generated a script, team members reviewed and ran the script. Output results (both numerical and images) were compared with the published results in the original paper, and discrepancies were discussed. Remaining differences and data availability limitations are documented in the report.
Specifically, Claude uses generated the replication scripts for Figures 1 - 4 using the original paper, the authors’ processed datasets released on Zenodo, and their publicly available code. This included code for data processing, statistical analyses, and figure generation. The data includes:
- Gridded satellite-derived nutrient stress scores (Θ′) and sea surface temperature, at 1° spatial resolution and 8-day intervals.
- A table of 1,137 metagenomic samples containing sampling locations and dates, nutrient stress scores (Ω), Prochlorococcus core-gene coverage, and ship-based measurements.
- Satellite Θ′ values matched to metagenomic samples across six spatial areas and 15 time windows.
- NOAA Oceanic Niño Index (ONI) and Pacific Decadal Oscillation (PDO) index records.
Claude assisted with data processing, statistical analyses, and figure generation. Follow-up prompts were used to revise the scripts, address errors, and investigate differences between reproduced and published results. Remaining discrepancies and data availability limitations are documented in the report.

## Prompt Log and Code Validation

Below documents some important prompts used for Claude:

### Prompt 1:
I want to reproduce the computational results of the attached paper for a student capstone project. First, help me assess feasibility of this paper replication and choose a manageable starting scope.
Please:

1. Briefly explain the paper’s main question and computational workflow.
2. Identify the required data, code, dependencies, and computing resources. Check availability if you have browsing access; otherwise, tell me what links or files you need.
3. Recommend one main figure, table, or experiment to reproduce first, and explain why.
4. Outline the steps from raw data to that result, distinguishing details explicitly stated in the paper from missing details or assumptions.
5. Explain how we would compare our output with the published result.

Start with this assessment before implementing the analysis. Ask me only for information that you cannot obtain from the provided materials. Never invent data, results, or successful execution; clearly distinguish proposed code from code you have actually run.

### Prompt 2:
https://www.science.org/doi/10.1126/sciadv.aed8089?referrer=https%3A%2F%2Fwww.google.com%2F Replicate this research paper as thoroughly and accurately as possible. Follow all the methods, procedures, and experiments described in the paper.
If you encounter any questions, uncertainties, or difficulties during the replication process, document them clearly, but continue working whenever possible. Only stop if you encounter a problem that makes further replication impossible.
Your goal is to replicate everything described in the paper to the greatest extent possible. Do not make up any information, results, methods, or assumptions. If something is unclear or missing, explicitly point it out rather than guessing.

### Prompt 3:
1. data files: I have attached `data_insitu_v2_25-04.csv and `the six `genomic_OCmatchups_v3_medians_*.csv` files ending in 1p, 3p, 5p, 1d, 2d, 5d; please verify if those are the correct files you need. And let me know if you need other files
2. MATLAB: I have matlab locally on my laptop for later analysis. For the later MATLAB analyses, approximately how much disk space will I need for the downloaded data, intermediate files, and outputs? Please distinguish between using the authors’ processed data and reprocessing the raw satellite data.

### Prompt 4:
Compare our reproduced Figures 4 A-C and their numerical results with the original paper. Identify which panels and statistics agree, which differ, and which could not be reproduced. Report the size of numerical discrepancies where possible.


