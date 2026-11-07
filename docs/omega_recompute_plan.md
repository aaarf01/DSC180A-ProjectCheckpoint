# Recomputing the DNA stress index (Ω) from raw reads: plan and open issues

Status: preparation done, no raw reads downloaded yet. Waiting on network access and a go-ahead on size.

## Sources for the method
- Martiny et al. 2026 Methods: Ω "calculated as described previously (8)" = Ustick et al. 2021, Science 372:287.
- Ustick 2021 supplement (Materials and Methods, table S3, data S1, data S2), read through the user's browser.
- Authors' code: github.com/ljustick/metagenome-cluster-coverage (only sums coverage over the genes in a cluster).

## Pipeline as stated
1. Trimmomatic 0.35: quality and adapter trimming. **Settings not given.**
2. Bowtie2 2.2.7 against 115 reference genomes (table S3), flags `--no-unal --local -D 15 -R 2 -L 15 -N 1 --gbar 1 --mp 3`.
3. samtools 1.3 sort and index.
4. Anvi'o v5 profile; pangenome of reference ORFs by BLAST + MCL; pick target clusters (data S2) and single-copy core genes (SCCG).
5. Gene coverage summed over all members of a cluster (cluster_coverage.m).
6. f_i = gene coverage / sum over ecotype clades of (mean SCCG coverage in that clade).
7. Drop samples with total Pro SCCG coverage < 5X. If HNLC clade SCCG coverage < 1X, set Fe-high f_i to 0.
8. z_i = (f_i - mean) / sd across the whole dataset; Ω_s = mean z over genes in class s.

## Checks already done (no raw data needed)
- Steps 7-8 reproduce the published Ω exactly from published f_i (data S1), using the sample sd (ddof=1). Max difference 0.
- Ω in the 2026 Zenodo file is identical to Ω in the 2021 data S1 for all 1,137 samples. So checking Ω = checking the 2021 pipeline.
- Implication: z-scores depend on the mean and sd of all 1,137 samples. A subset cannot produce its own Ω on the same scale.
  Plan: compare recomputed f_i directly with published f_i, and place recomputed f_i on the published scale using the published mean and sd.

## Tool versions installed
Trimmomatic 0.35, Bowtie2 2.2.7 (binaries run from package folder; the conda build is mis-packaged), samtools 1.3. BLAST 2.15, MCL, Prodigal for clustering.
Anvi'o v5 not installed (2018 release); its gene-coverage step (mean depth across a gene) will be reimplemented. To be validated against published f_i.

## Open issues (not stated anywhere; each will be flagged in results)
- Trimmomatic parameters.
- Exact reference genome versions. Several strains have 2014 drafts and 2025 complete genomes; pre-2020 versions will be used and checked with data S2 coordinates and sequences.
- 14 of 117 listed names not found in public databases: Prochlorococcus RS50, UH18301; Pelagibacter A9spades, B4spades, F4spades, AAA298_D23, HIMB122, HTCC1016, HTCC7214, HTCC7217, HTCC9022, IMCC1322, PRT004; Synechococcus NKBG042902. MIT9313 is listed twice. SAR86B is listed under Pelagibacter.
- Gene caller used for ORFs (anvi'o v5 used Prodigal by default), BLAST/MCL settings (anvi'o defaults assumed: minbit 0.5, MCL inflation 2).
- Which genes counted as SCCGs, and how the 41 Prochlorococcus genomes were assigned to ecotype clades. Neither is published.
- How the HLIII/IV "Fe high" index was computed from HNLC clade coverage. Data S2 only says "HLIII-IV core genes".
- Anvi'o v5 default minimum contig length (1,000 bp) would drop genes on short contigs in draft genomes.
- Whether the small MiSeq runs (about 50,000 reads) in each BioSample were included. Plan: shotgun NovaSeq/HiSeq runs only.

## Fig. 1D
- Trait genes (2026 table S1) not in the 2021 gene list, so not computable from anything public: nitC, nitD, hypA, hypB, cynB, cynD, nirA2, narX, amt.
- "nirA" in 2021 data S2 is a 162 bp PAC1 gene classed N-medium and correlates -0.27 with Θ′. Published nitrite box is about +0.21 to +0.27.
  The 2026 paper likely uses nirA for a different gene. Unresolved.
- Partial Fig 1D from published f_i: omega/fig1D_published_fi.png, omega/fig1D_from_published_fi.csv.

## Sizes
- Proposed subset: 12 Bio-GO-SHIP samples, stratified by dominant stress type (Fe high, Fe medium, N high, P high), all with Θ′: 13.6 GB, 104 M read pairs (omega/subset.csv).
  Processed one at a time, so peak disk use is about 5 GB.
- Reference genomes: about 103 assemblies, roughly 70 MB compressed.
- All 908 Bio-GO-SHIP samples: about 1.6 TB. Not possible here (25 GB free disk, 2 CPU cores, 7 GB RAM).
- Runtime: unknown until benchmarked on the first sample.
