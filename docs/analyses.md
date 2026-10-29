# Main analyses

This document maps the analyses described in the manuscript to the original scripts used to generate them. Some upstream scripts require additional inputs, software or configuration and have not been fully revalidated in the current repository structure. For reproducible figure generation from released data, see `release_analysis/figures/`.

These scripts are provided to document the computational workflow used in the study; they are not all intended as independently supported command-line tools.

## Sequence read processing

- [Trimming](../analysis/process_sequence_data/pipelines/trim.snakefile) and [Castanet](../analysis/process_sequence_data/pipelines/castanet.snakefile) rules are combined in the [read-processing Snakefile](../analysis/process_sequence_data/Snakefile).
- [Castanet finalisation](../analysis/figures/scripts/finalise_castanet_results.py) contains batch-specific filtering of results.
- [Index-hopping comparisons](../analysis/figures/scripts/index_hopping_by_batch.py) assesses the batch specific effects.

## Assessing enrichment of targets

### Coverage of mitochondrial genome

- [Alignment and coverage computation](../analysis/enrichment/pipelines/map_to_genome.snakefile) performed by a snakemake file.
- [Coverage collation and summaries](../analysis/enrichment/scripts/assess_enrichment/process.py).
- [Mitochondrial coverage plotting](../analysis/figures/scripts/enrichment_and_estimation_figure.py) shown in Figure 2A by the function `plot_mitochondrial_coverage`.

### Coverage of virus genomes

- [Castanet](../analysis/process_sequence_data/pipelines/castanet.snakefile) generates virus results.
- [Coverage extraction](../analysis/virus_detection/scripts/process.py), function `build_virus_coverage_plotting_data`, prepares positional depth.
- [Coverage plotting](../analysis/figures/scripts/enrichment_and_estimation_figure.py) uses `plot_virus_coverage_paneled`.
- [Ct–coverage modelling](../analysis/figures/scripts/virus_scatters_and_heatmap.py) uses `plot_nc2_curve_panels`.

### Modelling Ct vs COSMOS read count

- [`plot_virus_scatter` and `bfv_supp_figure`](../analysis/figures/scripts/virus_scatters_and_heatmap.py) fit the surveillance and BFV dilution relationships (Figure 3A and S3A).
- [Shared regression functions](../analysis/common/regression.py) calculate OLS fits and confidence intervals.

## Validation of mosquito species abundance estimation

- The separate repo [COSMOS-abundance](https://github.com/MultipathogenGenomics/COSMOS-abundance) implements species assignment, reference preparation and abundance estimation.
- [Control processing](../analysis/validation/scripts/process.py) prepares validation summaries.
- [Diagnostic-site analysis](../analysis/validation/scripts/check_pipiens_australicus_pools.py) assesses the morphologically similar Culex pools shown in Figure S4.

### Concordance correlation coefficient

- [Surveillance validation](../analysis/validation/scripts/surveillance_validation/process.py) matches surveillance counts to COSMOS samples and calculates Pearson correlation and CCC in `compute_correlation_values` (Figure 3C).

### Tests using in-silico mixtures of reads

- [Read extraction](../analysis/validation/scripts/in_silico_mixtures/get_mosquito_reads.py).
- [Mixture generation](../analysis/validation/scripts/in_silico_mixtures/mix.py).
- [KMA/Kallisto assignment](../analysis/validation/scripts/in_silico_mixtures/estimate_abundance.py).
- [Error calculation](../analysis/validation/scripts/in_silico_mixtures/abundance_error.py).
- [Plotting](../analysis/validation/scripts/in_silico_mixtures/plot.py).

## Mosquito community composition

### Computing Simpson's diversity index

- [`compute_simpsons_diversity_index`](../analysis/mosquito_community_composition/scripts/process.py) calculates diversity.
- [Diversity modelling](../analysis/mosquito_community_composition/scripts/PCA_loadings_with_virus.py) uses `plot_box_plot_of_location_type_vs_D1_with_regline` (Figure S7).

### Computing Principal Components on mosquito abundance

- [Input preparation](../analysis/mosquito_community_composition/scripts/linregress_abundance_viruses_coast.py), function `prepare_community_analysis_data`, includes lineage components.
- [PCA code](../analysis/mosquito_community_composition/scripts/PCA_loadings_with_virus.py) applies `hellinger_transform` and `run_PCA_on_mosquito_abundance`.

### Assessing the correlates of PCs

- [`associate_PCs_with_variables`](../analysis/mosquito_community_composition/scripts/PCA_loadings_with_virus.py) fits site-random-intercept models and applies Benjamini–Hochberg correction.
- [Figure 5](../release_analysis/figures/figure5.py) plots the released results.

## Cryptic population structure of Culex annulirostris

- [Diagnostic SNP/kmer discovery](../analysis/cryptic_mosquito_population_structure/Cx_annulirostris/scripts/determine_subtypes/) and [lineage-kmer counting](../analysis/cryptic_mosquito_population_structure/Cx_annulirostris/scripts/count_subtype_kmers/).
- [Allele counts and proportions](../analysis/cryptic_mosquito_population_structure/Cx_annulirostris/scripts/compute_subtype_snp_ratios/), [read-pair linkage](../analysis/cryptic_mosquito_population_structure/Cx_annulirostris/scripts/check_reads_for_linkage/) and [sample summaries](../analysis/cryptic_mosquito_population_structure/Cx_annulirostris/scripts/collate_results/).
- [Coastal associations and supplementary plots](../analysis/cryptic_mosquito_population_structure/Cx_annulirostris/scripts/plotting/); [map grouping](../analysis/common/annu_plot_groups.py) and [main-text map](../analysis/figures/scripts/results_figure.py).

## Virus phylogenies

- To be added in a later version.

# Supplementary analyses

## Insect-specific viruses

- Initial assembly and Kraken analysis: to be added in a later version.
- [`wuhan6_supp_figure`](../analysis/figures/scripts/virus_scatters_and_heatmap.py) handles downstream Castanet summaries, the virus heatmap and Wuhan Mosquito Virus 6 association.

## Cryptic population structure in Aedes notoscriptus

- [Sequence clustering](../analysis/cryptic_mosquito_population_structure/Ae_notoscriptus/scripts/tree/).
- [Diagnostic sites](../analysis/cryptic_mosquito_population_structure/Ae_notoscriptus/scripts/subtypes/).
- [Allele proportions](../analysis/cryptic_mosquito_population_structure/Ae_notoscriptus/scripts/compute_subtype_snp_ratios/), reusing Cx. annulirostris counting functions.
- [Linkage](../analysis/cryptic_mosquito_population_structure/Ae_notoscriptus/scripts/check_reads_for_linkage/), reusing Cx. annulirostris linkage functions.
- [Collation](../analysis/cryptic_mosquito_population_structure/Ae_notoscriptus/scripts/collate_results/).
- [Shared plotting code](../analysis/cryptic_mosquito_population_structure/Cx_annulirostris/scripts/plotting/plot.py) includes the latitude associations.
