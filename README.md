# COSMOS-paper

Code accompanying the manuscript *COSMOS: integrated vector-pathogen metagenomics for mosquito-borne disease surveillance*, by Forster et al.

This repository contains two complementary sets of scripts:

- [`release_analysis/`](release_analysis/) reproduces the computational panels of Figures 2–5 from the archived data release. Figure 4 includes panels A and B; virus phylogenies are generated separately. Sequencing reads are not required to run these figure scripts.
- [`analysis/`](analysis/) contains the code used for the analyses presented in the paper. These scripts have not been fully tested in the current repository structure and may require adjustments to run. Further details of each analysis can be found in `docs/`.

Processed data and reference sequences are archived on [Zenodo](https://doi.org/10.5281/zenodo.23310135). Sequencing reads are deposited in the [European Nucleotide Archive, project PRJEB128300](https://www.ebi.ac.uk/ena/browser/view/PRJEB128300). Mosquito species assignment and relative abundance estimation are implemented in the separate [COSMOS-abundance repository](https://github.com/MultipathogenGenomics/COSMOS-abundance).

## Installation

Clone the repository and create the Conda environment:

```bash
git clone https://github.com/MultipathogenGenomics/COSMOS-paper.git
cd COSMOS-paper
conda env create -f environment.yml
conda activate COSMOS_paper
```

Run all commands below from the repository root. Python packages and command-line software available through Conda are listed in [`environment.yml`](environment.yml).

Additional software used by upstream analyses includes Castanet (v9.2/v9.3, depending on the run), Iorek `tabulate-alignments`, the Python 3 adaptation of ngsutils `removeclipping_py3.py`, and CliqueSNV for the optional haplotype analysis. These are not needed to reproduce Figures 2–5 from released data. Supply their executable or script locations in `analysis/common/software.py` when running the relevant upstream analyses. COSMOS-abundance has its own installation instructions.

## Download the data

Download and extract the [Zenodo data release](https://doi.org/10.5281/zenodo.23310135) into `data/`. The analysis directories must sit directly inside `data/`, without an additional enclosing archive directory:

```text
data/
├── samples/
├── enrichment/
├── virus_detections/
├── mosquito_species_calls/
├── controls_and_validation/
├── cryptic_mosquito_subpopulations/
├── PCA/
└── references/
```

The figure scripts use the released tables and reference data in this layout. Use the data-release version associated with the code revision being run. Downloading ENA sequencing reads is necessary only for upstream processing and abundance estimation.

Figure 4 also requires the [Australian Bureau of Statistics 2025 digital boundary GeoPackage](https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs/edition-3-july-2021-june-2026/access-and-downloads/digital-boundary-files/ASGS_Ed3_Non_ABS_Structures_GDA94_updated_2025.zip). Extract the ZIP file and retain the path to the `.gpkg` file. The globe inset (which in the final manuscript was moved to Figure 1) uses the Natural Earth land dataset through `geodatasets`; an internet connection may be needed to download it on first use.

## Reproduce the figures

```bash
python -m release_analysis.figures.figure2
python -m release_analysis.figures.figure3
python -m release_analysis.figures.figure4A_B \
    --map-file "/path/to/ASGS_Ed3_Non_ABS_Structures_GDA94_updated_2025.gpkg"
python -m release_analysis.figures.figure5
```

| Script | Output |
| --- | --- |
| `figure2` | Mitochondrial enrichment, viral depth profiles and the Ct–genome coverage relationship |
| `figure3` | Viral read recovery, ASMMP abundance comparisons and mosquito pool controls |
| `figure4A_B` | Community composition, virus detections and the Culex annulirostris lineage map (panels A and B) |
| `figure5` | PCA scores, loadings and environmental/temporal associations |

Each command writes a PNG and PDF to `results/`, named `figure2`, `figure3`, `figure4A-B` and `figure5`. Figures 2 and 3 also write statistical summaries alongside the figures. The scripts calculate regression lines from the released observations; precomputed regression lines are not required.

Use `--data-dir` to specify another release-data directory, `--output` to change the output filename, or `--help` to see script-specific options. For example:

```bash
python -m release_analysis.figures.figure3 --help
```

## Upstream processing

The read-processing workflow has an [example configuration](analysis/process_sequence_data/config.example.yaml). Copy and edit it to supply sample FASTQ paths, adapter sequences and any required software locations. The configuration comments describe how to invoke the parent [Snakefile](analysis/process_sequence_data/Snakefile).

Upstream analyses generate intermediate files under `results/`; they do not need to be rerun to draw the figures from released data.

Software path settings and expected upstream reference files are described in
[upstream inputs](analysis/UPSTREAM_INPUTS.md).

## Mosquito species abundance estimation

Follow the installation instructions in [COSMOS-abundance](https://github.com/MultipathogenGenomics/COSMOS-abundance). Use the released COX1 reference panel and its species mapping, and populate a sample manifest file with paths to the trimmed paired-end FASTQ files used for assignment. A template for the sample manifest is provided: `data/mosquito_species_calls/template_sample_manifest_for_COSMOS_abundance_repo.tsv`.

The sample sheet requires tab-separated columns named `sample`, `read1` and `read2`. Relative FASTQ paths are resolved against the directory containing the sample sheet. The example below assumes the species mapping has been placed at `data/mosquito_species_calls/species_mapping_for_use_with_COSMOS_abundance_repo.tsv` and the completed sample sheet at the repository root:

```bash
cosmos-abundance run \
    --fasta data/references/cox1_reference_panel.fasta \
    --species data/mosquito_species_calls/species_mapping_for_use_with_COSMOS_abundance_repo.tsv \
    --samples samples.tsv \
    --out results/abundance \
    --cores 8
```

The existing reference FASTA is used directly, without repeating reference preparation. Consult the separate repository for input formats, parameters and outputs. Record the COSMOS-abundance release or commit used for a run.
