# Experimental warming restructures soil microbial and plant community composition in subalpine meadows

## Overview

This repository contains the data and R code supporting the manuscript:

> Experimental warming restructures soil microbial and plant community composition in subalpine meadows

Authors: Lorah Seltzer, Patrick O. Sorensen, Eoin L. Brodie, Aud Halbritter,
Ulas Karaoz, Jocelyn Navarro, Richard Telford, Alex Zafiro, Vigdis Vandvik, and
Brian J. Enquist

Study based on a whole community turf transplant experiment established in 2017
in the subalpine meadows of the Elk Mountains near Crested Butte, Colorado, USA
(near the Rocky Mountain Biological Laboratory or RMBL). This study examines how
plant, fungal, and bacterial/archaeal communities respond to environmental
change through whole community transplant of intact plant–soil communities
across a 400 m elevation gradient. The repository includes principal response
curve (PRC) analyses of plant, bacterial/archaeal, and fungal community
trajectories; comparisons of PRC taxon weights across the transplant experiment
and natural elevation gradient; microbial differential-abundance analyses;
non-metric multidimensional scaling (NMDS) analyses; and code that generates the
associated manuscript figures, statistical results, and publication tables.

## Relevant Methods

In 2017, we established three study sites along an elevational gradient at High
(“Monument”, 3300 m), Mid (“Pfeiler”, 3150 m), and Low (“Upper Montane”, 2900 m)
elevation in the subalpine meadows of the Elk Mountains near Crested Butte,
Colorado, USA. Each site location had a southern aspect and approximate 10%
slope. Mean growing season daytime air temperature was 2.6 degrees Celsius
higher and the gravimetric water content was 0.17 g H2O per g dry soil lower at
the Low compared to High site. Snowmelt, which initiates the start of the
growing season, occurred an average of 21 days earlier at the Low site than at
the High site, and 10 days earlier at the Mid than at the High site. Thus, the
Mid and Low sites were progressively warmer, drier and had an earlier snowmelt
date compared to the High site.

Each site was divided into five blocks of five 0.5 m x 0.5 m turfs arranged
perpendicularly to the hillslope (Figure 1). We focused on grassland
(grass/forb) vegetation and excluded shrubs when selecting turfs. Treatments
were assigned randomly to each turf within each block and included: two
‘untouched controls’ (UT), a control that was locally transplanted (LT) to
account for any effects of the transplanting, and two climate change
transplants, relocated either one or two sites higher or lower, with specific
climate change treatments varying between origin sites (W1, W2, C1, C2, Figure
1a). Each transplanted turf was moved intact and contained the entire plant and
soil community and an average depth of 23 cm (SD = 4.8 cm, Figure 1b) of intact
soil. Immediately after the transplant, turfs were watered from local streams
and covered by a shade-cloth canopy for one month. Fences were built around the
sites to exclude cattle grazing and trampling. In anticipation of natural
disturbance from gophers, a sixth block of transplanted turfs (5 at the High
site, 4 at the Mid site, and 6 at the Low site) was added in 2018, bringing the
total to 90 turfs. Nets were installed around each turf (1 m x 1 m) in the
growing season months to prevent cross-pollination from the transplanted plants
to the surrounding plant communities. One untouched control was left un-netted
in each block to test for the effect of nets.

Soil sampling was conducted in three control turfs and three transplanted turfs
at each of the sites. A total of 81 soil cores (diameter = 2 cm, depth = 5 cm)
were collected at three seasonal time points: late in the growing season in
2018, immediately after spring snowmelt in 2019, and at peak summer plant
biomass in 2019 (n = 9 samples per site and time point). Total soil DNA was
extracted from all soil samples in duplicate using the Qiagen PowerSoil DNA
extraction kit. During this process, 10 samples were mislabeled and treatments
could not be identified; thus, 71 samples remained. In order to assess responses
of microbial taxa to the experiment, PCR amplification of the bacterial/archaeal
16S Small subunit ribosomal RNA (16S rRNA, V-4 region) and fungal internal
transcribed spacer region 1-2 (ITS1-2) biomarkers was used following previously
described PCR methods (Sorensen et al. 2020). The PCR amplicons were pooled and
sequenced using MiSeq V3 300-bp paired-end sequencing. Microbial sequences were
given taxonomic assignments by comparing sequence identity to sequence databases
for bacteria, archaea, and fungi. The data were then quality-filtered to include
only amplicon sequence variants (ASVs) with at least 30 read counts and in at
least three samples from the same origin site. Total copy number of bacterial
and archaeal 16S rRNA genes were quantified to measure microbial abundance and
to correct ASV relative abundance to ASV absolute abundance in each sample.
Quantitative PCR was performed following methods that have been described
previously (Wang et al. 2021). The absolute abundance of all bacterial and
archaeal taxa was calculated from the relative abundance of each ASV multiplied
by the total absolute abundance in the sample (from qPCR), then divided by the
16S copy number. 16S copy numbers were obtained from the ribosomal RNA operon
copy number database (rrnDB version 5.7, Stoddard et al. 2015). When no copy
number was available, the uncorrected absolute abundance was calculated.

## Repository structure

```text
.
├── data/                 Input data used by the analyses
├── output/
│   ├── figures/          Figures generated by the analysis scripts
│   └── tables/           Analysis results and manuscript-ready tables
├── 01_plant_PRC_analysis.R
├── 02_microbial_PRC_analysis.R
├── 03_cross_community_PRC_analysis.R
├── 04_microbial_diff_abund_analysis.R
├── 05_microbial_NMDS_analysis.R
├── 06_make_manuscript_tables.R
└── README.md
```

## Data files

### `data/plant_abundance_2017.csv`

Abundance data were collected in 2017 (pre-transplant). Baseline plant species
abundance was recorded by counting individuals in one turf (out of five) per
block at each site during peak plant biomass.

### `data/plant_community_cover_2018_2023.csv`

Species cover data were collected from 2018-2023 (excluding 2020) during peak
summer biomass on all turfs using a 50 x 50 cm frame with 10 x 10 cm subplots.
First, percent cover of each plant species and ground cover type over the entire
turf was visually estimated to the nearest 1%. Second, species presence was
noted in each of the 25 subplots. Each subplot was assigned its own number,
working left to right (Figure 1, 1-1 through 5-5), and every plant species and
ground cover type present in each subplot was recorded. Ground cover types were:
Rock (greater than 5cm wide), litter (unrooted, dead plant material), bare soil
(soil or small gravel <5cm with no plants), and moss. We noted for plants if
each was: sterile (S), fertile (F), sterile and dominant (SD), fertile and
dominant (FD), or seedling (SE). We noted for ground cover types if each was:
present (P) or dominant (PD). A plant species or ground cover type was
considered dominant if it covered the most horizontal and vertical space in the
subplot. Only one plant species or ground cover type was marked as dominant per
subplot. Plants were marked as present if they covered a subplot, even if they
were rooted in another subplot. Plants were only considered fertile if there
were reproductive parts in the subplot. Third, plant and ground cover presence
were recorded for the 5 neighboring subplots to the right side of the turf when
looking uphill. Neighbors represent the plant community native to the
destination site, closest to the transplanted turf.

### `data/plant_plot_metadata.csv`

Contains the unique plot ID (`turfID`), destination site, origin site, treatment
including control type, treatment and origin-site combination
(`treatmentOrigin`), treatment group (`treatmentOriginGroup`), and origin-site
block.

### `data/plant_taxonomy.csv`

Plant species and families found in turfs.

### `data/bacteria_archaea_asv_table.csv`

Quality-filtered amplicon sequence variant (ASV) data for bacteria/archaea.

### `data/fungal_asv_table.csv`

Quality-filtered amplicon sequence variant (ASV) data for fungi.

## Packages

The clean reproducibility run was completed with R 4.3.1 (2023-06-16) and the
following package versions:

- `vegan` 2.6.4
- `ggvegan` 0.1.999
- `tidyverse` 2.0.0
- `ggplot2` 3.5.1
- `cowplot` 1.1.3
- `ggpubr` 0.6.0
- `DESeq2` 1.42.0
- `ComplexHeatmap` 2.18.0
- `circlize` 0.4.15
- `broom` 1.0.5
- `devtools` 2.4.5
- `precrec` 0.14.4

The `grid` package is included with R.

## Running the analysis

Run the scripts from the repository root. Scripts use repository-relative paths,
so the working directory should be the top-level repository directory.

The recommended order is:

1. `03_cross_community_PRC_analysis.R`
2. `04_microbial_diff_abund_analysis.R`
3. `05_microbial_NMDS_analysis.R`
4. `06_make_manuscript_tables.R`

Script 03 sources scripts 01 and 02, so scripts 01 and 02 do not need to be run
separately when reproducing the complete workflow.

The analysis scripts set fixed random seeds before permutation tests and NMDS
random starts so that results are reproducible across runs with the documented
software versions.

## Script descriptions

### `01_plant_PRC_analysis.R`

Performs principal response curve (PRC) analysis of the plant community and
generates Figure 2.

### `02_microbial_PRC_analysis.R`

Performs principal response curve (PRC) analysis of the bacterial/archaeal and
fungal communities and generates Figure 3.

### `03_cross_community_PRC_analysis.R`

Runs scripts 01 and 02 and performs comparisons of PRC taxon weights across
communities.

### `04_microbial_diff_abund_analysis.R`

Performs microbial differential-abundance analysis and generates Figure 4.

### `05_microbial_NMDS_analysis.R`

Performs microbial non-metric multidimensional scaling and generates Figures S1
and S2.

### `06_make_manuscript_tables.R`

Converts complete analysis outputs into manuscript-ready Tables 1–3 and S3–S6.

## Outputs

Generated figures are written to `output/figures/`. Generated statistical
results and manuscript-ready tables are written to `output/tables/` and
`output/tables/manuscript/`.

## Data and metadata

The complete data, detailed metadata, and supporting documentation are available
through the [OSF
project](https://osf.io/9hkc7/overview?view_only=b04a695af1004c0c8c25a259b12db1eb).

## Contact
Lorah Seltzer: [lorahseltzer@gmail.com], Patrick Sorensen:
[patrick.sorensen@uri.edu], Brian Enquist: [benquist@arizona.edu]
