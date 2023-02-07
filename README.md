Single-cell heterogeneity of EGFR and CKD4 co-amplification is linked to immune infiltration in glioblastoma
--------

This repository contains code, data, tables and plots to support analyses and reproduce results from the paper [Single-cell heterogeneity of EGFR and CKD4 co-amplification is linked to immune infiltration in glioblastoma](https://), currently provisionally accepted for publication as of February 2023.


Abstract
--------
Glioblastoma (GBM) is the most aggressive brain tumor with a median survival of ~15 months. Targeted approaches have not been successful in this tumor type due to the large extent of intratumor heterogeneity. Mosaic amplification of oncogenes suggests that multiple genetically distinct clones are present in each tumor. To uncover the relationships between genetically diverse subpopulations of GBM cells and their native tumor microenvironment, we employ highly multiplexed spatial protein profiling, coupled with single-cell spatial mapping of fluorescence in situ hybridization (FISH) for EGFR, CDK4, and PDGFRA. Single-cell FISH analysis of a total of 35,843 single nuclei reveals that tumors in which amplifications of EGFR and CDK4 more frequently co-occur in the same cell exhibit higher infiltration of CD163+ immunosuppressive macrophages. Our results suggest that high throughput assessment of genomic alterations at the single cell level could provide a measure for predicting the immune state of GBM.


Content
-------
* `/code/`: R scripts to reproduce analyses of single cell RNA sequencing data
* `/data/`: path folder to read in cellranger outputs, as well as preprocessed RDS data objects to be used as input for the scripts in this repo. This data needs to be downloaded from [Zenodo](https://zenodo.org/record/7454907#.Y97iXS-B2NF), and placed in the folder 'data'
* `/plots/`: output plots from the single cell RNA sequencing analyses, organized by self-explanatory folder names.
* `/tables/`: relevant gene lists (e.g. signatures and marker genes) used throughout the analyses, together with coordinates of cells after dimensionality reduction


Data
-------
The single cell RNA sequencing data generated in this study has been deposited in the NCBI GEO database under accession number [GSE203536](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203536).

Preprocessed RDS objects can be downloaded from [Zenodo](https://zenodo.org/record/7454907#.Y97iXS-B2NF).


## Contact
For any questions, comments or suggestions regarding the analyses or interpretation of the single cell RNA sequencing data in this manuscript, please feel to reach out to me via [GitHub](https://github.com/csimona), [Twitter](https://twitter.com/simocristea) or email (scristea@jimmy.harvard.edu).