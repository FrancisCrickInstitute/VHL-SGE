# VHL-SGE
Code for processing of SGE data sets and all analyses included in "Saturation Genome Editing Resolves the Functional Spectrum of Pathogenic VHL Alleles".

Please see the manuscript for a description of the data sets and analyses performed. This depository's purpose is to make available all scripts and files used to analyse SGE data for VHL.

The SGE analysis pipeline in contained in:  **SGE_main_analysis**. The order of scripts used in SGE data analysis data is:

1.  **20230408_VHL_final_pipeline.sh** - A custom pipeline for converting paired-end, Illumina-based sequencing data into variant counts, extracting editing outcomes, and annotating variants by genomic position using CADD data.

2. **230415_VHL_all_no_plots.Rmd** - Processing all SGE regions anaylzed by replicate and analyzing replicates together in R, leading to a single dataframe for each set of replicates (i.e. the same SGE region in the same cell line performed under the same experimental conditions).

3.  **VHL_global_20240419.Rmd** - Merging all data together, applying thresholds to exclude variants not reliably detected, importing additional data sets for analysis (ClinVar, cBioPortal, GeneBass, computational predictors, etc), and generating all SGE data figures. Output includes Supplementary Table 1 which contains all SGE "function scores" and all "RNA scores" derived.

Several input files are required to run the .sh pipeline. These are found in **sh_input_files**. Additionally, the pipeline calls custom scripts written in Python (v2.7.16), which are located in **custom_scripts_for_SGE_pipeline**.

Inputs for VHL_global_20230504.Rmd are included in **rmd_input_files**. Data sets from VHLdb and cBioPortal were parsed using custom python scripts, available in **scripts_to_parse_external_data** with required inputs (outputs are included in rmd_input_files). 

Separate folders contain:
1.)  DepMap data and R code to produce the analysis of DepMap Cell lines:  **DepMap_analysis**.
2.)  The inputs and code required for generating VHL SGE oligo pools for each library:  **Mutagenize**.

