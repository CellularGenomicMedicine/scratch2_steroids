# Scratch2_steroids

In this project we have compared the levels of steroids as well as steroid-metabolizing genes in the endometrium of 34 women who did and did not become pregnant after a first failed in vitro fertilization (IVF) procedure. Steroid concentrations in endometrium biopsies were measured with liquid chromatography-mass spectrometry. RNA was also extracted from the biopsies and paired-end sequenced (2x50bp) on a NovaSeq 6000 S2 flow cell (Illumina). 

The align_count_workflow directory contains a snakemake workflow for the mapping to the reference genome and counting and the scripts directory contains scripts used for data analysis and vizualisation. 

## Data analysis & visualization:

1. [Combining expression count matrix for each sample in one joint expression count matrix](https://github.com/darinaobukhova/scratch2_steroids/blob/main/scripts/01_count_matrix_generation.R).
2. [Expression count matrix processing] (https://github.com/darinaobukhova/scratch2_steroids/blob/main/scripts/02_count_matrix_and_metadata_preprocessing.R). 
   - removal of rows with all 0s and values >=5% in >75% samples
   - retrieval of gene names in the place of Ensembl identifiers
   - correcting for LH day using ComBat_seq() from sva package
3. Principal component analysis
4. Differential gene expression analysis using DESeq2 package
5. Performing descriptive statistical analysis and visualisation of tissue steroid levels as violin plots
6. Analysis and visualisation of raw and variance-stabilized steroid-metabolizing gene expression profiles as boxplots
   
