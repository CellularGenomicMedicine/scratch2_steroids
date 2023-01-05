# Scratch2_steroids

In this project we have compared the levels of steroids as well as steroid-metabolizing genes in the endometrium of 34 women who did and did not become pregnant after a first failed in vitro fertilization (IVF) procedure. Steroid concentrations in endometrium biopsies were measured with liquid chromatography-mass spectrometry. RNA was also extracted from the biopsies and paired-end sequenced (2x50bp) on a NovaSeq 6000 S2 flow cell (Illumina). 

The align_count_workflow directory contains a snakemake workflow for the mapping to the reference genome and counting and the scripts directory contains scripts used for data analysis and vizualization. 

## Data analysis & visualization:

1. Combining expression count matrix for each sample in one joint expression count matrix
2. Filtering of the expression matrix: 
   - removal of rows with all 0s and values >=5% in >75% samples
