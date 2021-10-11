<img src="../docs/logo_1_transparent.png" width="75">

# deepC - Tutorials
Tutorials for the basic deepC workflows.

-------------------------------------------------------------------------------

### Description

These tutorials follow the workflow we use for deepC predictions, analysis and visualizations.

The deep learning framework and the deepC predictions itself are run via python scripts.

The downstream analysis and visualizations is done within R using tidyverse and custom functions for working with deepC predictions and HiC data in general.

The tutorials use the *current_version* scripts for tensorflow 1.14. If you are running on 1.8 just switch the script links to the legacy version.

In addition, we wrote wrapper scripts that perform standard preprocessing and downstream visualization tasks for deepC from the command line. These are written in R, software and data requirements are outlined below.

## Access

Checkout the html files from your repo clone or download the html files. It has some plots that lead to a bigger file size so gitHub can't visualize it straight online. Alternatively run the R markdown files within e.g. Rstudio.


### Tutorials

1) [Predicting and plotting chromatin interactions](./tutorial_predict_and_plot.html)

2) [Formatting HiC data to skeleton data for deepC](./tutorial_format_HiC_data_for_deepC.html)

3) [Train a deepC network](./tutorial_train_a_model.md)


### Wrapper scripts

#### 1) Wrapper for Hi-C data pre-processing and formatting

##### Requirements

* recent version of R with the following libraries:
  * optparse
  * tidyverse
  * cowplot
  * Rcolorbrewer
* working version of perl installed and accessible for command line calls within R.

###### Inputs

* a chromosome sizes file (two columns chrom size) indicating the chromosome
size for every chromosome named in the Hi-C data used. For example, as downloaded from the [UCSC Sequence and Annotation Downloads](https://hgdownload.soe.ucsc.edu/downloads.html)
* Hi-C data as sparse contact matrices per chromosome. Contact frequencies can be raw or normalised. Two matrix input formats are supported:
  * Either a single sparse contact matrix per chromosome, requiring three tab separated columns (position_of_bin1  position_of_bin2 contact_frequency). Genomic position reference the left most (first) basepair position of a genomic bin.
  * Or two files per chromosome in [HiC-Pro](https://github.com/nservant/HiC-Pro) output style. One three column tab separated matrix file (id_of_bin1 id_of_bin2 contact_frequency) and a coordinate file with three columns tab separated (chr pos id) linking the ids used in the matrix file to genomic positions. Genomic positions per id refer to the leftmost first basepair position of the genomic bin.
  * To split a whole genome sparse interaction matrix, subset it by chromosome name. To subset a HiC-Pro style matrix with coordinate ids and an id-bed file we wrote a perl helper script [../helper_for_preprocessing_and_analysis/extract_hicpro_matrix.pl](../helper_for_preprocessing_and_analysis/extract_hicpro_matrix.pl) run it with the `--help` flag to see usage.
    * Example to extract the chr2 intra chromosomal contact matrix:
    * `perl ./deepC/helper_for_preprocessing_and_analysis/extract_hicpro_matrix.pl --bin insitu_k652_5000_abs.bed --matrix insitu_k652_5000_iced.matrix --outmatrix insitu_k652_5000_iced.chr20.matrix --outbed insitu_k652_5000_abs.chr20.bed -chr chr2`

##### Usage

Run from command line pointing to the hic input, the chromosome sizes file and
the directory of the deepC helper scripts. Indicate the **bin.size** (resolution) in
which the Hi-C data have been analysed and that will determine the deepC resolution.
Also indicate the the **window.size** which limits the linear distance of the
deepC model to train. For all published models we used 1 Mb + 1x bin.size. Note
that the due to the vertical zigzag pole interaction encoding window.size / bin.size
must be an odd number. If you experience errors in the formatting script,
try adding 1x bin.size to you chosen window.size.

The wrapper will:
1. Extract all Hi-C interactions at a linear distance relevant within
the specified window size.
2. Convert the sparse matrix to vertical zigzag pole encoding.
3. Impute 0 frequency interactions with the mean interaction value of a 5x5 neighbourhood.
4. Transform interaction frequencies to skeleton interactions.
5. Format them to a raw text single chromosome training dataset for deepC.
6. Optional. If indicated the wrapper will also produce an example plot of the Hi-C and/or Skeleton transformed Hi-C data over a specified example region.

To construct a full training dataset simply concatinate all chromosome output
files e.g. via `cat coords_and_hic_skeleton_5kb_chr*_IMR90.bed >training_set_IMR90.txt`

*  Example command:
  ```
  Rscript ./wrapper_preprocess_hic_data.R \
    --hic.matrix=chr20_5kb.contacts.KRnorm.matrix \
    --chromosome.sizes=hg19_chrom_sizes.txt \
    --sample=IMR90 \
    --bin.size=5000 \
    --window.size=1005000 \
    --chrom=chr20 \
    --helper=../../repositories/deepC/helper_for_preprocessing_and_analysis \
    --plot.hic \
    --plot.skeleton \
    --plot.start=1e+06 \
    --plot.end=7000000 \
    --plot.height=8 \
    --plot.width=10
    ```
* Use `--help` flag for detailed parameter explanations.
