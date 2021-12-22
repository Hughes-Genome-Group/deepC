<img src="../docs/logo_1_transparent.png" width="75">

# deepC - Tutorials
Tutorials for the basic deepC workflows.

-------------------------------------------------------------------------------

### Tutorials

These tutorials follow the workflow we use for deepC predictions, analysis and visualizations.

The deep learning framework and the deepC predictions itself are run via python scripts.

The downstream analysis and visualizations are done within R using tidyverse and custom functions for working with deepC predictions and Hi-C data in general.

The tutorials use the *current_version* scripts for tensorflow 1.14. If you are running on 1.8 just switch the script links to the legacy version.

Checkout the html files from your repo clone or download the html files. It has some plots that lead to a bigger file size so gitHub can't visualize it straight online. Alternatively run the R markdown files within e.g. RStudio.

1) [Predicting and plotting chromatin interactions](./tutorial_predict_and_plot.html)

2) [Formatting HiC data to skeleton data for deepC](./tutorial_format_HiC_data_for_deepC.html)

3) [Train a deepC network](./tutorial_train_a_model.md)


### Example bash scripts

1. Example bash script for [deepC training](./example_script_deepc_train.sh).

2. Example bash script for [deepC predictions](./example_script_deepc_predict.sh). Supply a bed like file chr start end seq (0-base half open genomic coordinates) see [./example_variant.bed](./example_variant.bed). A '.' is used to indicate that the genomic region is to be deleted. Alternatively, supply a DNA sequence in uppercase letters that is to replace the genomic window indicated. Lengths do not have to match and this can be a long sequence of several kb or even Mbs. For multiple modifications applied to the same sequence supply a replacement DNA sequence with all desired modifications or use [run_deploy_shape_combination_deepCregr.py](tensorflow2.1plus_compatibility_version/run_deploy_shape_combination_deepCregr.py) to apply a list of modifications.


### Wrapper scripts

In addition, we wrote wrapper scripts that perform standard preprocessing and downstream visualization tasks for deepC from the command line. These are written in R, software and data requirements are outlined below.

#### 1) Wrapper for Hi-C data pre-processing and formatting

[../helper_for_preprocessing_and_analysis/wrapper_preprocess_hic_data.R](../helper_for_preprocessing_and_analysis/wrapper_preprocess_hic_data.R)

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

Example command:

```
Rscript ./deepC/helper_for_preprocessing_and_analysis/wrapper_preprocess_hic_data.R \
  --hic.matrix=gm12878_primary_chr17_5kb.contacts.KRnorm.matrix\
  --chromosome.sizes=hg19_chrom_sizes.txt \
  --sample=GM128878 \
  --bin.size=5000 \
  --window.size=1005000 \
  --chrom=chr20 \
  --helper=./deepC/helper_for_preprocessing_and_analysis \
  --plot.hic \
  --plot.skeleton \
  --plot.start=2e+06 \
  --plot.end=2000000 \
  --plot.height=6 \
  --plot.width=8
```

Use `--help` flag for detailed parameter explanations.

#### 2) Wrapper for plotting predictions

[../helper_for_preprocessing_and_analysis/wrapper_preprocess_hic_data.R](../helper_for_preprocessing_and_analysis/wrapper_preprocess_hic_data.R)


##### Requirements

* recent version of R with the following libraries:
  * optparse
  * tidyverse
  * cowplot
  * Rcolorbrewer
  * rtracklayer - for processing 1D bigwig signals (requires GenomicRanges)
* if formatting Hi-C from sparse matrix for plotting a working version of
  perl installed and accessible for command line calls within R is required

###### Inputs

* To plot Hi-C data lot wrapper accepts Hi-C data either
  * preprocessed deepC input (use preprocess wrapper with `--no.transform`) or
  * sparse contact matrix input either only a matrix file or a matrix file and bin coordinate file (see preprocessing wrapper.) Providing the preprocessed Hi-C data speeds up the plotting script significantly.
    * if providing sparse contact matrices to plot Hi-C data a chromosome sizes file is also required
* To plot skeleton data the wrapper accepts the skeleton deepC format file produced
from the preprocessing wrapper (or otherwise).
* To plot deepC predictions the outputs of running `run_deploy_shape_deepCregr.py` or `run_deploy_shape_combination_deepCregr.py` are required.
  * A reference prediction and a variant prediction can be supplied.
  * You may use a variant prediction for reference to visualise the difference between two variants
  * If the reference prediction spans more genomic positions then the variant
  prediction, the missing variant prediction position can be filled from the reference prediction (`--fill.deepc.var`) to produce a larger variant plot. This assumes that there are no differences in the variant from the reference plot, so only use for genomic
  positions more then window.size away from any DNA alterations in the variant.
  * For examples see [example reference](./tutorials/test_predict_out/class_predictions_predict_provided_1_chr17_71000000_71999999.txt) and [example variant](./tutorials/test_variant_out/class_predictions_predict_variant_provided_1_chr17_71706322_71706671.txt) files in [./tutorials/test_predict_out](./tutorials/test_predict_out) and
[./tutorials/test_variant_out](./tutorials/test_variant_out) respectively.
* To plot bigwig tracks, e.g. DNase-seq or CTCF ChIP-seq overlapping the region of intersect
up to three bigwig files can be supplied to the wrapper.

###### Usage

Add a name tag to the plot file using `--sample` set `--out.dir` if not the current working directory. Specify the `--bin.size` and `--window.size` of the Hi-C data and deepC model.

Set he region to plot with `--chrom`, `--plot.start`, `--plot.end` making sure that the correct preprocessed and/or raw Hi-C data input for the selected chromosome is supplied.

Select which plots to produce using the flags `--plot.hic`, `--plot.skeleton`,
`--plot.deepc.ref`, `--plot.deepc.var`, `--plot.deepc.diff` & `--plot.tracks`.

Set the `--plot.width` and `--plot.height` and if desired set the relative hieghts of each sub plot with a comma separated string via `--rel.heights`.

Link to the appropriate input files: Hi-C either via `--hic.preprocessed` or via `--hic.matrix` (and `--hic.coords` for HiC-Pro style input);
Skeleton via `--skeleton.input` or is calculated from Hi-C input if desired to plot but no input available;
DeepC predictions via `--deepc.ref.input` & `--deepc.var.input`;
Bigwig tracks via `--track.input.1`, `--track.input.2`, `--track.input.3` up to three supported supply in increasing order how many desired to plot.

Set plot titles, track names and track colours as needed. And link to the helper deepC helper script directory: [../helper_for_preprocessing_and_analysis](../helper_for_preprocessing_and_analysis)

Example command:
```
Rscript wrapper_plot_deepc_predictions.R --sample=gm12878_test \
        --out.dir='.' \
        --bin.size 5000 \
        --window.size 1005000 \
        --chrom=chr17 \
        --plot.start=70500000 \
        --plot.end=72500000 \
        --plot.width=12 \
        --plot.height=16 \
        --rel.heights='0.75,0.75,0.75,0.75,1' \
        --plot.hic \
        --plot.skeleton \
        --plot.deepc.ref \
        --plot.deepc.var \
        --fill.deepc.var \
        --hic.preprocessed=coords_and_hic_skeleton_5kb_chr17_GM12878_no_transform.bed \
        --skeleton.input=coords_and_hic_skeleton_5kb_chr17_GM12878.bed \
        --deepc.ref.input=test_predict_out/class_predictions_predict_provided_1_chr17_71000000_71999999.txt \
        --deepc.var.input=test_variant_out/class_predictions_predict_variant_provided_1_chr17_71706322_71706671.txt \
        --plot.tracks \
        --track.input.1=./dnase_gm12878_encode_uw_merged_w50.bw \
        --track.input.2=./ctcf_gm12878_encode_broad_merged_w50.bw \
        --track.input.3=./dnase_gm12878_encode_uw_merged_w50.bw \
        --track.colour.1=#756bb1 \
        --helper=../../repositories/deepC/helper_for_preprocessing_and_analysis/
```

Use `--help` flag for detailed parameter explanations.
