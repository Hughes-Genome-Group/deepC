<img src="../docs/logo_1_transparent.png" width="75">

-------------------------------------------------------------------------------

### Description

Helper scripts for pre-processing Hi-C to deepC input, as well as helper for
visualization and analysis.

### Contents

* `example_head_10lines_data_K562_5kb_regression.txt` - example file including 10 lines demonstrating how to prepare Hi-C data for input to the deepC training script.
  * Format: chr start stop percentile_tags (comma separated)
  * each region is a 1 Mb + 1x bin_size region
  * the comma separated percentile tags correspond to the normalized Hi-C interaction values over the center of that region
  * regions are incremented in bin_size steps.
* `functions_for_deepC` - R functions for preprocessing Hi-C data to deepC format
and helper for analysing and visualizing deepC predictions e.g.:
  compare reference vs. variant predictions
  * virtual4C
* `functions_for_HiC.R` - general helper functions for loading, processing and visualising Hi-C data (copy from [here](https://github.com/rschwess/RonsUtilityBox/))
* [reference_manual](./reference_manual_deepC_R_utility.pdf) detailing all available functions.
* `match_query_table.pl` - helper perl script for processing Hi-C data
* `match_query_table.pl` - helper perl script for processing Hi-C data


### Notes

Some R functions require a working version of bedtools being installed and available in your standard path as queried from within R. Uses "bedtools makewindows".

Some R functions require a working version of perl being installed and available in your standard path as queried from within R.

Requires a recent version of "tidyverse" and "RColorBrewer".
