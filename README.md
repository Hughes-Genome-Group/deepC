<img src="docs/logo_1_transparent.png" width="150">

# deepC
A Tensorflow DL framework for predicting Hi-C chromatin interactions using megabase scale DNA sequence.

-------------------------------------------------------------------------------

### Description



### Requirements

  * python 3.5 +
  * tensorflow
    * GPU support is preferable for predictions and essential for training
    * developed under 1.8.0, use legacy version
    * compatible with latest version 1.14.0, use current version
    * the legacy versions should work with every tf version 1.8.0 or higher bearing with deprecation warnings
  * additional python modules:
    * numpy (v1.16.4 or above)
    * pysam (v 0.15.2)
    * pybedtools and a compatible version of bedtools installed
  
  * R version 3.4.4 +
    * packages:
      * tidyverse (v1.2.1 or above)
      * RColorBrewer (v1.1-2 or above)
      * cowplot (v0.9.2 or above)
      * for plotting 1D tracks (e.g. DNase, ChIP-seq) rtracklayer rtracklayer (v1.38.3 or above) and dependencies are required

  * some processing helper scripts require perl (v5.26.0 or above)

### Required Ressources

  * training deepC models requires running with GPU support for several hours (up to days depending on the dataset and epochs aimed for)
  * running predictions is feasible without but runs significantly faster with a GPU
  * for example predicting the impact of a variant as in the tutorial provided requires ~5 mins with GPU support and ~ 2h on CPU.

### Installation

Clone the repository. In python script, import from directory as `import deepCregr`.

### Tutorials

Find tutorials [here](./tutorials).

### Trained Models

Download links to trained models are provided under `./models`. See the README
file there for details.

### Publication

Please refer to the bioRxiv preprint [here](https://www.biorxiv.org/content/10.1101/724005v1)

### Acknowledgements

Implementation of dilated convolutions was adapted from [wavenet](https://github.com/ibab/tensorflow-wavenet).
