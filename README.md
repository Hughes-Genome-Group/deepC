<img src="docs/logo_1_transparent.png" width="150">

# deepC
A Tensorflow DL framework for predicting Hi-C chromatin interactions using megabase scale DNA sequence.

-------------------------------------------------------------------------------

### Description

This repository contains the core deepC python code, R scripts and functions for downstream analysis as well as tutorials and links to example data.

The core code is implemented in python (v3.5+) and tensorflow (v1). For downstream analysis and visualizations we use R and custom functions for handling HiC data
and deepC predictions.

### Requirements

  * python 3.5 +
  * tensorflow (tensorflow-gpu)
    * GPU support is preferable for predictions and essential for training
  * additional python modules:
    * numpy (v1.16.4 or above)
    * pysam (tested with v0.15.2)
    * pybedtools and a compatible version of bedtools installed

  * R version 3.4.4 +
    * packages:
      * tidyverse (v1.2.1 or above)
      * RColorBrewer (v1.1-2 or above)
      * cowplot (v0.9.2 or above)
      * for plotting 1D tracks (e.g. DNase, ChIP-seq) rtracklayer rtracklayer (v1.38.3 or above) and dependencies are required
    * Rstudio (not required but recommended)

  * some processing helper scripts require perl (v5.26.0 or above)

### Installation

* Make sure python 3.5-3.7 as supported by tensorflow is installed.

* Install [tensorflow](https://www.tensorflow.org/install) preferably with [GPU support](https://www.tensorflow.org/install/gpu).
  * We recommend tensorflow 2.1 but deepC was developed under v1.8 and supports (v1.8, 1.14 and 2.1 other versions have not been tested).
  * The tensorflow docker containers are the easiest way to set up tensorflow with GPU and come with the correct CUDA and cuDNN versions packaged.
  * If installing CUDA, cuDNN and tensorflow separately make sure to follow the [compatibility advice](https://www.tensorflow.org/install/source#linux)
  * To install an older version e.g. tensorflow 1 follow [this route](https://www.tensorflow.org/install/pip)

* Install additional python library (pysam and pybedtools) using e.g. pip or bioconda
  * `pip install pybedtools`
  * `pip install pysam`

* Clone the **deepC** github repository
* Check which version of tensorflow you have installed and choose the apropriate compatibility version of deepC

| tensorflow version |  CUDA version | deepC version  |
| ------------------ |:-------------:| --------------:|
| 2.1+               | 10.1          | [tensorflow2.1plus_compatibility_version](./tensorflow2.1plus_compatibility_version) |
| 2.0               | 10          | [tensorflow2.0_compatibility_version](./tensorflow2.0_compatibility_version)* |
| 1.14               | 10          | [tensorflow1_version](./tensorflow1_version) |
| 1.8               | 9          | [legacy_version_tf1.8](./legacy_version_tf1.8) |

*Compatibility with v2.0 not yet tested.

### Required Ressources

  * training deepC models requires running with GPU support for several hours (up to days depending on the dataset and epochs aimed for)
  * running predictions is feasible without but runs significantly faster with a GPU
  * for example predicting the impact of a variant as in the tutorial provided requires ~5 mins with GPU support and ~ 2h on CPU.

### Installation

Clone the repository. Make sure all dependencies are available.
To use from within a python script import as `import deepCregr`.

### Tutorials

Find tutorials [here](./tutorials).

### Trained Models

Download links to trained models are provided under `./models`. See the README
file there for details.

### Publication

Please refer to the bioRxiv preprint [here](https://www.biorxiv.org/content/10.1101/724005v1)

### Acknowledgements

Implementation of dilated convolutions was adapted from [wavenet](https://github.com/ibab/tensorflow-wavenet).
