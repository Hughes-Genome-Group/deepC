<img src="docs/logo_1_transparent.png" width="150">

# deepC
A Tensorflow DL framework for predicting Hi-C chromatin interactions using megabase scale DNA sequence.

-------------------------------------------------------------------------------

### Description



### Requirements

  * python 3.5+
  * tensorflow
    * developed under 1.8.0, use legacy version
    * compatible with latest version 1.14.0, use current version
    * the legacy versions should work with every tf version 1.8.0 or higher bearing with deprecation warnings
  * additional python modules:
    * h5py
    * numpy
    * pysam
    * pybedtools and a compatible version of bedtools installed

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
