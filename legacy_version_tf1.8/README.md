<img src="../docs/logo_1_transparent.png" width="75">

-------------------------------------------------------------------------------

### Description

Tensorflow 1.8.0 implementation of deepC. This is version of tensorflow under
which we developed deepC. This version should be compatible with all
later versions of tensorflow, bearing with deprecation warnings. Only
tensorflow 1.14.0 has been tested though.

### Content

*  **deepCregr.py** model implementation of deepC. Flexible number of convolutional
and dilated convolutional layers. Residuals (in the dilated layers) and batch normalization can be turned on.

* **deepCregr_utility.py** model implementation with more flexible intermediate outputs mainly for saliency computation

* **run_training_deepCregr.py** script for training a deepC model, requires formated/pre-rpocssed data such as the provided ones and a link to the matching reference genome.fa and .fai file

* **run_deploy_shape_deepCregr.py** script to run prediction from sequence.
Requires a trained deepC model and a bed like file with chrom start end replace
in a tab separated file with bed 0-based coordinate encoding, *replacer* being the sequence you want to exchange the respective genomic window for. Use reference if you want to run on the reference sequence.

* **run_deploy_shape_combination_deepCregr.py** same as above but applies all variants listed in the input file to the sequence before running the prediction. For multiple variants.

* **run_get_saliency.py** script for calculating the saliency with respect to input.


### Note

If help messages for command line arguments don't display, please, for now, have a look at the script. I will switch tf.app.flags for argparse in the future.
