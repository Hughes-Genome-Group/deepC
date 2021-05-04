<img src="../docs/logo_1_transparent.png" width="75">

-------------------------------------------------------------------------------

### Description

Download links to fully trained deepC models here. These are all models used for publication.

All human Hi-C data was retrieved from Rao et al. 2014 - [GSE63525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525)

Mouse Hi-C data was retrieved from Bonev et al. 2017 - [GSE96107](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96107)

Download and extract a model directory of interest. Every directory contains a
hyperparameter file listing the exact parameters the model has been trained with.
It also contains three `model.*` files that togetehr comprise a tensorflow checkpoint
saving the graph (layout) and weights of the model network. Direct the respective
deepC script to use the model as `/full_path/to_my_model/model`. Tensorflow will
recognize the model this way.

------------------------------------------

### Links

* 5 kb resolution
  * [GM12878_primary](http://datashare.molbiol.ox.ac.uk/public/rschwess/deepC/models/model_deepCregr_5kb_GM12878_primary.tar.gz)
  * [GM12878_combined](http://datashare.molbiol.ox.ac.uk/public/rschwess/deepC/models/model_deepCregr_5kb_GM12878_combined.tar.gz)
  * [K562](http://datashare.molbiol.ox.ac.uk/public/rschwess/deepC/models/model_deepCregr_5kb_K562.tar.gz)
  * mouse [mES](http://datashare.molbiol.ox.ac.uk/public/rschwess/deepC/models/model_deepCregr_5kb_mouse_ES.tar.gz)

* 10 kb resolution
  * [GM12878_primary](http://datashare.molbiol.ox.ac.uk/public/rschwess/deepC/models/model_deepCregr_10kb_GM12878_primary.tar.gz)
  * [GM12878_combined](http:/datashare.molbiol.ox.ac.uk/public/rschwess/deepC/models/model_deepCregr_10kb_GM12878_primary.tar.gz)
  * [K562](http://datashare.molbiol.ox.ac.uk/public/rschwess/deepC/models/model_deepCregr_10kb_K562.tar.gz)
  * [KBM7](http://datashare.molbiol.ox.ac.uk/public/rschwess/deepC/models/model_deepCregr_10kb_KBM7.tar.gz)
  * [HMEC](http://datashare.molbiol.ox.ac.uk/public/rschwess/deepC/models/model_deepCregr_10kb_HMEC.tar.gz)
  * [HUVEC](http://datashare.molbiol.ox.ac.uk/public/rschwess/deepC/models/model_deepCregr_10kb_HUVEC.tar.gz)
  * [IMR90](http://datashare.molbiol.ox.ac.uk/public/rschwess/deepC/models/model_deepCregr_10kb_IMR90.tar.gz)
  * [NHEK](http://datashare.molbiol.ox.ac.uk/public/rschwess/deepC/models/model_deepCregr_10kb_NHEK.tar.gz)
