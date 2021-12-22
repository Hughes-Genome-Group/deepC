<img src="../docs/logo_1_transparent.png" width="75">

# Links to formatted example data


-------------------------------------------------------------------------------

### Links

* convolutional filter weights for transfer learning obtained from training a [deepHaem](https://github.com/rschwess/deepHaem) CNN
  * [human](http://userweb.molbiol.ox.ac.uk/datashare/rschwess/deepC/data_links/saved_conv_weights_human_deepc_arch.npy.npz) trained on 932 chromatinfeatures
  * [mouse](http://userweb.molbiol.ox.ac.uk/datashare/rschwess/deepC/data_links/saved_conv_weights_mouse_deepc_arch.npy.npz) trained on 1022 chromatin features


* example HiC [skeleton](http://userweb.molbiol.ox.ac.uk/datashare/rschwess/deepC/data_links/example_skeleton_gm12878_5kb_chr17.bed) chr17 5kb GM12878 primary

* example GM12878 [HiC sparse matrix](http://userweb.molbiol.ox.ac.uk/datashare/rschwess/deepC/data_links/gm12878_primary_chr17_5kb.contacts.KRnorm.matrix.gz) KRnorm (Rao et al.)

* formatted data Hi-C skelton data ready for deepC training
  * [GM12878 at 5kb resolution](http://userweb.molbiol.ox.ac.uk/datashare/rschwess/deepC/data_links/data_GM12878_5kb_regression.txt.tar.gz)
  * [K562 at 5kb resolution](http://userweb.molbiol.ox.ac.uk/datashare/rschwess/deepC/data_links/data_K562_5kb_regression.txt.tar.gz)
  * [IMR90 at 5kb resolution](http://userweb.molbiol.ox.ac.uk/datashare/rschwess/deepC/data_links/data_IMR90_5kb_regression.txt.gz)
  * [Mouse mES at 5kb resolution](http://userweb.molbiol.ox.ac.uk/datashare/rschwess/deepC/data_links/data_mES_5kb_regression.txt.gz)


* a [minimal training set](http://userweb.molbiol.ox.ac.uk/datashare/rschwess/deepC/data_links/minimal_training_set_example_IMR90.txt.gz) using only the first 100 lines of IMR90 5kb chromatin interactions in skeleton format. Recommended for fast testing of the training procedure.

* [hg19 chromosome sizes](http://userweb.molbiol.ox.ac.uk/datashare/rschwess/deepC/data_links/hg19_chrom_sizes.txt) listing the size in bp for each chromosome
* [hg19 chr17](http://userweb.molbiol.ox.ac.uk/datashare/rschwess/deepC/data_links/hg19_chr17_fasta_for_test.tar.gz) fasta and index for tutorial
* [hg19](http://userweb.molbiol.ox.ac.uk/datashare/rschwess/deepC/data_links/hg19_ref_genome.tar.gz) fasta and index for whole genome

* example [DNase-seq](http://userweb.molbiol.ox.ac.uk/datashare/rschwess/deepC/data_links/dnase_gm12878_encode_uw_merged_w50.bw) and [CTCF ChIP-seq](http://userweb.molbiol.ox.ac.uk/datashare/rschwess/deepC/data_links/ctcf_gm12878_encode_broad_merged_w50.bw) bigwig tracks from ENCODE for tutorial
