<img src="../docs/logo_1_transparent.png" width="75">

# Links to formatted example data


-------------------------------------------------------------------------------

### Links

* convolutional filter weights for transfer learning obtained from training a [deepHaem](https://github.com/rschwess/deepHaem) CNN
  * [human](https://datashare.molbiol.ox.ac.uk/public/project/fgenomics/rschwess/deepC/data_links/saved_conv_weights_human_deepc_arch.npy.npz) trained on 932 chromatin features
  * [mouse](https://datashare.molbiol.ox.ac.uk/public/project/fgenomics/rschwess/deepC/data_links/saved_conv_weights_mouse_deepc_arch.npy.npz) trained on 1022 chromatin features


* example HiC [skeleton](https://datashare.molbiol.ox.ac.uk/public/project/fgenomics/rschwess/deepC/data_links/example_skeleton_gm12878_5kb_chr17.bed) chr17 5kb GM12878 primary

* example GM12878 [HiC sparse matrix](https://datashare.molbiol.ox.ac.uk/public/project/fgenomics/rschwess/deepC/data_links/gm12878_primary_chr17_5kb.contacts.KRnorm.matrix.gz) KRnorm (Rao et al.)

* formatted data Hi-C skeleton data ready for deepC training
  * [GM12878 at 5kb resolution](https://datashare.molbiol.ox.ac.uk/public/project/fgenomics/rschwess/deepC/data_links/data_GM12878_5kb_regression.txt.tar.gz)
  * [K562 at 5kb resolution](https://datashare.molbiol.ox.ac.uk/public/project/fgenomics/rschwess/deepC/data_links/data_K562_5kb_regression.txt.tar.gz)
  * [IMR90 at 5kb resolution](https://datashare.molbiol.ox.ac.uk/public/project/fgenomics/rschwess/deepC/data_links/data_IMR90_5kb_regression.txt.gz)
  * [Mouse mES at 5kb resolution](https://datashare.molbiol.ox.ac.uk/public/project/fgenomics/rschwess/deepC/data_links/data_mES_5kb_regression.txt.gz)


* a [minimal training set](https://datashare.molbiol.ox.ac.uk/public/project/fgenomics/rschwess/deepC/data_links/minimal_training_set_example_IMR90.txt.gz) using only the first 100 lines of IMR90 5kb chromatin interactions in skeleton format. Recommended for fast testing of the training procedure.

* [hg19 chromosome sizes](https://datashare.molbiol.ox.ac.uk/public/project/fgenomics/rschwess/deepC/data_links/hg19_chrom_sizes.txt) listing the size in bp for each chromosome
* [hg19 chr17](https://datashare.molbiol.ox.ac.uk/public/project/fgenomics/rschwess/deepC/data_links/hg19_chr17_fasta_for_test.tar.gz) fasta and index for tutorial
* [hg19](https://datashare.molbiol.ox.ac.uk/public/project/fgenomics/rschwess/deepC/data_links/hg19_ref_genome.tar.gz) fasta and index for whole genome

* example [DNase-seq](https://datashare.molbiol.ox.ac.uk/public/project/fgenomics/rschwess/deepC/data_links/dnase_gm12878_encode_uw_merged_w50.bw) and [CTCF ChIP-seq](https://datashare.molbiol.ox.ac.uk/public/project/fgenomics/rschwess/deepC/data_links/ctcf_gm12878_encode_broad_merged_w50.bw) bigwig tracks from ENCODE for tutorial


### File Archive

All files are also archived at [Zenodo](https://zenodo.org/records/5785805).
