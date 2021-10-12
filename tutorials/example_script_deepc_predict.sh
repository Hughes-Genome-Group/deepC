#!/bin/bash

# # explanation
# python ../tensorflow1_version/run_deploy_shape_deepCregr.py --input example_region_short.bed \  # input deepC variant bed-like file
#   --out_dir ./test_predict_out \  #output directory
#   --name_tag predict \  # name tag to add to ouput files
#   --model ./model_deepCregr_5kb_GM12878_primary/model \  # trained deepC model downloaded and extracted
#   --genome ./hg19_chr17_fasta_for_test/hg19_chr17.fa \  # link to whole genome or chromosome wise fasta file (needs a fasta index) or test chr17 fasta file dowloaded and extracted
#   --use_softmasked=False \  #Specify if to include base pairs soft masked in the fasta file (lower case) default=False
#   --bp_context 1005000 \  # bp context (1 Mb + bin.size)
#   --add_window 500000 \  # how much bp to add to either side of the specified window
#   --num_classes 201 \  # The number of classes corresponds to the number of outputs (output bins of the vertical pole) (201 for 5kb models; 101 for 10kb models)
#   --bin_size 5000 \  # bin size matching to the model and input data selected
#   --run_on gpu  # specify to run on gpu or cpu (cpu takes significantly longer)

# actually run
python ../tensorflow1_version/run_deploy_shape_deepCregr.py --input example_region_short.bed \
  --out_dir ./test_predict_out \
  --name_tag predict \
  --model ./model_deepCregr_5kb_GM12878_primary/model \
  --genome ./hg19_chr17_fasta_for_test/hg19_chr17.fa  \
  --use_softmasked=False  \
  --bp_context 1005000 \
  --add_window 500000 \
  --num_classes 201 \
  --bin_size 5000 \
  --run_on gpu

# Run in terminal
python ../tensorflow1_version/run_deploy_shape_deepCregr.py --input example_variant.bed \
  --out_dir ./test_variant_out \
  --name_tag predict_variant \
  --model ./model_deepCregr_5kb_GM12878_primary/model \
  --genome ./hg19_chr17_fasta_for_test/hg19_chr17.fa  \
  --bp_context 1005000 \
  --add_window 500000 \
  --num_classes 201 \
  --bin_size 5000 \
  --run_on gpu
