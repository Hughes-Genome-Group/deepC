#!/bin/bash

# GLOBAL RUN OPTIONS ==========================================================
SCRIPT_PATH="./deepC/tensorflow2.1plus_compatibility_version"
DATA_FILE="./minimal_training_set_example_IMR90.txt"

# Select Test and Validation Chromosomes
# Test chromosomes will be checked after each epoch
# or the specified number of chromosomes trained on
# Validation chromosomes will be hold out entirely
test_chromosomes='chr12,chr13'
validation_chromosomes='chr16,chr17'

# Settings ===================================================================
report_every='10'   # how often to report training loss every X steps
num_classes='201'  # number of classes (output vector entries), 201 for 5kb models
# 101 for 10 kb models
bp_context='1005000'  # bp context processed (1 Mb + 1x bin_size)

keep_prob_inner='0.8'  # keep probability (1-dropout rate for first conv module)
max_epoch='1'  # maximum number of epochs to train
max_chroms='18'  # maximum number of chroms to train on
save_every_chrom='3'  # how often to evaluate test performance and save checkpoint

learning_rate='0.0001'  # initial learning rate
epsilon='0.1'  # ADAM epsilon
l2_strength='0.001'
batch_size='1'  # usually restricted to 1 for memory

# specifications for first conv. module
# if using seeding / transfer learning must have the same or
# bigger dimensions then  then the architecture of the first phase model
# below are the default values for the provided first phase trained network
# the last conv. layers needs to be the same dimension as the dilation units
conv_layers='6'  # number of conv. layers in first module
hidden_units_scheme='300,600,600,900,900,100' # number of filters
kernel_width_scheme='8,8,8,4,4,1' # width of filters
max_pool_scheme='4,5,5,5,2,1'  # max pooling widths

# specifications for dilated conv. module
# dilation scheme should be chosen so that the model reaches the full sequence context
dilation_scheme='2,4,8,16,32,64,128,256,1'  # dilation rates
dilation_units='100'  # dilation units/filters throughout
dilation_width='3'
dilation_residual='True'  # if to use residual connections in the dil layers

# Transfer learning settings
seed_weights=True  # use seeding /transfer learning at all
seed_scheme='1,1,1,1,1,0'  # specify which layers to seed (1: seed, 0: not seed)
seed_file='./saved_conv_weights_human_deepc_arch.npy.npz' #trained filters phase I download from gitHub link

# Other
shuffle=True
store_dtype='bool'  # how to store the sequence
whg_fasta='./hg19.fa'  # link to whole genome fasta file for retrieving the sequences has to be indexed
use_softmasked=False  # specify if to use soft masked bases from the fasta file (lowercase). Default=False

# if multiple GPUs present select a single one to run training on
# and not block the remaining
GPU=0

train_dir='./minimal_imr90_training'

# Run ==========================================================================
python ${SCRIPT_PATH}/run_training_deepCregr.py \
        --data_file ${DATA_FILE} \
        --train_dir ${train_dir} \
        --test_chroms ${test_chromosomes} \
        --validation_chroms ${validation_chromosomes} \
        --report_every ${report_every} \
        --num_classes ${num_classes} \
        --bp_context ${bp_context} \
        --learning_rate ${learning_rate} \
        --l2_strength ${l2_strength} \
        --max_epoch ${max_epoch} \
        --max_chroms ${max_chroms} \
        --save_every_chrom ${save_every_chrom} \
        --keep_prob_inner ${keep_prob_inner} \
        --batch_size ${batch_size} \
        --conv_layers ${conv_layers} \
        --hidden_units_scheme ${hidden_units_scheme} \
        --kernel_width_scheme ${kernel_width_scheme} \
        --max_pool_scheme ${max_pool_scheme} \
        --dilation_scheme ${dilation_scheme} \
        --dilation_units ${dilation_units} \
        --dilation_width ${dilation_width} \
        --dilation_residual=${dilation_residual} \
        --epsilon ${epsilon} \
        --seed_weights=${seed_weights} \
        --seed_scheme ${seed_scheme} \
        --seed_file ${seed_file} \
        --shuffle=${shuffle} \
        --store_dtype ${store_dtype} \
        --whg_fasta ${whg_fasta} \
        --use_softmasked=${use_softmasked} \
        --gpu ${GPU}

# To continue training from a previous checkpoint use the flags:
# --model "./my_run/best_checkpoint-10000" \
# --reload_model=True
