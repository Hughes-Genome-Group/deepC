"""Trains and Evaluates deepCregr network.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os.path
import time
import sys

import numpy as np
import pysam

import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()

import deepCregr


# Basic model parameters as external flags -------------------------------------
flags = tf.app.flags
FLAGS = flags.FLAGS

#flags.DEFINE_boolean('help', False, 'For help.')

flags.DEFINE_string('data_file', '', 'Input data: pseudo bed format: chr start end and comma separated classes')

# TRAININGS SETTINGS
flags.DEFINE_string('test_chroms', 'chr12, chr13', 'Comma seperated list of test chromosoems to use ...')
flags.DEFINE_string('validation_chroms', 'chr16, chr17', 'Comma seperated list of validation chromosoems to use ...')
# flags.DEFINE_float('learning_rate_decay_steps', 5000, 'Steps to parameterize the exponential learning rate decay: LR will be LR * 0.96 every X steps.')
flags.DEFINE_integer('max_epoch', 2, 'Number of epoch through train data to run trainer.')
flags.DEFINE_integer('max_chroms', 18, 'Max number of training chromosomes to run through.')
flags.DEFINE_integer('save_every_chrom', 6, 'Save every X\'th chromosome.')
flags.DEFINE_float('keep_prob_inner', 0.8, 'Keep probability for dropout')
flags.DEFINE_float('keep_prob_outer', 0.8, 'Keep probability for dropout. LEGACY Option not used in current model implementation.')
flags.DEFINE_integer('batch_size', 1, 'Batch size.')
flags.DEFINE_float('l2_strength', 0.0001, 'L2 regularization strength.')
flags.DEFINE_boolean('shuffle', True, 'If to shuffle the trainset at the start of each epoch.')
# ARCHITECHTURE
# CONVOLUTIONAL STACK OPTIONS
flags.DEFINE_integer('conv_layers', 3, 'Number of convolutional layers.')
flags.DEFINE_string('hidden_units_scheme', '300,600,900,20', 'Comma seperated hidden units scheme. Must have length of number of conv layers specified!')
flags.DEFINE_string('kernel_width_scheme', '20,8,8,1', 'Comma seperated kernel width scheme. Must have length of number of conv layers specified!')
flags.DEFINE_string('max_pool_scheme', '5,5,5,1', 'Comma seperated max pool scheme. Must have length of number of conv layers specified!')
# DILATIONAL
flags.DEFINE_string('dilation_scheme', '2,4,8', 'Comma seperated dilation scheme to use..')
flags.DEFINE_integer('dilation_units', 20, 'Dilation Units (Filter).')
flags.DEFINE_integer('dilation_width', 3, 'Dilation Width (Only 2 supported at the moment).')
flags.DEFINE_boolean('dilation_batch_norm', False, 'If to apply batchnorm propagate residuals through the dilated layer stacks.')
# # RESIDUAL AND SKIP CONNECTIONS
flags.DEFINE_boolean('dilation_residual', False, 'If to propagate residuals through the dilated layer stacks.')
flags.DEFINE_boolean('dilation_residual_dense', False, 'If the residual/dilated layer should have a dense 1x1 convolution build in.')
# OPTIMIZER (ADAM) options
flags.DEFINE_float('learning_rate', 0.0001, 'Initial learning rate.')
flags.DEFINE_float('beta1', 0.9, 'ADAM: beta1.')
flags.DEFINE_float('beta2', 0.999, 'ADAM: beta2.')
flags.DEFINE_float('epsilon', 1e-08, 'ADAM: epsilon.')
# TRAIN LOCATION
flags.DEFINE_string('train_dir', 'training_run_data', 'Directory to put the training data.')
flags.DEFINE_string('reload_model', "False", 'If to reload a checkpoint/model file use string: \"False\" if not to reload (default); \"continue\" if to continue a training process on the same data or \"transfer\" if to load the entire model but strat training from scratch.')
flags.DEFINE_string('model', None, 'Path to checkpoint/model file.')
# Options for preseeding with pretreined weights
flags.DEFINE_boolean('seed_weights', False, 'Select if to pre seed weights with numpy array stored weights. Convolutional only ...')
flags.DEFINE_string('seed_scheme', '0,0,0', 'Specify which layers are preseeded with the weights provided. [format: 1,1,0]')
flags.DEFINE_string('seed_file', None, 'Path to saved numpy file with saved weights. Weight and bias dimensions must match with the ones specified as hyper params for this run!')
# machine options
flags.DEFINE_integer('gpu', 0, 'Select a single available GPU and mask the rest. Default 0.')
# Log Options
flags.DEFINE_integer('report_every', 100, 'Set interval of batch steps o when to report raining loss and log progress, losses and weights etc.')
# flag if to train with boolean values stored (labels and sequence)
flags.DEFINE_string('store_dtype', 'bool', 'Indicate that sequence where stored as bools rather then integers. Will convert automatically.')
# Whole genome fasta file for accessing genome sequences
flags.DEFINE_string('whg_fasta', None, 'Path to whole genome fasta file for reading the genomic sequence')
flags.DEFINE_integer('seed', 1234, 'Random seed for tensorflow (graph level).')

# GLOBAL Options ---------------------------------------------------------------
flags.DEFINE_integer('bp_context', 1000000, 'Number of classes to classify. Default 600.')
flags.DEFINE_integer('num_classes', 50, 'Number of classes to classify. Default 182.')
BP_CONTEXT = FLAGS.bp_context
NUM_CLASSES = FLAGS.num_classes

#if FLAGS.help:
#    print(FLAGS.__dict__['__flags'])

# SET RANDOM SEED --------------------------------------------------------------
np.random.seed(FLAGS.seed)  # use same seed for numpy --> for shuffeling

# Process/Prepare Train Test Valid Options ---------------------------------------
test_chromosomes = [x.strip() for x in FLAGS.test_chroms.split(',')]
test_chromosomes = list(test_chromosomes)
validation_chromosomes = [x.strip() for x in FLAGS.validation_chroms.split(',')]
validation_chromosomes = list(validation_chromosomes)

# Process/Prepare some Dilation Options ---------------------------------------
# Make list out of the passed dilation scheme string
dilation_scheme = [x.strip() for x in FLAGS.dilation_scheme.split(',')]
dilation_scheme = list(map(int, dilation_scheme))
hidden_units_scheme = [x.strip() for x in FLAGS.hidden_units_scheme.split(',')]
hidden_units_scheme = list(map(int, hidden_units_scheme))
kernel_width_scheme = [x.strip() for x in FLAGS.kernel_width_scheme.split(',')]
kernel_width_scheme = list(map(int, kernel_width_scheme))
max_pool_scheme = [x.strip() for x in FLAGS.max_pool_scheme.split(',')]
max_pool_scheme = list(map(int, max_pool_scheme))
seed_scheme = [x.strip() for x in FLAGS.seed_scheme.split(',')]
seed_scheme = list(map(int, seed_scheme))
# residual channels/units must be the same as dilation
residual_units = FLAGS.dilation_units

# TEST
if FLAGS.reload_model == "continue":
    print("Will restore previous model checkpoint for contiuning training ...")
elif FLAGS.reload_model == "transfer":
    print("Will restore previous model checkpoint for transfer learning ...")

# Assert length of schemes all match specified number of conv layers
if len(hidden_units_scheme) != FLAGS.conv_layers:
    print("Hidden Units Scheme does not have the number of entries expected from 'conv_layers' ...")
    sys.exit()
if len(kernel_width_scheme) != FLAGS.conv_layers:
    print("Hidden Width Scheme does not have the number of entries expected from 'conv_layers' ...")
    sys.exit()
if len(max_pool_scheme) != FLAGS.conv_layers:
    print("Max Pool Scheme does not have the number of entries expected from 'conv_layers' ...")
    sys.exit()
if FLAGS.seed_weights and len(seed_scheme) != FLAGS.conv_layers:
    print("Seed Scheme does not have the number of entries expected from 'conv_layers' ...")
    sys.exit()

# STARTING ---------------------------------------------------------------------
print("STARTING: ...")

# HELPER FUNCTIONS -------------------------------------------------------------
def placeholder_inputs(batch_size, dtype):
  """Generate placeholder variables to represent the input tensors.

  These placeholders are used as inputs by the rest of the model building
  code and will be fed from the downloaded data in the .run() loop, below.

  Args:
    batch_size: The batch size will be baked into both placeholders.
    dtype: dtaype in which seq and labels are/ will be stored

  Returns:
    seqs_placeholder: Sequences (hot coded) placeholder.
    labels_placeholder: Labels placeholder.
  """
  # sess = tf.InteractiveSession()
  if dtype == 'bool':
      seqs_placeholder = tf.compat.v1.placeholder(tf.bool, [None, BP_CONTEXT, 4], name='seqs')
      labels_placeholder = tf.compat.v1.placeholder(tf.uint8, shape=[None, NUM_CLASSES], name='labels')
  if dtype == 'uint8':
      seqs_placeholder = tf.compat.v1.placeholder(tf.uint8, [None, BP_CONTEXT, 4], name='seqs')
      labels_placeholder = tf.compat.v1.placeholder(tf.uint8, shape=[None, NUM_CLASSES], name='labels')
  else:
      seqs_placeholder = tf.compat.v1.placeholder(tf.int32, [None, BP_CONTEXT, 4], name='seqs')
      labels_placeholder = tf.compat.v1.placeholder(tf.int32, shape=[None, NUM_CLASSES], name='labels')
  # Note that the shapes of the placeholders match the shapes
  return seqs_placeholder, labels_placeholder
def load_chromosomes(genome_file):
  """ Load genome segments from either a FASTA file or
          chromosome length table. -- from Basenji: Kelly et al 2018"""

  # is genome_file FASTA or (chrom,start,end) table?
  file_fasta = (open(genome_file).readline()[0] == '>')

  chrom_segments = {}

  if file_fasta:
    fasta_open = pysam.Fastafile(genome_file)
    for i in range(len(fasta_open.references)):
      chrom_segments[fasta_open.references[i]] = [(0, fasta_open.lengths[i])]
    fasta_open.close()

  else:
    for line in open(genome_file):
      a = line.split()
      chrom_segments[a[0]] = [(0, int(a[1]))]

  return(chrom_segments)

def get_hot_coded_seq(sequence, dtype="int32"):
    """Convert a 4 base letter sequence to 4-row x-cols hot coded sequence"""
    # initialise empty
    hotsequence = np.zeros((len(sequence),4), dtype=dtype)
    # set hot code 1 according to gathered sequence
    for i in range(len(sequence)):
        if sequence[i] == 'A':
            hotsequence[i,0] = 1
        elif sequence[i] == 'C':
            hotsequence[i,1] = 1
        elif sequence[i] == 'G':
            hotsequence[i,2] = 1
        elif sequence[i] == 'T':
            hotsequence[i,3] = 1
    # return the numpy array
    return hotsequence
def get_hot_coded_seq_bool(sequence):
    """Convert a 4 base letter sequence to 4-row x-cols hot coded sequence"""
    # initialise empty
    hotsequence = np.zeros((len(sequence),4), dtype='bool')
    # set hot code 1 according to gathered sequence
    for i in range(len(sequence)):
        if sequence[i] == 'A':
            hotsequence[i,0] = True
        elif sequence[i] == 'C':
            hotsequence[i,1] = True
        elif sequence[i] == 'G':
            hotsequence[i,2] = True
        elif sequence[i] == 'T':
            hotsequence[i,3] = True
    # return the numpy array
    return hotsequence

def split_bins_regr(label, num_classes):
  '''Helper Function to create a binarized numpy array for representing
  the respective class associations'''
  # init regr bin representatons -----------------------------------------------
  # Go through labels per seq/position and fill regr bin values
  label_bin = np.zeros((len(label), num_classes), dtype="uint32")
  for j in range(len(label)):
      l = label[j].split(",")  # split by comma
      for i in range(len(l)):
          label_bin[j,i] = l[i]

  return(label_bin)
def read_train_file(fi, num_classes, test_chromosomes, validation_chromosomes, dtype='int32'):
    '''helper method: read taining file: store chrom start end - dict,
    for each chrom store first and last coordinate of each chromosome to extract
    sequence later on'''
    with open(fi, "r") as f:
        chroms = []
        starts = []
        ends = []
        label = []
        chrom_dict = {}  # init chromosome dict
        for i,l in enumerate(f):
            l = l.rstrip()
            l = l.split("\t")
            if l[0] == "chr":
                continue
            chrom = l[0]
            chroms.append(chrom)
            start = int(l[1])
            starts.append(start)
            end = int(l[2])
            ends.append(end)
            # store labels and all labels temporaryly for setting up binary labels later
            label.append(l[3])
            if chrom in chrom_dict:
                if start < chrom_dict[chrom]['start']:
                    chrom_dict[chrom]['start'] = start
                if end > chrom_dict[chrom]['end']:
                    chrom_dict[chrom]['end'] = end
            # check if same chromosome
            else:
                chrom_dict[chrom] = {}
                chrom_dict[chrom]['start'] = start
                chrom_dict[chrom]['end'] = end

    # convert chrom start and stops to numpy arrays
    chroms = np.array(chroms)  # chromomes
    # and [start, end] x input_lines matrix
    position = np.array((starts, ends), dtype='i')
    position = np.transpose(position)

    # create representation
    regr_bin = split_bins_regr(label, num_classes)

    # split all up into test validaiton and training as specified
    # chrom dict
    test_chrom_d = {}
    valid_chrom_d = {}
    train_chrom_d = {}
    for c in chrom_dict.keys():
        if c in test_chromosomes:
            test_chrom_d[c] = chrom_dict[c]
        elif c in validation_chromosomes:
            valid_chrom_d[c] = chrom_dict[c]
        else:
            train_chrom_d[c] = chrom_dict[c]

    # get line indeces for splitting the rest
    train_index = []
    test_index = []
    valid_index = []
    for i in range(len(chroms)):
        if chroms[i] in test_chromosomes:
            test_index.append(i)
        elif chroms[i] in validation_chromosomes:
            valid_index.append(i)
        else:
            train_index.append(i)
    # chroms
    test_chroms = chroms[test_index]
    valid_chroms = chroms[valid_index]
    train_chroms = chroms[train_index]
    # positions
    test_position = position[test_index,:]
    valid_position = position[valid_index,:]
    train_position = position[train_index,:]
    # regr_bin
    test_regr_bin = regr_bin[test_index,:]
    valid_regr_bin = regr_bin[valid_index,:]
    train_regr_bin = regr_bin[train_index,:]

    return(train_chrom_d,
      train_chroms,
      train_position,
      train_regr_bin,
      test_chrom_d,
      test_chroms,
      test_position,
      test_regr_bin,
      valid_chrom_d,
      valid_chroms,
      valid_position,
      valid_regr_bin)
def get_chr_seq(whg_fasta, chromosome, dtype = "uint8"):
    '''Get entire chromosome sequence as hot encoded sequence to query from'''
    with pysam.Fastafile(whg_fasta) as fa:
        seq = fa.fetch(reference = chromosome)
        # convert to one hot encoding
        if dtype == 'bool':
            seq = get_hot_coded_seq_bool(seq)
        elif dtype == 'uint8':
            seq = get_hot_coded_seq(seq, dtype)
        else:
            seq = get_hot_coded_seq(seq, dtype)
    return(seq)
def fill_seqs(chr_seq, positions_stored, indeces, batch_size, seq_length, dtype = "uint8"):
    # init
    if dtype == 'bool':
        feed_seqs = np.zeros((batch_size, seq_length, 4), dtype = 'bool')
    elif dtype == 'uint8':
        feed_seqs = np.zeros((batch_size, seq_length, 4), dtype = 'uint8')
    else:
        feed_seqs = np.zeros((batch_size, seq_length, 4), dtype = 'int')
    # fill from temp chromosome sequence stored
    for i in range(batch_size):
        feed_seqs[i,] = chr_seq[positions_stored[indeces[i],0]:positions_stored[indeces[i],1]]

    return(feed_seqs)

def make_numpy_seed(
    seed_weights_list_loaded,
    seed_scheme,
    hidden_units_scheme,
    kernel_width_scheme):
    '''Make a list of numpy arrays to preseed the conv_max pool part of the net.
    If the new dimensions (hidden units, kernel_width, ...) match with the new run
    this will pre seed them as they are provided in the numpy file.
    If the dimenstions
    are larger, this will sample the excess weights to use from the distribution of
    weights to preseed in the respective same layer.
    If the new dimensions are
    smaller then the previous one (throw an error).
    Arguments:
        seed_weights_list_loaded: numpt list of arrays in right dimension.
        seed_scheme: list of seed scheme 1s and 0s
    Returns:
        list of numpy arrays (arr_0 ... arr_9) depending on the layers to preseed,
        where 0,2,4,... are the 3D weigths and 1,3,5,... are the 1D biases'''
    # get how many layers to preseed
    seed_sum = sum(seed_scheme)
    # 0,2,4, ... is weights /// 1,3,5 ... is biases
    # initialise empty arrays according to shapes
    seed_weights_list = {}
    excess_layer_count = 0
    for i in range(seed_sum):  #TODO get number of layers to preload from preloading scheme
        w = i * 2
        b = w + 1
        weights_load_string = 'arr_' + str(w)
        biases_load_string = 'arr_' + str(b)
        # get respective dimensions
        hu = hidden_units_scheme[i]
        kw = kernel_width_scheme[i]
        if i == 0:
            indim = 4
        else:
            indim = hidden_units_scheme[i-1]
        init_weights = np.zeros((kw, indim, hu), dtype='float32')# initialise weights
        init_biases = np.zeros((hu), dtype='float32')  # and biases
        # get respective preloadable seeds
        seed_weights = seed_weights_list_loaded[weights_load_string]
        seed_biases = seed_weights_list_loaded[biases_load_string]
        seed_weights_shape = seed_weights.shape
        # get excess dimension range if present
        excess_neurons_in = indim - seed_weights_shape[1]
        excess_neurons_out = hu - seed_weights_shape[2]
        excess_neurons_width = kw - seed_weights_shape[0]

        # check for negative dimensions abort
        if excess_neurons_in < 0 and excess_neurons_out and excess_neurons_width < 0:
            print('Specified Convolutional Dimensions are smaller then what you are trying to preseed with. Only equal or larger allowed ...')
            sys.exit()

        if excess_neurons_width > 0 or excess_neurons_in > 0 or excess_neurons_out > 0:
            print('New specified dimensions in layer %s are larger then seed provided. Sampling excess weights from that distribution per layer!' % (i+1))
            excess_layer_count += 1
        # seed with whats available
        init_weights[0:seed_weights_shape[0], 0:seed_weights_shape[1], 0:seed_weights_shape[2]] = seed_weights

        # if excess present sample excess weights from seeded (maintain distribution)
        if excess_neurons_out > 0:  # A) excess in hidden units dimension
            to_sample = kw * indim * excess_neurons_out
            sampled = np.random.choice(seed_weights.flatten(), size = to_sample, replace = True)
            # fill array
            init_weights[0:kw, 0:indim, seed_weights_shape[2]:hu] = np.reshape(
                sampled, (kw, indim, excess_neurons_out))

        if excess_neurons_in > 0:  # B) excess in hidden units previous dimension (or input)
            to_sample = kw * excess_neurons_in * seed_weights_shape[2]
            sampled = np.random.choice(seed_weights.flatten(), size = to_sample, replace = True)
            # fill array
            init_weights[0:kw, seed_weights_shape[1]:indim, 0:seed_weights_shape[2]] = np.reshape(
                sampled, (kw, excess_neurons_in, seed_weights_shape[2]))

        if excess_neurons_width > 0:  # B) excess in hidden units previous dimension (or input)
            to_sample = excess_neurons_width * seed_weights_shape[1] * seed_weights_shape[2]
            sampled = np.random.choice(seed_weights.flatten(), size = to_sample, replace = True)
            # fill array
            init_weights[seed_weights_shape[0]:kw, 0:seed_weights_shape[1], 0:seed_weights_shape[2]] = np.reshape(
                sampled, (excess_neurons_width, seed_weights_shape[1], seed_weights_shape[2]))

        # fill with seeds
        init_biases[0:seed_weights_shape[2]] = seed_biases
        if excess_neurons_out > 0:
            to_sample = excess_neurons_out
            sampled = np.random.choice(seed_biases.flatten(), size = to_sample, replace = True)
            init_biases[seed_weights_shape[2]:hu] = sampled

        seed_weights_list[weights_load_string] = init_weights
        seed_weights_list[biases_load_string] = init_biases

    return(seed_weights_list, excess_layer_count)

def do_eval(sess,
            # eval_correct,
            eval_loss,
            seqs_placeholder,
            labels_placeholder,
            chrom_d_test,
            chroms_test,
            positions_test,
            labels_test,
            keep_prob_inner_placeholder,
            keep_prob_outer_placeholder
            ):
    """Runs one evaluation against the full epoch of test data.
    Return test accuracy and mean test loss per batch.

    Args:
    sess: The session in which the model has been trained.
    # eval_correct: The Tensor that returns the number of correct predictions.
    seqs_placeholder: The sequences placeholder.
    labels_placeholder: The labels placeholder.
    keep_prob_pl: placeholder for the keep probability
    lines: Opend lines object of training data file
    cases: number of lines/cases in input file
    """
    # And run one epoch of eval.
    true_count = 0  # Counts the number of correct predictions.
    test_loss = 0
    cases = labels_test.shape[0]
    test_step = 0

    # pick chromosomes in training set
    test_chroms = list(chrom_d_test.keys())
    # set up test index
    test_cases = np.shape(labels_test)[0]
    test_index = np.asarray(range(test_cases))

    # RUN THROUGH TRAIN CHROMOSOMES
    for current_chr in test_chroms:

        # pick sub training index
        sub_test_index = test_index[chroms_test == current_chr,]
        sub_test_cases = len(sub_test_index)
        # preload the respective chromosome sequence
        tmp_chr_seq = get_chr_seq(FLAGS.whg_fasta, current_chr, FLAGS.store_dtype)
        sub_step = 0
        test_loss = 0

        while (sub_step + 1) * FLAGS.batch_size < sub_test_cases:
            sub_step += 1
            test_step += 1
            # Fill a feed dictionary with the respective set of seqs and labels
            test_batch_start = sub_step * FLAGS.batch_size
            test_batch_end = sub_step * FLAGS.batch_size+FLAGS.batch_size
            test_batch_range = range(test_batch_start, test_batch_end)
            feed_dict = {
                  seqs_placeholder: fill_seqs(tmp_chr_seq, positions_test, test_batch_range, FLAGS.batch_size, BP_CONTEXT, dtype=FLAGS.store_dtype),
                  labels_placeholder: labels_test[test_batch_range,],
                  keep_prob_inner_placeholder: 1.0,
                  keep_prob_outer_placeholder: 1.0
                  }

            tmp_test_loss = sess.run([eval_loss], feed_dict=feed_dict)
            # print(tmp_test_loss)

            # true_count += tmp_true_count
            test_loss = test_loss + tmp_test_loss[0]

    test_loss = test_loss / test_step

    print('Num examples: %d Test Loss: %0.04f' %
    (cases, test_loss))

    return test_loss  # return precision

def run_training():
  """Train for a number of steps."""

  # Get the sets labels and chromosome positions
  # split into train test and validation -------------------------------

  train_chrom_d, train_chroms, train_position, train_regr_bin, test_chrom_d, test_chroms, test_position, test_regr_bin, valid_chrom_d, valid_chroms, valid_position, valid_regr_bin = read_train_file(FLAGS.data_file, NUM_CLASSES, test_chromosomes, validation_chromosomes, FLAGS.store_dtype)

  print(test_chroms)
  print(valid_chroms)

  # # test a input
  # print('Test chromosomes')
  # print(train_chroms)
  # print('Test positions:')
  # print(train_position)
  # print('Train classes_bin:')
  # print(train_regr_bin)
  # sys.exit()

  print('training_data')
  print(np.shape(train_regr_bin))
  training_cases = np.shape(train_regr_bin)[0]

  print('test_data')
  print(test_regr_bin.shape)

  print('validation_data')
  print(valid_regr_bin.shape)

  # make indices for shuffeling the data
  training_index = np.asarray(range(training_cases))

  # load seed weights if specified
  if FLAGS.seed_weights:
    print("Loading saved weights ...")
    seed_weights_list = np.load(FLAGS.seed_file)
    # get numpy seed array augment excess dimension with sampling where needed
    seed_weights_list, excess_count = make_numpy_seed(seed_weights_list, seed_scheme, hidden_units_scheme, kernel_width_scheme)
  else:
    seed_weights_list = ""

  # Write configuration/parameters to specified train dir
  # WRITE HYPER PARAMETERS and DATA PROPERTIES TO RUN LOG FILE
  current_time = time.localtime()
  timestamp = str(current_time[0]) + str(current_time[1]) + str(current_time[2]) + str(current_time[3]) + str(current_time[4])

  if not os.path.exists(FLAGS.train_dir):  # make train dir
      os.makedirs(FLAGS.train_dir)

  param_log_file = open(FLAGS.train_dir + '/hyperparameters_' + timestamp + '.log', 'w')
  param_log_file.write("# Hyperparameters for dilationDHS run at: " + str(current_time))
  param_log_file.write("\n\nInput: " + str(BP_CONTEXT) + " bp")
  param_log_file.write("\nTest Chromosomes: " + str(FLAGS.test_chroms))
  param_log_file.write("\nValidation Chromosomes: " + str(FLAGS.validation_chroms))
  param_log_file.write("\n\nArchitechture:")
  for i in range(FLAGS.conv_layers):
    j = i + 1
    param_log_file.write("\nHidden Units Layer %s: %s" % (j, hidden_units_scheme[i]))
    param_log_file.write("\nKernel Width Layer %s: %s" % (j, kernel_width_scheme[i]))
    param_log_file.write("\nMax Pool Layer %s: %s" % (j, max_pool_scheme[i]))
  param_log_file.write("\nDilatio Scheme: " + str(FLAGS.dilation_scheme))
  param_log_file.write("\nDialtion Units: " + str(FLAGS.dilation_units))
  param_log_file.write("\nDilation Width: " + str(FLAGS.dilation_width))
  param_log_file.write("\nDilation Batch Norm: " + str(FLAGS.dilation_batch_norm))
  param_log_file.write("\nDilation with Residuals: " + str(FLAGS.dilation_residual))
  param_log_file.write("\nDilation with dense Residuals: " + str(FLAGS.dilation_residual_dense))
  if FLAGS.seed_weights:
      param_log_file.write("\nPre-seeding with saved weights: " + FLAGS.seed_file)
      if excess_count > 0:
          param_log_file.write("\nSome new dimensions are larger then the pre seeded ones;\n sampling the excess weights from same distribution per layer.")
  param_log_file.write("\n\nTraining Parameters:")
  param_log_file.write("\nBatch size: " + str(FLAGS.batch_size))
  param_log_file.write("\nDropout Keep Probability Inner: " + str(FLAGS.keep_prob_inner))
  param_log_file.write("\nDropout Keep Probability Outer: " + str(FLAGS.keep_prob_outer) + " # LEGACY option not used in final models published.")
  param_log_file.write("\nLearning Rate (Intital): " + str(FLAGS.learning_rate))
  param_log_file.write("\n(ADAM) Beta 1: " + str(FLAGS.beta1))
  param_log_file.write("\n(ADAM) Beta 2: " + str(FLAGS.beta2))
  param_log_file.write("\n(ADAM) Epsilon: " + str(FLAGS.epsilon))
  param_log_file.write("\nMaximum Epoch Number: " + str(FLAGS.max_epoch) + "\n")
  param_log_file.write("\nL2 Regularizer Strength: " + str(FLAGS.l2_strength) + "\n")
  param_log_file.write("\nShuffle Training Set for each CHromosome Call: " + str(FLAGS.shuffle) + "\n\n")

  param_log_file.write("____________________________________________________\n")
  param_log_file.close()

  # Tell TensorFlow that the model will be built into the default Graph.
  with tf.Graph().as_default():

    # set seed
    tf.compat.v1.set_random_seed(FLAGS.seed)

    # Generate placeholders for the seqs and labels (and dropout prob).
    seqs_placeholder, labels_placeholder = placeholder_inputs(FLAGS.batch_size, FLAGS.store_dtype)
    keep_prob_inner_placeholder = tf.compat.v1.placeholder(tf.float32, name='keep_prob_inner')
    keep_prob_outer_placeholder = tf.compat.v1.placeholder(tf.float32, name='keep_prob_outer')

    # Building the Graph -------------------------------------------------------
    # Create a variable to track the global step.
    global_step = tf.Variable(0, name='global_step', trainable=False)

    # Ops to calc regression_score
    regression_score = deepCregr.inference(
        seqs_placeholder,
        FLAGS.conv_layers,
        hidden_units_scheme,
        kernel_width_scheme,
        max_pool_scheme,
        dilation_scheme,
        FLAGS.dilation_units,
        FLAGS.dilation_width,
        FLAGS.dilation_residual,
        FLAGS.dilation_residual_dense,
        FLAGS.dilation_batch_norm,
        NUM_CLASSES,
        FLAGS.batch_size,
        keep_prob_inner_placeholder,
        keep_prob_outer_placeholder,
        FLAGS.seed_weights,
        seed_scheme,
        seed_weights_list
        )
    tf.compat.v1.add_to_collection("regression_score", regression_score)

    # Add to the Graph the Ops for loss calculation.
    loss = deepCregr.loss(regression_score, labels_placeholder, FLAGS.l2_strength, FLAGS.batch_size)
    loss_test = deepCregr.loss_test(regression_score, labels_placeholder, FLAGS.batch_size)
    tf.compat.v1.add_to_collection("loss_test", loss_test)

    # Add to the Graph the Ops that calculate and apply gradients.
    train_op = deepCregr.training(
        loss,
        FLAGS.learning_rate,
        FLAGS.beta1,
        FLAGS.beta2,
        FLAGS.epsilon,
        global_step)
    tf.compat.v1.add_to_collection("train_op", train_op)

    # # Add the Ops to compare the regression_score to the labels during evaluation.
    # eval_op = deepCregr.evaluation(regression_score, labels_placeholder)
    # tf.compat.v1.add_to_collection("eval_op", eval_op)

    # Build the summary Tensor based on the TF collection of Summaries.
    summary = tf.compat.v1.summary.merge_all()

    # Create a saver for writing training checkpoints.
    saver = tf.compat.v1.train.Saver(max_to_keep=5)

    # init op
    init = tf.compat.v1.global_variables_initializer()

    # Create a session for running Ops on the Graph.
    config = tf.compat.v1.ConfigProto();
    config.gpu_options.visible_device_list = str(FLAGS.gpu)
    config.allow_soft_placement = True
    sess = tf.compat.v1.Session(config = config)

    # Instantiate a SummaryWriter to output summaries and the Graph.
    summary_writer = tf.compat.v1.summary.FileWriter(FLAGS.train_dir, graph=sess.graph)

    # print(tf.get_collection(tf.GraphKeys.GLOBAL_VARIABLES))
    # sys.exit()

    # reload model if specified or initialize ----------------------------------
    if FLAGS.reload_model == "continue":
        print("Restoring previous model checkpoint for contiuning training ...")
        saver.restore(sess, FLAGS.model)
        # reload global step
        total_step = tf.train.global_step(sess, global_step)
    elif FLAGS.reload_model == "transfer":
        print("Restoring previous model checkpoint for transfer learning ...")
        saver.restore(sess, FLAGS.model)
        # reset global step
        total_step = 0
    else:
        sess.run(init)
        total_step = 0

    # Start the TRAINING loop ==================================================
    epoch = 0
    step = 0
    best_loss = 5000000.0

    # chromosome step counter
    chrom_counter = 0

    # pick chromosomes in training set
    to_train_chroms = list(train_chrom_d.keys())
    to_train_chroms = sorted(to_train_chroms)
    print(to_train_chroms)

    while epoch < FLAGS.max_epoch:

        # RUN THROUGH TRAIN CHROMOSOMES
        for current_chr in to_train_chroms:

            print("Training on Chromosome %s" % current_chr)

            sub_step = 0
            chrom_counter += 1

            # pick sub training index
            sub_training_index = training_index[train_chroms == current_chr,]
            sub_training_cases = len(sub_training_index)
            # preload the respective chromosome sequence
            tmp_chr_seq = get_chr_seq(FLAGS.whg_fasta, current_chr, FLAGS.store_dtype)

            # :q! training_index if flagged
            if FLAGS.shuffle:
                np.random.shuffle(sub_training_index)

            while ((sub_step + 1) * FLAGS.batch_size) < sub_training_cases:
                sub_step += 1
                step += 1
                total_step += 1
                start_time = time.time()

                # Fill a feed dictionary with the respective set of seqs and labels
                batch_start = (sub_step)*FLAGS.batch_size
                batch_end = (sub_step)*FLAGS.batch_size+FLAGS.batch_size
                batch_range=range(batch_start, batch_end)
                # convert to shuffled traiing index
                batch_range = sub_training_index[batch_range]
                batch_range = np.sort(batch_range)
                batch_range = batch_range.tolist()
                # create feed dictionaries
                feed_dict = {
                      seqs_placeholder: fill_seqs(tmp_chr_seq, train_position,
                        batch_range, FLAGS.batch_size, BP_CONTEXT, dtype = FLAGS.store_dtype),
                      labels_placeholder: train_regr_bin[batch_range,],
                      keep_prob_inner_placeholder: FLAGS.keep_prob_inner,
                      keep_prob_outer_placeholder: FLAGS.keep_prob_outer
                      }

                # Run one step of the model.  The return values are the activations
                # from the `train_op` (which is discarded) and the `loss`
                _, loss_value, regr = sess.run([train_op, loss, regression_score], feed_dict=feed_dict)

                # # TEST PRINTS
                # print('labels')
                # print(train_regr_bin[batch_range,])
                # print('regression score')
                # print(regr)
                # print('loss value')
                # print(loss_value)

                # Write the summaries and print and overview every X steps =============
                if step % FLAGS.report_every == 0:
                    duration = time.time() - start_time # get step time
                    # Print status to stdout.
                    print('Step %d: loss = %.2f (%.3f sec)' % (step, loss_value, duration))
                    # Update the events file.
                    summary_str = sess.run(summary, feed_dict=feed_dict)
                    summary_writer.add_summary(summary_str, total_step)
                    summary_writer.flush()

            # Every 6 Training Chromsomes evaluate test loss on Test Chromosomes and potentially save check point =============
            if chrom_counter % FLAGS.save_every_chrom == 0:
                print('Trained on %s Chromosomes' % chrom_counter)
                print('Test Data Accuracy Eval:')
                test_loss = do_eval(sess,
                        loss_test,
                        seqs_placeholder,
                        labels_placeholder,
                        test_chrom_d,
                        test_chroms,
                        test_position,
                        test_regr_bin,
                        keep_prob_inner_placeholder,
                        keep_prob_outer_placeholder
                        )
                test_loss_summary = tf.Summary(value=[tf.Summary.Value(tag="test_loss", simple_value=test_loss)])
                summary_writer.add_summary(test_loss_summary, total_step)
                # Save Checkpoint if test loss is smaller then the previous best ===
                if best_loss > test_loss:
                    checkpoint_file = os.path.join(FLAGS.train_dir, 'best_checkpoint')
                    saver.save(sess, checkpoint_file, global_step=total_step)
                    best_loss = test_loss

            if chrom_counter >= FLAGS.max_chroms:
                # end training
                print('Reached max chroms finishing up ...')
                sys.exit()


        # EPOCH ends after running through all training chromosomes ones =======
        epoch += 1  # count up epoch
        step = 0  # reset step
        print("Epoch %s done:" % epoch)

def main(_):
  run_training()

if __name__ == '__main__':
  tf.compat.v1.app.run()
