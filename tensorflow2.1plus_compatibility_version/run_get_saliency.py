'''DESCRIPTION:
Run Get Saliency Scores:  Gradient of Output with respect to input and/or intermediate layers.
# FORMAT
chr start end replacer --> will work like variant deploy in terms of applying variants
# NOTE 0 based indexed bed like format --> need to implemenet switch between 1 and 0 based
# fetching 0 based coordinates, reporting 0 based bedlike coordinates
'''

from __future__ import absolute_import, division, print_function
import os.path
import time
import sys
import re
import numpy as np
from math import log
from itertools import islice, cycle
import pysam

import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()

# helper custom rounding to arbitrary base
def customround(x, base=5):
    return int(base * round(float(x)/base))
def customfloor(x, base=5):
    f = int(base * round(float(x)/base))
    if f > x:
        f = f - base
    return f
def customceil(x, base=5):
    f = int(base * round(float(x)/base))
    if f < x:
        f = f + base
    return f

# Basic model parameters as external flags -------------------------------------
flags = tf.app.flags
FLAGS = flags.FLAGS
flags.DEFINE_string('dlmodel', 'deepCregr_utility', 'Specifcy the DL model file to use e.g. <endpoolDeepHaemElement>.py')
# RUN SETTINGS
flags.DEFINE_integer('batch_size', 1, 'Batch size.')
flags.DEFINE_string('out_dir', '.', 'Directory to store the predicted results')
flags.DEFINE_string('name_tag', 'sal', 'Nametag to add to filenames')
flags.DEFINE_boolean('mutate_sequence', False, 'Seelect if the sequence specified should be mutated. If False will just extract the sequence specified as a full bed window. [Default  False]')
flags.DEFINE_string('report_format', 'long', 'How to report the saliency scores [\'long\' or \'broad\']')
flags.DEFINE_string('report_summary', 'abs_sum', 'What to report in long format: \'sum\' of all columns (bases); \'abs_sum\' the absolute sum; or all 4 \'columns\' one per base (A,C,G,T)')
flags.DEFINE_boolean('times_seq', True, 'Multiply gradients/saliency with hot encoded sequence element wise.')
flags.DEFINE_boolean('report_wig', False, '[True False] select if to also write a wiggle track. Currently only works for report_summary abs_sum or sum. If selected, bp_context must amtch the input sequence window (bed format).')
flags.DEFINE_string('wig_name_tag', 'sal', 'Nametag for wiggle file if specified to report this way.')
flags.DEFINE_string('wig_description', 'saliency per bp', 'Description for wiggle track if specified o report.')
flags.DEFINE_boolean('filter_zeros', True, '[True False] select if to filter out value zero entries (onlf for txt output).')
# Network Architechture
flags.DEFINE_boolean('report_conv_hidden_state', False, '[True False] Select if only to report the hidden state after convolution.')
flags.DEFINE_integer('conv_layers', 6, 'Number of convolutional layers.')
flags.DEFINE_string('hidden_units_scheme', '300,600,600,900,900,100', 'Comma seperated hidden units scheme. Must have length of number of conv layers specified!')
flags.DEFINE_string('kernel_width_scheme', '8,8,8,4,4,1', 'Comma seperated kernel width scheme. Must have length of number of conv layers specified!')
flags.DEFINE_string('max_pool_scheme', '4,5,5,5,2,1', 'Comma seperated max pool scheme. Must have length of number of conv layers specified!')

flags.DEFINE_string('dilation_scheme', '2,4,8,16,32,64,128,256,1', 'Comma seperated dilation scheme to use..')
flags.DEFINE_integer('dilation_units', 100, 'Dilation Units (Filter).')
flags.DEFINE_integer('dilation_width', 3, 'Dilation Width (Only 2 supported at the moment).')
flags.DEFINE_boolean('dilation_batch_norm', False, 'If to apply batchnorm propagate residuals through the dilated layer stacks.')
# # RESIDUAL AND SKIP CONNECTIONS
flags.DEFINE_boolean('dilation_residual', True, 'If to propagate residuals through the dilated layer stacks.')
flags.DEFINE_boolean('dilation_residual_dense', False, 'If the residual/dilated layer should have a dense 1x1 convolution build in.')
# Seeding
flags.DEFINE_boolean('seed_weights', False, 'Select if to pre seed weights with numpy array stored weights. Convolutional only ...')
flags.DEFINE_string('seed_scheme', '0,0,0,0,0,0', 'Specify which layers are preseeded with the weights provided. [format: 1,1,0]')
flags.DEFINE_string('seed_file', None, 'Path to saved numpy file with saved weights. Weight and bias dimensions must match with the ones specified as hyper params for this run!')
# EXTERNAL files
flags.DEFINE_string('input', '', 'Must be a variant file specifying the mutations to apply to the reference (custom made format for now)!')
flags.DEFINE_string('model', './model', 'Checkpoint of model file to be tested. (Full path to model without suffix!)')
flags.DEFINE_string('genome', 'hg19.fasta', 'Full path to fasta reference genome of interest to extract the sequence from.')
flags.DEFINE_string('padd_ends', 'none', 'Specify if to padd with half times bp_context N\'s to predict over chromosome ends [left, right, none, both].')
# Data Options
flags.DEFINE_integer('bp_context', 1010000, 'Basepairs per feature.')
flags.DEFINE_integer('add_window', 0, 'Basepairs to add around variant of interest for prediction and hence visualization later.')
flags.DEFINE_integer('num_classes', 101, 'Number of classes.')
flags.DEFINE_integer('bin_size', 10000, 'Bin size to apply when running over the new sequence.')
flags.DEFINE_string('store_dtype', 'float32', 'Indicate that sequence where stored as bools rather then integers. Will convert automatically.')
flags.DEFINE_boolean('use_softmasked', False, 'Include soft masked sequences (lower case). If False will set them to Ns. Default = False')

# machine options
flags.DEFINE_string('run_on', 'gpu', 'Select where to run on (cpu or gpu)')
flags.DEFINE_integer('gpu', 0, 'Select a single available GPU and mask the rest. Default 0.')

# PREPARATION ------------------------------------------------------------------
# import dl model architechture selected
dlmodel = __import__(FLAGS.dlmodel)

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

half_bp_context = FLAGS.bp_context/2

# GLOBAL OPTIONS ---------------------------------------------------------------

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
        seqs_placeholder = tf.compat.v1.placeholder(tf.bool, [None, FLAGS.bp_context, 4], name='seqs')
        labels_placeholder = tf.compat.v1.placeholder(tf.uint8, shape=[None, FLAGS.num_classes], name='labels')
    if dtype == 'uint8':
        seqs_placeholder = tf.compat.v1.placeholder(tf.uint8, [None, FLAGS.bp_context, 4], name='seqs')
        labels_placeholder = tf.compat.v1.placeholder(tf.uint8, shape=[None, FLAGS.num_classes], name='labels')
    if dtype == 'float32':
        seqs_placeholder = tf.compat.v1.placeholder(tf.float32, [None, FLAGS.bp_context, 4], name='seqs')
        labels_placeholder = tf.compat.v1.placeholder(tf.float32, shape=[None, FLAGS.num_classes], name='labels')
    else:
        seqs_placeholder = tf.compat.v1.placeholder(tf.int32, [None, FLAGS.bp_context, 4], name='seqs')
        labels_placeholder = tf.compat.v1.placeholder(tf.int32, shape=[None, FLAGS.num_classes], name='labels')
    return seqs_placeholder, labels_placeholder

# Helper get hotcoded sequence
def get_hot_coded_seq(sequence, use_soft=False):
    """Convert a 4 base letter sequence to 4-row x-cols hot coded sequence"""
    # initialise empty
    hotsequence = np.zeros((len(sequence),4), dtype = 'uint8')

    # transform to uppercase if using softmasked sequences
    if use_soft:
        sequence = sequence.upper()

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


def map_saliency(sess,
    saliency_op,
    seqs_placeholder,
    seqs,
    keep_prob_inner_placeholder,
    keep_prob_outer_placeholder
    ):
    """Make predictions --> get sigmoid output of net per sequence and class"""
    cases = seqs.shape[0]
    line_counter = 0
    batches_to_run = cases // FLAGS.batch_size
    # cover cases where remainder cases are left
    remaining = cases - FLAGS.batch_size * batches_to_run
    saliencies = np.zeros((cases, FLAGS.bp_context, 4))  # init empty predictions array

    for step in range(batches_to_run):
        line_counter += 1
        test_batch_start = step * FLAGS.batch_size
        test_batch_end = step * FLAGS.batch_size + FLAGS.batch_size
        test_batch_range=range(test_batch_start, test_batch_end)
        feed_dict = {
              seqs_placeholder: seqs[test_batch_range],
              keep_prob_inner_placeholder: 1.0,
              keep_prob_outer_placeholder: 1.0
              }
        # print(seqs[test_batch_range])
        # print(seqs[test_batch_range].shape)
        tmp_saliencies_score = sess.run(saliency_op, feed_dict=feed_dict)
        tmp_saliencies_score = np.asarray(tmp_saliencies_score)

        # print("Shape out of session ...")
        # print(tmp_saliencies_score.shape)
        #
        # print("Shape after Squeeze ...")
        # print(tmp_saliencies_score.shape)
        tmp_saliencies_score = np.squeeze(tmp_saliencies_score)

        # add to the empty saliency scores array
        saliencies[step*FLAGS.batch_size:step*FLAGS.batch_size+FLAGS.batch_size,] = tmp_saliencies_score

        if line_counter % 10 == 0:
            print('%s lines done ...' % line_counter)

    # # handle remaining cases
    # if remaining > 0:
    #     test_batch_range=range(cases-remaining, cases)
    #     # workaround for single value prediction
    #     if remaining == 1:
    #         test_batch_range=range(cases-remaining-1, cases)
    #     feed_dict = {
    #           seqs_placeholder: seqs[test_batch_range],
    #           keep_prob_inner_placeholder: 1.0,
    #           keep_prob_outer_placeholder: 1.0
    #           }
    #     tmp_saliencies_score = sess.run(input_gradients, feed_dict=feed_dict)
    #     tmp_saliencies_score = np.asarray(tmp_saliencies_score)
    #     tmp_saliencies_score = np.squeeze(tmp_saliencies_score)
    #     # workaround for single value prediction (only use last remaining corresponding predicitons)
    #     saliencies[-remaining:,] = tmp_saliencies_score[-remaining:]

    return saliencies

''' START '''

# check if existent --> else create out_dir and Init Output File ---------------
if not os.path.exists(FLAGS.out_dir):
    os.makedirs(FLAGS.out_dir)

# Tell TensorFlow that the model will be built into the default Graph.
with tf.Graph().as_default():

    # Generate placeholders for the seqs and labels (and dropout prob).
    seqs_placeholder, labels_placeholder = placeholder_inputs(FLAGS.batch_size, FLAGS.store_dtype)
    keep_prob_inner_placeholder = tf.compat.v1.placeholder(tf.float32, name='keep_prob_inner')
    keep_prob_outer_placeholder = tf.compat.v1.placeholder(tf.float32, name='keep_prob_outer')

    # Building the Graph -------------------------------------------------------

    # Ops to calc regression_score
    regression_score = dlmodel.inference(
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
        FLAGS.num_classes,
        FLAGS.batch_size,
        keep_prob_inner_placeholder,
        keep_prob_outer_placeholder,
        FLAGS.seed_weights,
        seed_scheme,
        "",
        FLAGS.report_conv_hidden_state
        )

    # define saliency
    saliency_op = tf.gradients(regression_score[:,], seqs_placeholder)
    print("Saliency OP:")
    print(saliency_op[0])

    # SET SAVER ---------------------------------------------------
    saver = tf.train.Saver()

    # init op
    init = tf.compat.v1.global_variables_initializer()

    # Load Model ---------------------------------------------------------------
    # Create a session
    config = tf.compat.v1.ConfigProto();
    if FLAGS.run_on == 'gpu':
        config.gpu_options.visible_device_list = str(FLAGS.gpu)
    config.allow_soft_placement = True

    # Launch Session and retrieve stored OPs and Variables
    with tf.compat.v1.Session(config = config) as sess:

        # init all variables
        sess.run(init)

        #  restore weights
        print("Restoring subset of weights ....... ")
        saver.restore(sess, FLAGS.model)

        # read in region file ====================================================
        region_counter = 0
        with open(FLAGS.input, "r") as rdf:

            for line in rdf:

                reference_flag = 0  # some flags for process structure
                deletion_flag = 0
                if re.match('^#', line):  # skip comment and header lines
                    continue
                region_counter += 1
                print('processing entry %s' % region_counter)
                chrom, start, end, replacer = line.split()
                start = int(start)
                end = int(end)

                if FLAGS.mutate_sequence == True:
                    '''Take and mutate sequence specified and extract enough sequence surrounding'''

                    # count bases specified
                    reference_length = end - start + 1
                    # decide on mutation mode
                    if re.match('reference', replacer):  # REPORT REFERENCE
                        reference_flag = 1
                        replacer_length = reference_length
                    elif re.match('\.', replacer):  # DELETION
                        deletion_flag = 1
                        replacer_length = 0
                    else:
                        replacer_length = len(replacer) # count bases in replacer

                    # get difference in bases ----------------------------------------------
                    length_difference = reference_length - replacer_length

                    # set new coordinates --------------------------------------------------
                    new_start = start
                    new_end = end - length_difference
                    if deletion_flag == 1:  # adjust for full deletions as the first base pair goes as well
                        new_start = new_start - 1
                        new_end = new_end - 1

                    # save coordinates that are plotted/analysed over respeective to the reference
                    relative_reference_start = start - FLAGS.add_window
                    relative_reference_start = end + FLAGS.add_window

                    # round new coordinates to full bins and set sequence window to extract -----
                    # add 990,000 bp to either side of the last bin start / endpool
                    patch_start = customfloor(start, FLAGS.bin_size) - (FLAGS.add_window) - half_bp_context
                    patch_end = customceil(end, FLAGS.bin_size) + (FLAGS.add_window) + half_bp_context
                    patch_new_end = customceil(end, FLAGS.bin_size) + length_difference + (FLAGS.add_window) + half_bp_context

                else:
                    '''Take start and stop from coordinates as specified : full bed coordinates'''
                    # correct bed coordinates
                    patch_start = start
                    patch_new_end = end

                # CONTINUE COMBINED ============================================
                # check start and end of range
                # set start_diff if sequence to query is over the chromosome ends --> ready to padd
                start_diff = 0
                seq_start = patch_start
                if patch_start < 0:
                    start_diff = abs(patch_start)
                    seq_start = 0 # cover over the border cases for sequence retrival

                # extract reference sequence -------------------------------------------
                with pysam.Fastafile(FLAGS.genome) as fa:
                        seq = fa.fetch(reference = chrom, start = seq_start, end = patch_new_end)

                # pad if specifiedand at end of chromosome
                if start_diff > 0:
                    if FLAGS.padd_ends in ['left', 'both']:
                        # padd with N's
                        print('padding with N\'s left wards')
                        seq = 'N' * start_diff + seq
                    else:
                        print('%s:%s-%s is smaller then bp_context and no padding specified ... skipping' % (chrom, start, end))
                        continue
                else:
                # END padding ... need chrom sizes
                    end_diff = patch_new_end - (len(seq) + patch_start)
                    if end_diff > 0:
                        if FLAGS.padd_ends in ['right', 'both']:
                            print('padding with N\'s right wards')
                            seq = seq + 'N' * end_diff
                        else:
                            print('%s:%s-%s is smaller then bp_context and no padding specified ... skipping' % (chrom, start, end))
                            continue

                # mutate reference sequence --------------------------------------------
                if FLAGS.mutate_sequence:
                    if deletion_flag == 1:
                        seq = seq[0:(start-patch_start)] + seq[(end-patch_start+1):]
                    elif reference_flag == 0:
                        # partial deletion or insertion or matching replacemant
                        seq = seq[:(start-patch_start)] + replacer + seq[(end-patch_start+1):]

                seq_length = len(seq)

                # bin the new patch -----------------------------------------------------
                i = 0
                # run_chroms = []
                run_starts = []
                run_ends = []
                run_seqs = []

                # select if to use half_bin_size  or not
                while i < (seq_length - FLAGS.bp_context + FLAGS.bin_size) / FLAGS.bin_size:
                    js = patch_start + i * FLAGS.bin_size
                    je = patch_start + i * FLAGS.bin_size + FLAGS.bp_context
                    jseq = seq[(i*FLAGS.bin_size):((i)*FLAGS.bin_size + FLAGS.bp_context)]
                    run_starts.append(js)
                    run_ends.append(je)
                    run_seqs.append(jseq)
                    i += 1

                # Predict ----------------------------------------------------------------
                # make hotcoded sequences
                hotseqs = []
                for seq in run_seqs:
                    seq = get_hot_coded_seq(seq, use_soft=FLAGS.use_softmasked)  # hot encode
                    hotseqs.append(seq)
                hotseqs = np.asarray(hotseqs)

                saliencies = map_saliency(
                    sess,
                    saliency_op,
                    seqs_placeholder,
                    hotseqs,
                    keep_prob_inner_placeholder,
                    keep_prob_outer_placeholder)

                # ROUND
                saliencies = np.round(saliencies, 4)

                # print("Run saliency mapping ...")
                # print(saliencies.shape)

                if FLAGS.times_seq == True:
                    # print("Performing element wise multiplication with hot coded sequence ...")
                    saliencies = np.multiply(saliencies, hotseqs)

                # Report -----------------------------------------------------------------------
                outfile_name = FLAGS.out_dir + "/" + 'saliencies_%s_%s_%s_%s_%s.txt' % (FLAGS.name_tag, region_counter, chrom, start, end)
                with open(outfile_name, "w") as fw:
                    fw.write('# Region Queried: %s' % line.rstrip())
                    fw.write('\n')
                    if FLAGS.report_format == 'broad':
                        # report in broad format
                        for i in range(len(run_starts)):
                            # get predictions per base
                            sal_A_out = '\t'.join(map(str, saliencies[i,:,0]))
                            sal_C_out = '\t'.join(map(str, saliencies[i,:,1]))
                            sal_G_out = '\t'.join(map(str, saliencies[i,:,2]))
                            sal_T_out = '\t'.join(map(str, saliencies[i,:,3]))

                            fw.write("%s\t%s\t%s\tA\t%s\n" % (chrom, run_starts[i], run_ends[i], sal_A_out))
                            fw.write("%s\t%s\t%s\tC\t%s\n" % (chrom, run_starts[i], run_ends[i], sal_C_out))
                            fw.write("%s\t%s\t%s\tG\t%s\n" % (chrom, run_starts[i], run_ends[i], sal_G_out))
                            fw.write("%s\t%s\t%s\tT\t%s\n" % (chrom, run_starts[i], run_ends[i], sal_T_out))
                    else:
                        # report in long format
                        for i in range(len(run_starts)):
                            sal_A_out = saliencies[i,:,0]
                            sal_C_out = saliencies[i,:,1]
                            sal_G_out = saliencies[i,:,2]
                            sal_T_out = saliencies[i,:,3]

                            j = 0
                            if FLAGS.report_summary == 'columns':
                                for p in range(run_starts[i], run_ends[i]):
                                    pos = p + 1
                                    fw.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, pos, sal_A_out[j], sal_C_out[j], sal_G_out[j], sal_T_out[j]))
                                    j = j + 1
                            elif FLAGS.report_summary == 'sum':
                                for p in range(run_starts[i], run_ends[i]):
                                    pos = p + 1
                                    sal_summary = sal_A_out[j] + sal_C_out[j] + sal_G_out[j] + sal_T_out[j]
                                    j = j + 1
                                    if FLAGS.filter_zeros == True and sal_summary == 0:
                                            continue
                                    fw.write("%s\t%s\t%s\t\n" % (chrom, pos, sal_summary))

                            elif FLAGS.report_summary == 'abs_sum':
                                for p in range(run_starts[i], run_ends[i]):
                                    pos = p + 1
                                    sal_summary = abs(sal_A_out[j]) + abs(sal_C_out[j]) + abs(sal_G_out[j]) + abs(sal_T_out[j])
                                    j = j + 1
                                    if FLAGS.filter_zeros == True and sal_summary == 0:
                                            continue
                                    fw.write("%s\t%s\t%s\t\n" % (chrom, pos, sal_summary))

                            else:
                                print("No valid report mode specfied for long format! Select: [\'abs_sum\', \'sum\', \'columns\']")

                    # Report Wig track if specfieid -----------------------------------------------------------------------
                    print(FLAGS.report_wig)
                    if FLAGS.report_wig == True:
                        if FLAGS.report_summary in ['abs_sum', 'sum']:
                            outwig_name = FLAGS.out_dir + "/" + 'saliencies_%s_%s_%s_%s_%s.wig' % (FLAGS.name_tag, region_counter, chrom, start, end)
                            wiggle_name = 'saliency_%s_%s_%s' % (FLAGS.report_summary, FLAGS.wig_name_tag, region_counter)
                            wiggle_description = '%s %s %s' % (FLAGS.wig_description, FLAGS.report_summary, region_counter)
                            with open(outwig_name, "w") as wig:
                                wig.write('track type=wiggle_0 name="%s" description="%s" visibility=full\n' % (wiggle_name, wiggle_description))
                                wig.write('fixedStep chrom=%s start=%s step=1 span=1\n' % (chrom, start + 1))

                                for i in range(len(run_starts)):
                                    sal_A_out = saliencies[i,:,0]
                                    sal_C_out = saliencies[i,:,1]
                                    sal_G_out = saliencies[i,:,2]
                                    sal_T_out = saliencies[i,:,3]

                                    j = 0
                                    if FLAGS.report_summary == 'sum':
                                        for p in range(run_starts[i], run_ends[i]):
                                            sal_summary = sal_A_out[j] + sal_C_out[j] + sal_G_out[j] + sal_T_out[j]
                                            wig.write("%s\n" % (sal_summary))
                                            j = j + 1
                                    elif FLAGS.report_summary == 'abs_sum':
                                        for p in range(run_starts[i], run_ends[i]):
                                            sal_summary = abs(sal_A_out[j]) + abs(sal_C_out[j]) + abs(sal_G_out[j]) + abs(sal_T_out[j])
                                            wig.write("%s\n" % (sal_summary))
                                            j = j + 1
                                    else:
                                        print("No valid report modefor wig reporting specfied! Select: [\'abs_sum\', \'sum\']")
