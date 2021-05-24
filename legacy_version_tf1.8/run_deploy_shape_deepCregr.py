'''DESCRIPTION:
Run predictions of class changes --> shape changes over mutations
# Important: all changes will be treated independendly so e.g. deletion in
mutation one will not affect the coordinates of mutation 2 and so on
# FORMAT
chr start end replace
# !!! NOTE 0 based indexed bed like format (half open) --> need to implemenet switch between 1 and 0 based
# fetching 0 based coordinates, reporting 0 based bedlike coordinates (half open)

replace columns specifies what to relace the specified window with --> can be different lengths
special cases:
    "." indicates a pure deletion
    "reference" indicates to just use the reference for this variant
'''

from __future__ import absolute_import, division, print_function
import os.path
import time
import sys
import re
import numpy as np
import tensorflow as tf
from math import log
from itertools import islice, cycle
import pysam

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
flags.DEFINE_string('dlmodel', 'deepCregr', 'Specifcy the DL model file to use e.g. <endpoolDeepHaemElement>.py')
# RUN SETTINGS
flags.DEFINE_integer('batch_size', 1, 'Batch size.')
flags.DEFINE_string('out_dir', '.', 'Directory to store the predicted results')
flags.DEFINE_string('name_tag', 'pred', 'Nametag to add to filenames')
# WHAT TO DO
flags.DEFINE_string('slize', 'all', 'Comma separated list of start and end position of columns to slice out (0) indexed. Will use all if unspecified.')
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
flags.DEFINE_integer('variant_counter', 0, 'Variant Counter variable to start from (will set +1 for first variable).')
flags.DEFINE_string('bin_steps', 'full', 'Specify if to predict in  half bin size steps of full bin size steps [hal Mainly for restarting again from a longer list. [default: 0]')
flags.DEFINE_boolean('use_softmasked', False, 'Include soft masked sequences (lower case). If False will set them to Ns. Default = False')

# machine options
flags.DEFINE_string('run_on', 'gpu', 'Select where to run on (cpu or gpu)')
flags.DEFINE_integer('gpu', 0, 'Select a single available GPU and mask the rest. Default 0.')

# PREPARATION ------------------------------------------------------------------
# import dl model architechture selected
dlmodel = __import__(FLAGS.dlmodel)

half_bp_context = int(FLAGS.bp_context/2)

half_bin_size = int(FLAGS.bin_size/2)

if FLAGS.bin_steps not in ['half', 'full']:
    print('Set bin_steps to \'half\' or \'full\'')
    sys.exit()

# prepare for column slizes if specified
if FLAGS.slize != 'all':
    slize_scheme = [x.strip() for x in FLAGS.slize.split(',')]
    slize_scheme = list(map(int, slize_scheme))

# GLOBAL OPTIONS ---------------------------------------------------------------

# HELPER FUNCTIONS -------------------------------------------------------------
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

def predict(sess,
    regression_score,
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
    predictions = np.zeros((cases, FLAGS.num_classes))  # init empty predictions array
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
        tmp_regression_score = sess.run(regression_score, feed_dict=feed_dict)
        tmp_regression_score = np.asarray(tmp_regression_score)
        tmp_regression_score = np.squeeze(tmp_regression_score)
        # add to the empty prediction scores array
        predictions[step*FLAGS.batch_size:step*FLAGS.batch_size+FLAGS.batch_size,] = tmp_regression_score
        if line_counter % 100 == 0:
            print('%s lines done ...' % line_counter)

    # handle remaining cases
    if remaining > 0:
        test_batch_range=range(cases-remaining, cases)
        # workaround for single value prediction
        if remaining == 1:
            test_batch_range=range(cases-remaining-1, cases)
        feed_dict = {
              seqs_placeholder: seqs[test_batch_range],
              keep_prob_inner_placeholder: 1.0,
              keep_prob_outer_placeholder: 1.0
              }
        tmp_regression_score = sess.run(regression_score, feed_dict=feed_dict)
        tmp_regression_score = np.asarray(tmp_regression_score)
        tmp_regression_score = np.squeeze(tmp_regression_score)
        # workaround for single value prediction (only use last remaining corresponding predicitons)
        predictions[-remaining:,] = tmp_regression_score[-remaining:]

    return predictions

''' START '''

# check if existent --> else create out_dir and Init Output File ---------------
if not os.path.exists(FLAGS.out_dir):
    os.makedirs(FLAGS.out_dir)

# Load Model -------------------------------------------------------------------
# Create a session
config = tf.ConfigProto();
if FLAGS.run_on == 'gpu':
    config.gpu_options.visible_device_list = str(FLAGS.gpu)
config.allow_soft_placement = True

# Launch Session and retrieve stored OPs and Variables
with tf.Session(config = config) as sess:
    # load meta graph and restore weights
    saver = tf.train.import_meta_graph(FLAGS.model + '.meta')
    saver.restore(sess, FLAGS.model)
    # get placeholders and ops ------------------------------------------------
    graph = tf.get_default_graph()
    seqs_placeholder = graph.get_tensor_by_name("seqs_1:0")
    labels_placeholder = graph.get_tensor_by_name("labels:0")
    keep_prob_inner_placeholder = graph.get_tensor_by_name("keep_prob_inner:0")
    keep_prob_outer_placeholder = graph.get_tensor_by_name("keep_prob_outer:0")
    regression_score = tf.get_collection("regression_score")[0]

    # read in mutation file ====================================================
    variant_counter = FLAGS.variant_counter
    with open(FLAGS.input, "r") as rdf:

        for line in rdf:

            reference_flag = 0  # some flags for process structure
            deletion_flag = 0
            if re.match('^#', line):  # skip comment and header lines
                continue
            variant_counter += 1
            print('processing entry %s' % variant_counter)
            chrom, start, end, replacer = line.split()
            start = int(start)
            end = int(end)
            # correct end for 0-based coordinates
            end = end - 1

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
            # TODO change extract mode HERE!!!!
            patch_start = customfloor(start, FLAGS.bin_size) - (FLAGS.add_window) - half_bp_context - half_bin_size
            patch_end = customceil(end, FLAGS.bin_size) + (FLAGS.add_window) + half_bp_context + half_bin_size
            patch_new_end = customceil(end, FLAGS.bin_size) + length_difference + (FLAGS.add_window) + half_bp_context + half_bin_size

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
            if deletion_flag == 1:
                seq = seq[0:(start-patch_start)] + seq[(end-patch_start+1):]
            elif reference_flag == 0:
                # partial deletion or insertion or matching replacemant
                seq = seq[:(start-patch_start)] + replacer + seq[(end-patch_start+1):]
                print("replacing with %s" % replacer)

            # print(seq)
            seq_length = len(seq)

            # bin the new patch -----------------------------------------------------
            i = 0
            # run_chroms = []
            run_starts = []
            run_ends = []
            run_seqs = []

            # select if to use half_bin_size  or not
            if FLAGS.bin_steps == 'full':
                while i < (seq_length - FLAGS.bp_context + FLAGS.bin_size) / FLAGS.bin_size:
                    js = patch_start + i * FLAGS.bin_size
                    je = patch_start + i * FLAGS.bin_size + FLAGS.bp_context
                    jseq = seq[(i*FLAGS.bin_size):((i)*FLAGS.bin_size + FLAGS.bp_context)]
                    run_starts.append(js)
                    run_ends.append(je)
                    run_seqs.append(jseq)
                    i += 1

            else:
                # half bin size step mode
                print("Using %s half in size" % half_bin_size)
                while i < (seq_length - FLAGS.bp_context + FLAGS.bin_size - half_bin_size) / half_bin_size:
                    js = patch_start + i * half_bin_size
                    je = patch_start + i * half_bin_size + FLAGS.bp_context
                    jseq = seq[(i*half_bin_size):((i)*half_bin_size + FLAGS.bp_context)]
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

            # TODO implement reporting predictions from sequence in chunks --> memory friendlier
            predictions = predict(
                sess,
                regression_score,
                seqs_placeholder,
                hotseqs,
                keep_prob_inner_placeholder,
                keep_prob_outer_placeholder)

            # round predictions
            predictions = np.round(predictions, 4)

            # Report -----------------------------------------------------------------------
            outfile_name = FLAGS.out_dir + "/" + 'class_predicitions_%s_%s_%s_%s_%s.txt' % (FLAGS.name_tag, variant_counter, chrom, start, end)
            with open(outfile_name, "w") as fw:
                fw.write('# Variant Queried: %s' % line.rstrip())
                fw.write('\n')
                fw.write('# Mapping to relative reference coordinates: %s %s %s' % (chrom, (patch_start + half_bp_context - half_bin_size), (patch_end - half_bp_context - half_bin_size)))
                fw.write('\n')
                fw.write('# Bp to adjust after Variant: %s' % length_difference)
                fw.write('\n')
                for i in range(len(run_starts)):
                    pred_out = '\t'.join(map(str, predictions[i,:]))
                    fw.write("%s\t%s\t%s\t%s" % (chrom, run_starts[i], run_ends[i], pred_out))
                    fw.write('\n')


# close up
