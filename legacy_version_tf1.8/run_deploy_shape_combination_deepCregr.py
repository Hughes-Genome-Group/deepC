'''DESCRIPTION:
Run predictions of class changes --> shape changes over mutations
# This Version will apply all variants to the same sequence extracted
# FORMAT
chr start end replace
chromosomes must be consistent
replace columns specifies what to relace the specified window with --> can be different lengths
# NOTE 1 based indexed but bed like one over end coordinates (need to address)
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
flags.DEFINE_string('dlmodel', 'deepC', 'Specifcy the DL model file to use e.g. <endpoolDeepHaemElement>.py')
# RUN SETTINGS
flags.DEFINE_integer('batch_size', 1, 'Batch size.')
flags.DEFINE_string('out_dir', 'predictions_dir', 'Directory to store the predicted results')
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
flags.DEFINE_boolean('use_softmasked', False, 'Include soft masked sequences (lower case). If False will set them to Ns. Default = False')

# machine options
flags.DEFINE_string('run_on', 'gpu', 'Select where to run on (cpu or gpu)')
flags.DEFINE_integer('gpu', 0, 'Select a single available GPU and mask the rest. Default 0.')

# PREPARATION ------------------------------------------------------------------
# import dl model architechture selected
dlmodel = __import__(FLAGS.dlmodel)

half_bp_context = int(FLAGS.bp_context/2)

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
    hotsequence = np.zeros((len(sequence),4), dtype=dtype)

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
              seqs_placeholder: np.expand_dims(seqs[test_batch_range][0], axis=0),
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
        if line_counter % 10 == 0:
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
    variant_counter = 0
    with open(FLAGS.input, "r") as rdf:

        reference_flag = []  # some flags for process structure
        deletion_flag = []
        chroms = []
        starts = []
        ends = []
        replacers = []
        length_difference = []

        for line in rdf:

            if re.match('^#', line):  # skip comment and header lines
                continue
            variant_counter += 1
            chrom, start, end, replacer = line.split()
            start = int(start)
            end = int(end)
            chroms.append(chrom)
            starts.append(start)
            ends.append(end)
            replacers.append(replacer)

            # count bases specified
            reference_length = end - start
            # decide on mutation mode
            if re.match('reference', replacer):  # REPORT REFERENCE
                reference_flag = 1
                replacer_length = reference_length
            elif re.match('\.', replacer):  # DELETION
                deletion_flag = 1
                replacer_length = 0
            else:
                replacer_length = len(replacer) # count bases in replacer

            # store difference in bases
            length_difference.append(reference_length - replacer_length)

        # through variant files
        num_variants = len(starts)
        # sum up bp difference over variants
        total_bp_difference = sum(length_difference)

        # get first desired sequence start (start add_window and half_bp_context)
        seq_start = min(starts) - FLAGS.add_window - half_bp_context
        seq_start = customfloor(seq_start, base = FLAGS.bin_size)

        if seq_start < 0:
            to_padd = abs(seq_start)
            seq_start = 0
        else:
            to_padd = 0

        # get end of desired sequence (add total_bp_difference to end up at bin_size divsible number
        seq_end = max(ends) + FLAGS.add_window + half_bp_context
        print(seq_end)
        seq_end = customceil(seq_end, base = FLAGS.bin_size)
        print(seq_end)
        seq_end += total_bp_difference
        print(seq_end)
        # # add difference from start 1 to end last
        # spanned_length = ends[-1] - starts[0]
        # seq_end += FLAGS.bin_size - spanned_length

        # extract reference sequence -------------------------------------------
        print("Extracting %s : %s - %s" % (chroms[0], seq_start, seq_end))
        with pysam.Fastafile(FLAGS.genome) as fa:
                seq = fa.fetch(reference = chrom, start = seq_start, end = seq_end)

        # Apply Variants / Mutations
        print("Applying %s variants to sequence" % num_variants)
        s = 0
        e = starts[0] - seq_start
        var_seq = seq[0:e]
        for i in range(len(starts)):
            s = starts[i] - seq_start
            e = ends[i] - seq_start
            if replacers[i] == 'reference':
                var_seq = var_seq + seq[s:e]
            elif replacers[i] == '.':
                var_seq = var_seq  # add nothing
            else:
                var_seq = var_seq + replacers[i]

            # add next segment
            s = e
            if i == (len(starts) - 1):
                e = len(seq)
            else:
                e = starts[i+1] - seq_start
            var_seq = var_seq + seq[s:e]

        var_seq_length = len(var_seq)

        # pad (after mutating) if specified and at end of chromosome
        if to_padd > 0:
            if FLAGS.padd_ends in ['left', 'both']:
                # padd with N's
                print('padding with N\'s left wards')
                seq = 'N' * to_padd + seq
            else:
                print('%s:%s-%s is smaller then bp_context and no padding specified ... stopping' % (chrom, start, end))
                sys.exit()

        # TODO implement END padding ... need chrom sizes

        # bin the new patch ----------------------------------------------------
        i = 0
        # run_chroms = []
        run_starts = []
        run_ends = []
        run_seqs = []
        while i < (var_seq_length - FLAGS.bp_context + FLAGS.bin_size) / FLAGS.bin_size:
            js = seq_start + i * FLAGS.bin_size
            je = seq_start + i * FLAGS.bin_size + FLAGS.bp_context
            jseq = var_seq[(i*FLAGS.bin_size):((i)*FLAGS.bin_size + FLAGS.bp_context)]
            run_starts.append(js)
            run_ends.append(je)
            run_seqs.append(jseq)
            i += 1

        # Predict --------------------------------------------------------------
        # make hotcoded sequences
        hotseqs = []
        for seq in run_seqs:
            seq = get_hot_coded_seq(seq, use_soft=FLAGS.use_softmasked)  # hot encode
            hotseqs.append(seq)
        hotseqs = np.asarray(hotseqs)

        predictions = predict(
            sess,
            regression_score,
            seqs_placeholder,
            hotseqs,
            keep_prob_inner_placeholder,
            keep_prob_outer_placeholder)

        # round predictions
        np.round(predictions, 4)

        # Report -----------------------------------------------------------------------
        outfile_name = 'class_predicitions_%s.txt' % (FLAGS.name_tag)
        with open(outfile_name, "w") as fw:
            fw.write('# Combined Variants Queried')
            fw.write('\n')
            fw.write('# Mapping to relative reference coordinates: %s %s %s' % (chrom, seq_start, seq_end))
            fw.write('\n')
            fw.write('# Bp to adjust after Variant: %s' % total_bp_difference)
            fw.write('\n')
            for i in range(len(run_starts)):
                pred_out = '\t'.join(map(str, predictions[i,:]))
                fw.write("%s\t%s\t%s\t%s" % (chrom, run_starts[i], run_ends[i], pred_out))
                fw.write('\n')


# close up
