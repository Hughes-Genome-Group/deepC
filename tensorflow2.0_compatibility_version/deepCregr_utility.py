"""Tensorflow implementation of DeepC using convolutional and dilated layers.
------------------------------------------------------------------

Note:
Operations especially the implementation of the causal_convolution in "2D" but
only diluting along 1 dimension is taken and adapted from gitHUB tensorflow-wavenet.

UTITLITY Version:
Rebuild the inference with more flexible outs.
1) Saliency mapping: > remove all tf.cast functions for salieny mapping (gradient of regression score with respect to input)
    need to feed sequence as float value!

Acknowledgement:
    Code for performing dilated convolutions has been adapted from
    https://github.com/ibab/tensorflow-wavenet

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import math
import re
import numpy as np

import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()

# OPERATIONS DEFINITION ============================================
def time_to_batch(value, dilation, name=None):
    with tf.name_scope('time_to_batch'):
        shape = tf.shape(value)
        pad_elements = dilation - 1 - (shape[1] + dilation - 1) % dilation
        padded = tf.pad(value, [[0, 0], [0, pad_elements], [0, 0]])
        reshaped = tf.reshape(padded, [-1, dilation, shape[2]])
        transposed = tf.transpose(reshaped, perm=[1, 0, 2])
        return tf.reshape(transposed, [shape[0] * dilation, -1, shape[2]])

def batch_to_time(value, dilation, name=None):
    with tf.name_scope('batch_to_time'):
        shape = tf.shape(value)
        prepared = tf.reshape(value, [dilation, -1, shape[2]])
        transposed = tf.transpose(prepared, perm=[1, 0, 2])
        return tf.reshape(transposed,
                          [tf.math.divide(shape[0], dilation), -1, shape[2]])

def dilated_conv(value, filter_, dilation, name='dilated_conv'):
    with tf.name_scope(name):
        filter_width = tf.shape(filter_)[0]
        if dilation > 1:
            transformed = time_to_batch(value, dilation)
            conv = tf.nn.conv1d(transformed, filter_, stride=1,
                                padding='SAME')
            restored = batch_to_time(conv, dilation)
        else:
            restored = tf.nn.conv1d(value, filter_, stride=1, padding='SAME')
        # Remove excess elements at the ends.
        # desired width
        out_width = tf.shape(value)[1]
        # get difference
        diff = out_width - tf.shape(restored)[1]
        # get half difference to remove SAME padding zeros
        diff = diff//2
        result = tf.slice(restored,
                          [0, 0, 0],
                          [-1, out_width, -1])

        return result

def _activation_summary(x):
  """Helper to create summaries for activations.
  Creates a summary that provides a histogram of activations.
  Creates a summary that measures the sparsity of activations.
  Args:
    x: Tensor
  Returns:
    nothing
  """
  # session. This helps the clarity of presentation on tensorboard.
  tensor_name = x.op.name
  tf.summary.histogram(tensor_name + '/activations', x)
  tf.summary.scalar(tensor_name + '/sparsity', tf.nn.zero_fraction(x))

def preload_variable(name, data):
    '''Create variable from numpy data.'''
    variable = tf.Variable(data, name=name)
    return variable

def create_variable(name, shape):
    '''Create a convolution filter variable with the specified name and shape,
    and initialize it using Xavier initialition.'''
    initializer = tf.keras.initializers.glorot_normal()
    variable = tf.Variable(initializer(shape=shape), name=name)
    # variable = tf.Variable(tf.truncated_normal(shape), name=name)
    return variable

def create_bias_variable(name, shape):
    '''Create a bias variable with the specified name and shape and initialize
    it to zero.'''
    initializer = tf.constant_initializer(value=0, dtype=tf.float32)
    variable = tf.Variable(initializer(shape=shape), name=name)
    return  variable

def slice_seqs(positions, remainder, bp_context, batch_size, chrom_seq, translation_tensor):
    with tf.name_scope('Get_Chrom_Seq'):

        # adjust positions to match 3 base coding -->

        # slice sequence from pre stored whole chrom sequence on gpu memory
        indices = (tf.range(bp_context) + positions[:,tf.newaxis])[...,tf.newaxis]
        print(indices.shape)
        seqs = tf.gather_nd(chrom_seq, indices)
        print('sequence shape after fetching from chrom')
        print(seqs.get_shape().as_list())

        seqs = tf.cast(seqs, 'int32')

        # TRANSLATE BACK FROM 3 base chunk uin8 ENCODING
        seqs = seqs[...,tf.newaxis]
        seqs = tf.gather_nd(translation_tensor, seqs)

        print('sequence shape fetching hot encoding')
        print(seqs.get_shape().as_list())

        seqs = tf.transpose(seqs, perm=[0,1,3,2])
        seqs = tf.reshape(seqs, [batch_size,(bp_context*3),5])

        # strip Ns / strip additional bases extracted from 3 base basis (remainder)
        seqs = seqs[:,:,1:6]

        print('sequence shape after translation')
        print(seqs.get_shape().as_list())
        # return sequences as batch size
        return(seqs)

# Define Dilated Convolutional Layer
def dilated_layer(name,
                  input_batch,
                  dilation,
                  dilation_width,
                  dilation_units,
                  residual=False,
                  dense_residual=False,
                  to_batch_norm=False):
    '''Create a dilation layer:
        INPUT:
        name: must be unique for graph purpose
        input_batch: 3D input tensor batch, length, channels/width
        Current implementation keeps the channels/hidden units intakt
        dialton: dialtion rate to apply
        dilation_width: with of the dilation filter (only 2 supported?)
        dilation_units: dilation_units or channels
        residual: True/False --> select if to propagate residual in that layer/stack
        to_batch_norm: True/False select of to perform batch norm at every layer
        RETURNS:
        3D tensor batch, length-dilation rate, channels width'''
    with tf.name_scope(name):
        # get shapes
        channels = input_batch.get_shape().as_list()[2]
        # create variables
        dilation_weights = create_variable('dilation_weights', [dilation_width, channels, dilation_units])
        dilation_biases = create_bias_variable('dilation_biases', [dilation_units])
        gate_weights = create_variable('gate_weights', [dilation_width, channels, dilation_units])
        gate_biases = create_bias_variable('gate_biases', [dilation_units])
        # redisual and skip
        if residual == True:
            if dense_residual == True:
                dense_weights = create_variable('dense_weights', [dilation_units, channels, dilation_units])
                dense_biases = create_bias_variable('dense_biases', [dilation_units])
                # skip_weights = create_variable('skip_weights', [1, dilation_units, skip_units])
                # skip_biases = create_bias_variable('skip_biases', [skip_units])

        # define convolutional steps
        dilated = tf.add(dilated_conv(input_batch, dilation_weights, dilation=dilation), dilation_biases)
        gated = tf.add(dilated_conv(input_batch, gate_weights, dilation=dilation), gate_biases)
        dilated_gated = tf.tanh(dilated) * tf.sigmoid(gated)

        if residual == True:
            # if dense residual connection desired make a 1x1 convolutiion before adding
            if dense_residual == True:
                # 1x1 dense convolution for residual
                transformed = tf.nn.conv1d(
                    dilated_gated, dense_weights, stride=1, padding="SAME", name="dense")
                transformed = transformed + dense_biases
                # add up residual to 1x1 transformed output
                out = input_batch + transformed
            else:
                # else just add input_batch shortcut to dilated/gated
                out = input_batch + dilated_gated
        else:
            # else dilated gated is out
            out = dilated_gated
        # # The 1x1 conv to produce the skip output
        # skip_cut = tf.shape(out)[1] - output_width
        # out_skip = tf.slice(out, [0, skip_cut, 0], [-1, -1, -1])
        # weights_skip = variables['skip']
        # skip_contribution = tf.nn.conv1d(
        #     out_skip, weights_skip, stride=1, padding="SAME", name="skip")
        # # batch norm
        if to_batch_norm == True:
            out = tf.layers.batch_normalization(out)

        # make summary histograms of weights
        tf.summary.histogram(name + '_dilation_weights', dilation_weights)
        tf.summary.histogram(name + '_dilation_biases', dilation_biases)
        tf.summary.histogram(name + '_gate_weights', gate_weights)
        tf.summary.histogram(name + '_gate_biases', gate_biases)
        if dense_residual == True:
            tf.summary.histogram(name + '_dense_weights', dense_weights)
            tf.summary.histogram(name + '_dense_biases', dense_biases)
        # tf.summary.histogram(name + '_skip_weights', skip_weights)
        # tf.summary.histogram(name + '_skip_biases', skip_biases)

        return out

# Define Convolutional Layer
def convolutional_layer(name,
                  input_batch,
                  units,
                  kernel_width,
                  pool_width,
                  keep_prob,
                  to_seed,
                  seed_weights,
                  seed_biases,
                  to_batch_norm=False):
    '''Create a convolutional layer:
        INPUT:
        name: must be unique for graph purpose
        input_batch: 3D input tensor batch, length, channels/width
        units: hidden units (# kernels)
        kernel_width: width of the convolutional kernels/filters
        pool_width: (max) pool width
        keep_prob: dropout keep probability
        to_seed: True / False if to pre seed weights and biases in this layer
        seed_weights: numpy array of seed weights
        seed_biases: numpy array of seed biases
        to_batch_norm: True/False select of to perform batch norm at every layer
        RETURNS:
        3D tensor batch, length/pool_width, channels width'''
    with tf.name_scope(name):
        # get shapes
        channels = input_batch.get_shape().as_list()[2]
        # create variables
        if to_seed:
            weights = preload_variable("weights", seed_weights)
            biases = preload_variable("biases", seed_biases)
        else:
            weights = create_variable('weights', [kernel_width, channels, units])
            biases = create_bias_variable('biases', [units])
        # define convolutional steps
        conv = tf.add(tf.nn.conv1d(input_batch, weights, stride=1, padding='SAME'), biases)
        conv = tf.nn.relu(conv)
        # make summary histograms of weights
        tf.summary.histogram(name + '_conv_weights', weights)
        tf.summary.histogram(name + '_conv_biases', biases)
        # activation summary
        _activation_summary(conv)
        # Max Pool
        conv = tf.layers.max_pooling1d(conv, pool_width, strides=pool_width, padding='same', name=str(name+'max_pool'))
        # Dropout
        out = tf.nn.dropout(conv, rate = 1 - keep_prob)
        # # batch norm
        if to_batch_norm == True:
            out = tf.layers.batch_normalization(out)
        return out

# INFERENCE ===================================================================
def inference(seqs,
              conv_layers,
              hidden_units_scheme,
              kernel_width_scheme,
              max_pool_scheme,
              dilation_scheme,
              dilation_units,
              dilation_width,
              dilation_residual,
              dilation_residual_dense,
              dilation_batch_norm,
              num_classes,
              batch_size,
              keep_prob_inner,
              keep_prob_outer,
              seed_weights,
              seed_scheme,
              seed_weights_list,
              report_conv_hidden_state
              ):
    """INFERENCE
    Args:

    Returns:
        regression score or hidden representation after convolutional but before dilated convolutional layer
    """

    print('seqs shape')
    print(seqs.get_shape().as_list())

    current_layer = seqs

    # Convolutional Stack with Max Pooling =====================================
    # run an inital dilated layer with dilation 1 to map to the dilational unit output
    with tf.name_scope('Convolutional_stack'):
        for i in range(conv_layers):
            j = i + 1
            k = i * 2
            if seed_weights and seed_scheme[i] == 1:
                weights_load_string = 'arr_' + str(k)
                biases_load_string = 'arr_' + str(k+1)
                print('Pre-seeding Layer: ' + str(j))
                current_layer = convolutional_layer(
                    'conv_layer{}'.format(j),
                    current_layer,
                    hidden_units_scheme[i],
                    kernel_width_scheme[i],
                    max_pool_scheme[i],
                    keep_prob_inner,
                    True,
                    seed_weights_list[weights_load_string],
                    seed_weights_list[biases_load_string],
                    to_batch_norm=False)
            else:
                current_layer = convolutional_layer(
                    'conv_layer{}'.format(j),
                    current_layer,
                    hidden_units_scheme[i],
                    kernel_width_scheme[i],
                    max_pool_scheme[i],
                    keep_prob_inner,
                    False,
                    "dummy",
                    "dummy",
                    to_batch_norm=False)
            print('Conv %s shape' % j)
            print(current_layer.get_shape().as_list())

    # Report only HIdden STate afte r Conv layers (optional) ===================
    if report_conv_hidden_state:
        ''''only report hte hidden state after the convolutional stacks'''
        hidden_conv_state = current_layer
        return(hidden_conv_state)

    # Dilational Layers stack ==================================================
    # run an inital dilated layer with dilation 1 to map to the dilational unit output
    with tf.name_scope('dilated_stack'):
        current_layer = dilated_layer(
            'dilated_layer1',
            current_layer,
            1,
            dilation_width,
            dilation_units,
            residual = dilation_residual,
            dense_residual = dilation_residual_dense,
            to_batch_norm = dilation_batch_norm)
        print('Dilated shape')
        print(current_layer.get_shape().as_list())
        for i, dilation in enumerate(dilation_scheme):
            i = i+1  # skipping 0 count as this is pre-established
            current_layer = dilated_layer(
                'dilated_layer{}'.format(i),
                current_layer,
                dilation,
                dilation_width,
                dilation_units,
                residual = dilation_residual,
                dense_residual = dilation_residual_dense,
                to_batch_norm=dilation_batch_norm)

            print('Dilated shape')
            print(current_layer.get_shape().as_list())

    # reshape for FC layer
    with tf.name_scope('reshape_layer'):
        fully_connected_width = current_layer.get_shape().as_list()[1] * dilation_units
        current_layer = tf.reshape(current_layer, [batch_size, fully_connected_width])
        print('fully connection reshaped')
        print(current_layer.get_shape().as_list())

    # Final full connection(s) into logits
    with tf.name_scope('final_dense'):
        weights = create_variable('weights', [fully_connected_width, num_classes])
        biases = create_bias_variable('biases', [num_classes])
        regression_score = tf.add(tf.matmul(current_layer, weights), biases)
        print('Regression score shape')
        print(regression_score.get_shape().as_list())
        _activation_summary(regression_score)
        tf.summary.histogram('final_dense_weights', weights)

    # return regression_score, into_dilation
    return regression_score

def loss(regression_score, labels, l2_regularization_strength, batch_size):
  """Calculates the loss from the logits and the labels.

  Args:
    regression_score
    labels: Labels tensor, int32 - [batch_size, NUM_CLASSES].

  Returns:
    loss: Loss tensor of type float.
  """
  with tf.name_scope('Loss'):
      # labels = tf.to_float(labels) # old version
      labels = tf.cast(labels, dtype="float32")
      # MSQE
      mean_squared_error = tf.losses.mean_squared_error(labels, regression_score, reduction=tf.losses.Reduction.SUM)
      # mean over batch size
      mean_squared_error = mean_squared_error/batch_size
      # add summarizers if training case (loss for test is reported per epoch)
      tf.summary.scalar('mean_squared_error', mean_squared_error)

      # add regularizer:
      if l2_regularization_strength == 0:
          return mean_squared_error
      else:
          # L2 regularization for all trainable parameters
          l2_loss = tf.add_n([tf.nn.l2_loss(v)
                              for v in tf.trainable_variables()
                              if not('bias' in v.name)])
          # Add the regularization term to the loss
          total_loss = (mean_squared_error + l2_regularization_strength * l2_loss)
          # add summarizers
          tf.summary.scalar('l2_loss', l2_loss)
          tf.summary.scalar('total_loss', total_loss)

          return total_loss

def loss_test(regression_score, labels, batch_size):
  """Calculates the loss from the logits and the labels for training case without summaries.
  Args:
    regression_score
    labels: Labels tensor, int32 - [batch_size, NUM_CLASSES].
  Returns:
    loss: Loss tensor of type float.
  """
  with tf.name_scope('Test_Loss'):
      # labels = tf.to_float(labels)
      labels = tf.cast(labels, dtype="float32")
      test_mean_squared_error = tf.losses.mean_squared_error(labels, regression_score, reduction=tf.losses.Reduction.SUM)
      test_mean_squared_error = test_mean_squared_error/batch_size
      return test_mean_squared_error

def training(loss, learning_rate, beta_1, beta_2, epsilon, global_step):
  """Sets up the training Operations.
  Creates a summarizer to track the loss over time in TensorBoard.
  Creates an optimizer and applies the gradients to all trainable variables.
  The Op returned by this function is what must be passed to the
  `sess.run()` call to cause the model to train.
  Args:
    loss: Loss tensor, from loss().
    learning_rate: The learning rate to use for gradient descent.
  Returns:
    train_op: The Op for training.
  """
  # with Learning Rate decay
  # learning_rate = tf.train.exponential_decay(learning_rate, global_step, learning_rate_decay_steps, 0.96)
  optimizer = tf.train.AdamOptimizer(
    learning_rate = learning_rate,
    beta1 = beta_1,
    beta2 = beta_2,
    epsilon = epsilon)
  trainables = tf.trainable_variables()
  train_op = optimizer.minimize(loss, var_list=trainables, global_step=global_step)

  return train_op

# def evaluation(logits, labels):
#   """Evaluate the quality of the logits at predicting the label.
#
#   Args:
#     logits: Logits tensor, float - [batch_size, NUM_CLASSES].
#     labels: Labels tensor, int32 - [batch_size], with values in the
#       range [0, NUM_CLASSES).
#
#   Returns:
#     A scalar int32 tensor with the number of examples (out of batch_size)
#     that were predicted correctly.
#   """
#   # For a classifier model, we can use the in_top_k Op.
#   # It returns a bool tensor with shape [batch_size] that is true for
#   # the examples where the label is in the top k (here k=1)
#   # of all logits for that example.
#   # correct = tf.nn.in_top_k(logits, tf.argmax(labels, 1), 1)
#   # return tf.reduce_sum(tf.cast(correct, tf.int32))
#   correct_prediction = tf.equal(tf.argmax(logits,1), tf.argmax(labels,1))
#   mean_correct = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
#
#   return mean_correct
