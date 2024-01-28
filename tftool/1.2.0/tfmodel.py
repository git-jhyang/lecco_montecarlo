import tensorflow as tf
import numpy as np
import pyutil.util as pu
from activation import activate, conv_to_str
from datatool import Structure_Data
from aenettool import update_network

class MODEL:
    def __init__(self, flot_prec, NNP,
                 nType, nLays, nNodes, ActFuncs, 
                 E_Atom, E_scale, E_shift, 
                 nStructure, test_ratio, dataset,
                 shuffle_data, random_state, 
                 optimizer, output_path):

        if flot_prec == 16:
            self.fdtype = tf.float16
        elif flot_prec == 64:
            self.fdtype = tf.float64
        else:
            self.fdtype = tf.float32

        self.NNP         = NNP
        self.output_path = output_path

        self.SD = Structure_Data(output_path  = self.output_path,
                                 nType        = nType,
                                 E_Atom       = E_Atom,
                                 nsfs         = [nNode[0] for nNode in nNodes],
                                 nStructure   = nStructure,
                                 dataset      = dataset,
                                 fdtype       = self.fdtype,
                                 random_state = random_state)

        self.SD.train_valid_split(test_ratio  = test_ratio, 
                                  shuffle     = shuffle_data)

        output_energies = [self._element_network(Type    = iType,
                                                 nLay    = nLays[iType],
                                                 nNodes  = nNodes[iType],
                                                 ActFunc = ActFuncs[iType])
                                                 for iType in range(nType)]
        self.summ_net = tf.summary.merge_all()        

        with tf.variable_scope('input/'):
            inv_E_scale = tf.constant(1/E_scale, dtype=self.fdtype)
            E_shift     = tf.constant(E_shift, dtype=self.fdtype)    
            E_atomic    = tf.placeholder(self.fdtype, shape=[1,], name='E_atomic')
            self.nAtom  = tf.placeholder(self.fdtype, shape=[1,], name='nAtom')

        with tf.name_scope('energy'):
            self.E_label  = tf.placeholder(self.fdtype, shape=[None,], name='E_label')
            self.E_pred   = tf.add_n(output_energies, name='E_pred')
            self.dE_out   = tf.multiply(tf.add_n([-self.E_label, self.E_pred]),
                                        1/self.nAtom, name='dE_out')
            self.E_coh    = tf.add(tf.multiply(self.E_label, inv_E_scale), 
                                   tf.multiply(E_shift, self.nAtom), name='E_coh')
            self.E_tot    = tf.add(self.E_coh, E_atomic, name='E_tot')
            self.E_coh_nn = tf.add(tf.multiply(self.E_pred, inv_E_scale), 
                                   tf.multiply(E_shift, self.nAtom), name='E_coh_nn')
            self.E_tot_nn = tf.add(self.E_coh_nn, E_atomic, name='E_tot_nn')


        with tf.name_scope('loss'):
            dE_in     = tf.placeholder(self.fdtype, shape=[None,], name='dE')

            self.MAE  = tf.squeeze(tf.multiply(tf.reduce_mean(tf.abs(dE_in)), 
                                               inv_E_scale), name='MAE')
            self.SSE  = tf.squeeze(tf.reduce_mean(tf.multiply(dE_in, dE_in)), name='SSE')
            self.RMSE = tf.squeeze(tf.multiply(tf.sqrt(self.SSE), inv_E_scale), name='RMSE')

        sum_dE   = tf.summary.histogram('dE', dE_in)
        sum_MAE  = tf.summary.scalar('MAE', self.MAE)
        sum_SSE  = tf.summary.scalar('SSE', self.SSE)
        sum_RMSE = tf.summary.scalar('RMSE', self.RMSE)
        self.summ_loss = tf.summary.merge([sum_dE, sum_MAE, sum_SSE, sum_RMSE])

        with tf.name_scope('optimizer'):
            opt      = self._get_optimizer(optimizer)
            self.opt = opt.minimize(tf.reduce_mean(tf.abs(self.dE_out)), name='minimize')

    def run(self, epoch_max, epoch_ckpt, batch_size, shuffle_batch):
        self.saver = tf.train.Saver(max_to_keep=None)
        log_ckpt   = int(epoch_ckpt/5)
        if log_ckpt == 0: log_ckpt = epoch_ckpt
        
        writer_train = tf.summary.FileWriter(self.output_path+'summary'+'/train')
        writer_valid = tf.summary.FileWriter(self.output_path+'summary'+'/valid')

        with tf.Session() as sess:
            writer_train.add_graph(sess.graph)
            sess.run(tf.global_variables_initializer())

            print()
            print('         {:30s} |   {:30s}'.format('Training set','Validation set'))
            print(' Epoch   {:14s}  {:14s} |   {:14s}  {:14s}'.format(
                    'MAE','RMSE','MAE','RMSE'))

            for iEpoch in range(epoch_max+1):
                nBatchT, _ = self.SD.prepare_batch(batch_size = batch_size,
                                                   shuffle    = shuffle_batch)

                for iBatch in range(nBatchT):
                    train_batch = self.SD.get_next_batch(iBatch=iBatch, isTrain=True)
                    sess.run(self.opt, feed_dict=train_batch)

                if np.mod(iEpoch, log_ckpt) == 0 or iEpoch == epoch_max:
                    nBatchT, nBatchV = self.SD.prepare_batch(batch_size = batch_size,
                                                             shuffle    = False)

                    summ_net, loss_train = self._save_state(sess    = sess, 
                                                            iEpoch  = iEpoch, 
                                                            nBatch  = nBatchT, 
                                                            isTrain = True)

                    _, loss_valid = self._save_state(sess    = sess,
                                                     iEpoch  = iEpoch,
                                                     nBatch  = nBatchV,
                                                     isTrain = False)

                    writer_train.add_summary(summ_net,   iEpoch)
                    writer_train.add_summary(loss_train, iEpoch)
                    writer_valid.add_summary(loss_valid, iEpoch)

                if np.mod(iEpoch, epoch_ckpt) == 0 or iEpoch == epoch_max:
                    nBatchT, nBatchV = self.SD.prepare_batch(batch_size = batch_size,
                                                             shuffle    = False)
                    self._print_state(sess     = sess, 
                                      iEpoch   = iEpoch, 
                                      nBatches = (nBatchT, nBatchV))

    def _save_state(self, sess, iEpoch, nBatch, isTrain):
        for iBatch in range(nBatch):
            batch = self.SD.get_next_batch(iBatch=iBatch, isTrain=isTrain)
            summ_net, dE_out = sess.run([self.summ_net, self.dE_out], feed_dict=batch)
            if iBatch == 0:
                dE = dE_out
            else:
                dE = np.hstack([dE, dE_out])
        summ_loss = sess.run(self.summ_loss, feed_dict={'loss/dE:0':dE})
        return summ_net, summ_loss

    def _print_state(self, sess, iEpoch, nBatches):
        def _write_energy(fobj, energies, nAtoms):
            for e in energies:
                fobj.write('{:12.6f} {:12.6f} {:12.6f} {:12.6f} {:14.9f} {:14.9f}\n'.format(
                    e[0], e[1], e[2], e[3], e[2]/nAtoms[0], e[3]/nAtoms[0]))

        MAE     = [0, 0]
        RMSE    = [0, 0]
        isTrain = [True,False]
        state   = ['train','valid']
    
        for i, nBatch in enumerate(nBatches):
            with open(self.output_path+'energy/{:06d}_{}.txt'.format(iEpoch, state[i]),'w') as f:
                f.write('# {:10s} {:12s} {:12} {:12} {:14} {:14}\n'.format(
                    'E_tot','E_tot_nn','E_coh','E_coh_nn','E_coh/at','E_coh_nn/at'))
                for iBatch in range(nBatch):
                    batch = self.SD.get_next_batch(iBatch=iBatch, isTrain=isTrain[i])
                    dE_out, E1, E2, E3, E4, nAtoms = sess.run(
                            [self.dE_out, self.E_tot, self.E_tot_nn, 
                            self.E_coh, self.E_coh_nn, self.nAtom], feed_dict=batch)
                    if iBatch == 0:
                        dE = dE_out
                    else:
                        dE = np.hstack([dE, dE_out])
                    _write_energy(f, np.stack([E1, E2, E3, E4], axis=1), nAtoms)
            MAE[i], RMSE[i] = sess.run([self.MAE, self.RMSE], feed_dict={'loss/dE:0':dE})

        print('{:6d}   {:.6e}    {:.6e}   |   {:.6e}    {:.6e}'.format(
                iEpoch, MAE[0], RMSE[0], MAE[1], RMSE[1]))

        data_structure = tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES)
        data_values    = sess.run(data_structure)

        self.saver.save(sess, self.output_path+'ckpt/model', global_step=iEpoch)
        self.NNP = update_network(self.NNP, iEpoch, data_structure, data_values)

    def _get_optimizer(self, optimizer):
        algo   = optimizer['algorithm']
        l_rate = optimizer['learning_rate']

        if algo == 'Adam':
            if l_rate is None:
                opt = tf.train.AdamOptimizer()
            else:
                opt = tf.train.AdamOptimizer(learning_rate=l_rate)
        elif algo in ('GradientDescent','GD'):
            if l_rate is None:
                opt = tf.train.GradientDescentOptimizer()
            else:
                opt = tf.train.GradientDescentOptimizer(learning_rate=l_rate)
        else:
            pu.call_error('_','not supported optimizer - {}'.format(algo))
        return opt

    def _element_network(self, Type, nLay, nNodes, ActFunc):
        def _get_random_weight(nNode1, nNode2, iLay, nLay):
            #if iLay != nLay-1:
            #    rand = tf.random_uniform(shape=[nNode1,nNode2], maxval=1, dtype=self.fdtype)
            #else:
            rand = tf.random_normal(shape=[nNode1,nNode2], dtype=self.fdtype)
            return tf.Variable(rand, name='W')

        def _get_random_bias(nNode2, iLay, nLay):
            #if iLay != nLay-1:
            #    rand = tf.random_uniform(shape=[nNode2], maxval=1, dtype=self.fdtype)
            #else:
            rand = tf.random_normal(shape=[nNode2], dtype=self.fdtype)
            return tf.Variable(rand, name='B')

        with tf.name_scope('input/'):
            layer = tf.placeholder(self.fdtype, shape=[None, None, nNodes[0]],
                                   name='Type_{:02d}'.format(Type))

        with tf.name_scope('Type_{:02d}'.format(Type)):
            for iLay in range(1,nLay):
                with tf.name_scope('layer_{}'.format(iLay)):
                    nNode1 = nNodes[iLay-1]
                    nNode2 = nNodes[iLay]

                    W      = _get_random_weight(nNode1, nNode2, iLay, nLay)
                    layer  = tf.tensordot(layer, W, axes=1, name='W_mult')
                    b      = _get_random_bias(nNode2, iLay, nLay),
                    layer  = tf.add(layer, b, name='b_add')

                    layer  = activate(layer, ActFunc[iLay-1])

                    tf.summary.histogram('0_W', W)
                    tf.summary.histogram('1_b', b)

            output = tf.squeeze(tf.reduce_sum(layer, axis=1), name='output_energy')

            return output
