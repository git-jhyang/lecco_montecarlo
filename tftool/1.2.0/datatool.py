import pyutil.util as pu
import numpy as np
import tensorflow as tf
from time import time

class Structure_Data:
    def __init__(self, output_path, nType, E_Atom, nsfs, nStructure, 
                 dataset, fdtype, random_state):
        t1 = time()
        self.output_path = output_path

        if fdtype == tf.float16:
            self.fdtype = np.float16
        elif fdtype == tf.float64:
            self.fdtype = np.float64
        else:
            self.fdtype = np.float32

        self.random_state = random_state

        self.nGroup     = np.max(dataset['Group']) + 1
        self.mask_group = []
        self.fns        = dataset['fn']

        self.data = {'energy/E_label:0':dataset['E'],
                     'input/E_atomic:0':np.zeros((self.nGroup),dtype=self.fdtype),
                     'input/nAtom:0':np.zeros((self.nGroup),dtype=self.fdtype)}
        
        for iType in range(nType):
            self.data['input/Type_{:02d}:0'.format(iType)] = np.ndarray((nStructure),dtype=np.object)

        for iGroup in range(self.nGroup):
            self.mask_group.append(dataset['Group'] == iGroup)
            nSpec = dataset['nSpec'][self.mask_group[iGroup]][0]
            nAtom = dataset['nAtom'][self.mask_group[iGroup]][0]
            self.data['input/E_atomic:0'][iGroup] = np.sum(np.multiply(nSpec,E_Atom))
            self.data['input/nAtom:0'][iGroup]    = nAtom
        
        for iStru in range(nStructure):
            for iType in range(nType):
                mask_type = dataset['Type'][iStru] == iType
                sfdata    = dataset['sfval'][iStru][mask_type]
                if sfdata.shape[0] == 0:
                    sfdata = np.zeros((0,nsfs[iType]), dtype=self.fdtype)
                else:
                    sfdata = np.stack(sfdata)
                self.data['input/Type_{:02d}:0'.format(iType)][iStru] = sfdata
        
        print(' Structure Dataset initialized : {:.2f} sec'.format(time()-t1))
#        for data in self.data_list:
#            for key, item in data.items():
#                print(key, type(item), item.shape)
#                print(item)

    def train_valid_split(self, test_ratio, shuffle):
        self.nTrain = np.zeros((self.nGroup), dtype=np.int32)
        self.nValid = np.zeros((self.nGroup), dtype=np.int32)
        self.train_list = {}
        self.valid_list = {}
        train_masks     = []

        for iGroup in range(self.nGroup):
            self.nValid[iGroup] = int(np.sum(self.mask_group[iGroup])*test_ratio)
            self.nTrain[iGroup] = np.sum(self.mask_group[iGroup]) - self.nValid[iGroup]
            mask    = np.ndarray(np.sum(self.mask_group[iGroup]), dtype=np.bool)
            mask[:] = False
            mask[:self.nTrain[iGroup]] = True

            if shuffle: 
                np.random.shuffle(mask)
            train_masks.append(mask)

            train_dict = {}
            valid_dict = {}

            for key, item in self.data.items():
                if item.shape[0] == self.nGroup:
                    train_dict[key] = np.array([item[iGroup]])
                    valid_dict[key] = np.array([item[iGroup]])
                else:
                    item_group = item[self.mask_group[iGroup]]
                    train_dict[key] = np.stack(item_group[mask])
                    valid_dict[key] = np.stack(item_group[~mask])

            self.train_list[iGroup] = train_dict
            self.valid_list[iGroup] = valid_dict

        with open(self.output_path+'list.train','w') as f:
            for iGroup in range(self.nGroup):
                for fn in self.fns[self.mask_group[iGroup]][train_masks[iGroup]]:
                    f.write('{}\n'.format(fn))

        with open(self.output_path+'list.valid','w') as f:
            for iGroup in range(self.nGroup):
                for fn in self.fns[self.mask_group[iGroup]][~train_masks[iGroup]]:
                    f.write('{}\n'.format(fn))

    def prepare_batch(self, batch_size, shuffle):
        def _get_batch_data(nGroup, batch_sizes, data, fdtype):
            batch = []
            for iGroup in range(nGroup):
                idx2 = 0
                for iBatch, batch_size in enumerate(batch_sizes[iGroup]):
                    idx1  = idx2
                    idx2 += batch_size
                    data_dict = {}
                    for key, item in data[iGroup].items():
                        if item.shape[0] == 1:
                            data_dict[key] = item
                        else:
                            data_dict[key] = item[idx1:idx2]
                    batch.append(data_dict)
            return batch

        max_batch_sizes = [np.min(np.hstack([self.nTrain, batch_size])),
                           np.min(np.hstack([self.nValid, batch_size]))]
        nBatch       = np.zeros((2,self.nGroup), dtype=np.int32)
        batch_remain = np.zeros((2,self.nGroup), dtype=np.int32)
        for iGroup in range(self.nGroup):
            nBatch[0,iGroup] = np.int32(self.nTrain[iGroup]/max_batch_sizes[0]) + 1
            nBatch[1,iGroup] = np.int32(self.nValid[iGroup]/max_batch_sizes[1]) + 1
            batch_remain[0,iGroup] = np.mod(self.nTrain[iGroup], max_batch_sizes[0])
            batch_remain[1,iGroup] = np.mod(self.nValid[iGroup], max_batch_sizes[1])

            for i in range(2):
                if batch_remain[i,iGroup] == 0:
                    nBatch[i,iGroup]      -= 1
                    batch_remain[i,iGroup] = max_batch_sizes[i]

        batch_sizes = [[[max_batch_sizes[i] for iBatch in range(nBatch[i,iGroup])] 
                            for iGroup in range(self.nGroup)] for i in range(2)]

        for i in range(2):
            for iGroup in range(self.nGroup):
                batch_sizes[i][iGroup][-1] = batch_remain[i][iGroup]

        self.batch_idx = np.array([np.arange(0, np.sum(nBatch[0]), dtype=np.int32),
                                   np.arange(0, np.sum(nBatch[1]), dtype=np.int32)])
        if shuffle:
            for i in range(2):
                np.random.shuffle(self.batch_idx[i])

        self.train_batch = _get_batch_data(nGroup      = self.nGroup, 
                                           batch_sizes = batch_sizes[0], 
                                           data        = self.train_list, 
                                           fdtype      = self.fdtype)

        self.valid_batch = _get_batch_data(nGroup      = self.nGroup, 
                                           batch_sizes = batch_sizes[1], 
                                           data        = self.valid_list, 
                                           fdtype      = self.fdtype)
        return np.sum(nBatch[0]), np.sum(nBatch[1])

    def get_next_batch(self, iBatch, isTrain):
        if isTrain:
            return self.train_batch[self.batch_idx[0][iBatch]]
        else:
            return self.valid_batch[self.batch_idx[1][iBatch]]
