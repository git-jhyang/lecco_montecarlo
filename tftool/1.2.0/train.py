
from input import Input, print_input
import aenettool as ant
from tfmodel import MODEL
import os

class TRAIN:
    def __init__(self, fn='train.in'):
        self.isinit = True
        self.initialize(fn)

    def initialize(self, fn):
        self.Inp = Input(fn)
        if os.path.isfile(self.Inp.train['TrnSet_fn']+'.scaled'):
            self.Data = ant.read_aenet_genout(self.Inp.train['TrnSet_fn']+'.scaled')
        else:
            Data = ant.read_aenet_genout(self.Inp.train['TrnSet_fn'])
            if not Data.Hdr['isNorm'][0]:
                self.Data = ant.normalize(Data, E_cut=self.Inp.train['E_cut'])
                ant.write_aenet_genout(self.Data)
        self.Inp.get_InpNetwork(Atoms = self.Data.Hdr['Atoms'],
                                nsfs  = [self.Data.Stps[iType]['nsf'] 
                                            for iType in range(self.Data.Hdr['nType'])])
        self.NNP = ant.build_nnp(network  = self.Inp.network,
                                 setup    = self.Data.Stps,
                                 trainset = self.Data.Hdr)
        print_input(self.Inp)

    def build_model(self):
        self.M = MODEL(flot_prec    = self.Inp.train['float_prec'],
                       NNP          = self.NNP,
                       nType        = self.Data.Hdr['nType'],
                       nLays        = [self.Inp.network['Nets'][iType]['nLay']
                                            for iType in range(self.Data.Hdr['nType'])],
                       nNodes       = [self.Inp.network['Nets'][iType]['nNodes']
                                            for iType in range(self.Data.Hdr['nType'])],
                       ActFuncs     = [self.Inp.network['Nets'][iType]['ActFunc']
                                            for iType in range(self.Data.Hdr['nType'])],

                       E_Atom       = self.Data.Hdr['E_Atom'],
                       E_scale      = self.Data.Hdr['E_scale'],
                       E_shift      = self.Data.Hdr['E_shift'],

                       nStructure   = self.Data.Hdr['nStru'],
                       test_ratio   = self.Inp.train['test_ratio'],
                       dataset      = self.Data.Strs,
                       shuffle_data = self.Inp.train['shuffle'],
                       random_state = self.Inp.train['random_state'],

                       optimizer    = self.Inp.train['method'],

                       output_path  = self.Inp.train['output_path'])

    def train_model(self):
        self.M.run(epoch_max     = self.Inp.train['epoch_max'],
                   epoch_ckpt    = self.Inp.train['epoch_ckpt'],
                   batch_size    = self.Inp.train['batch_size'],
                   shuffle_batch = self.Inp.train['shuffle'])
