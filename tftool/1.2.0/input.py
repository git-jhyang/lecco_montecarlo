
import numpy as np
import pyutil.util as pu
from activation import conv_to_idx, conv_to_str
import os
from datetime import datetime

class Input:
    def __init__(self, fn):
        self.fn      = fn
        pu.assert_file(self.fn,'input.get_inputs.InpTrain()')
        date_now = str(np.datetime64(datetime.now())).replace(':','').replace('-','')[:13]
        date_now = date_now.replace('T','_T_')
        self._get_defaults(date_now)
        with open(self.fn,'r') as f:
            trnset_fn = self._get_line(f,'TRAININGSET')
            if trnset_fn != None:
                self.train['TrnSet_fn'] = trnset_fn[0][1]
            else:
                pu.call_error('','')

            testrate = self._get_line(f,'TESTPERCENT')
            if testrate != None:
                self.train['TestRate']  = np.float32(testrate[0][1])
                if self.train['TestRate'] > 1: self.train['TestRate'] *= 0.01

            e_cut = self._get_line(f,'MAXENERGY')
            if e_cut != None:
                self.train['E_cut'] = np.float64(e_cut[0][1])

            method = self._get_line(f,'METHOD',1)
            if method != None:
                self.train['method']['algorithm'] = method[0][0]
                if len(method[0]) != 1:
                    for string in method[0][1:]:
                        comm = string.split('=')
                        if comm[0] == 'learning_rate':
                            self.train['method']['learning_rate'] = np.float32(comm[1])
                        # add training algorithm options

            epoch_max = self._get_line(f,'EPOCH.MAX')
            if epoch_max != None:
                self.train['epoch_max'] = np.int32(epoch_max[0][1])

            epoch_ckpt = self._get_line(f,'EPOCH.CKPT')
            if epoch_ckpt != None:
                self.train['epoch_ckpt'] = np.int32(epoch_ckpt[0][1])
                if self.train['epoch_ckpt'] < 1: self.train['epoch_ckpt'] = 1

            batch_size = self._get_line(f,'BATCH.SIZE')
            if batch_size != None:
                self.train['batch_size'] = np.int32(batch_size[0][1])

            float_type = self._get_line(f,'FLOAT.PREC')
            if float_type != None:
                self.train['float_prec'] = np.int32(float_type[0][1])

            shuffle = self._get_line(f,'SHUFFLE.DATA')
            if shuffle != None:
                self.train['shuffle'] = np.bool(shuffle[0][1])

            output_path = self._get_line(f,'OUTPUT.PATH')
            if output_path != None:
                self.train['output_path'] = '{}/{}'.format(output_path[0][1],date_now)
            self.train['output_path'] = self._bulid_directories(self.train['output_path'])

    def get_InpNetwork(self, Atoms=None, nsfs=None):
        if Atoms is None:
            pu.call_error('input error - input.get_inputs.InpNetwork()',
                '\'Atoms\' is not given - string array')
        Atoms = [Atom.strip() for Atom in Atoms]
        if nsfs is None:
            pu.call_error('input error - input.get_inputs.InpNetwork()',
                '\'nsfs\' is not given - integer array')
        if np.shape(Atoms) != np.shape(nsfs):
            pu.call_error('input error - input.get_inputs.InpNetwork()',
                'two inputs have differentT. shape - Atoms (%s), nsfs (%s)'%(
                np.shape(Atoms), np.shape(nsfs)))
        nType = len(Atoms)
        self.network['fns']   = np.ndarray((nType),dtype=np.object)
        self.network['Atoms'] = Atoms
        self.network['Nets']  = np.ndarray((nType),dtype=np.object)
        with open(self.fn,'r') as f:
            network_setting  = self._get_line(f,'NETWORKS',nType)
        if network_setting is None:
            pu.call_error('input error - input.get_inputs.InpNetwork()',
                'unable to read NETWORKS')
        for setting in network_setting:
            try:
                iType = Atoms.index(setting[0])
            except ValueError:
                pu.call_error('input error - input.get_inputs.InpNetwork()',
                    'Atom \'%s\' is missing in TRAININGSET'%(setting[0]))
            self.network['fns'][iType] = self.train['output_path']+'nnp/'+setting[1]
            Net = {}
            nsf = nsfs[iType]
            Net['nLay']       = np.int32(setting[2])+2
            Net['nW']         = np.int32(0)
            Net['nNode']      = np.int32(0)
            Net['nNodes']     = np.ones((Net['nLay']),dtype=np.int32)
            Net['nNodes'][0]  = nsf
            Net['nNodes_max'] = np.int32(0)
            Net['ActFunc']    = np.zeros((Net['nLay']-1),dtype=np.int32)
            Net['iW_max']     = np.zeros((Net['nLay']),dtype=np.int32)
            Net['iNode_max']  = np.zeros((Net['nLay']),dtype=np.int32)
            Net['W']          = np.ndarray((Net['nLay']-1),dtype=np.object)

            for iLay, Lay in enumerate(setting[3:1+Net['nLay']]):
                nNode, ActFunc        = Lay.split(':')
                Net['nNodes'][iLay+1] = nNode
                Net['ActFunc'][iLay]  = conv_to_idx(ActFunc=ActFunc)

            Net['nNodes_max'] = np.int32(np.max(Net['nNodes']))

            for iLay in range(Net['nLay']-1):
                nNode1 = Net['nNodes'][iLay]
                nNode2 = Net['nNodes'][iLay+1]
                Net['nNode']            += nNode1 + 1
                Net['nW']               += (nNode1 + 1)*nNode2
                Net['iNode_max'][iLay+1] = Net['nNode']
                Net['iW_max'][iLay+1]    = Net['nW']
                Net['W'][iLay]           = np.zeros((nNode1+1,nNode2),dtype=np.float64)
            
            Net['nNode']               += Net['nNodes'][-1]
            self.network['Nets'][iType] = Net

    def _get_defaults(self, date_now):
        self.network = {}
        self.train   = {'TrnSet_fn'     : None,
                        'output_path'   : 'outputs/{}'.format(date_now),
                        'test_ratio'    : 0.25,
                        'E_cut'         : 1e7,
                        'method'        : {'algorithm'    :'Adam',
                                           'learning_rate':None},
                        'epoch_max'     : 20000,
                        'epoch_ckpt'    : 100,
                        'batch_size'    : 100,
                        'float_prec'    : 32,
                        'shuffle'       : True,
                        'random_state'  : None}

    def _get_line(self, fo, string, Len=0):
        count = 0
        while True:
            count += 1
            if count > 1000:
                fo.seek(0)
                return None
            line = fo.readline()
            if string not in line.upper(): continue
            if line[0:1] in ('!','#','%'): continue
            break

        if Len == 0:
            fo.seek(0)
            return [line.split()]
        
        else:
            count = 0
            lines = []
            while len(lines) < Len:
                if count > 1000:
                    fo.seek(0)
                    return None
                count += 1
                line = fo.readline().split()
                if len(line) == 0: continue
                if line[0][0:1] in ('!','#','%'): continue
                lines.append(line)
            fo.seek(0)
            return lines

    def _bulid_directories(self, path_in):
        dirs = np.array(path_in.split('/'))
        path = '.'
        if dirs[0] == '': path = ''

        dirs = dirs[(dirs != '') & (dirs != '.')]
        for dir in dirs[:-1]:
            path += '/' + dir
            if not os.path.isdir(path): 
                os.mkdir(path)

        path += '/'+dirs[-1]
        if os.path.isdir(path):
            for i in range(99):
                if not os.path.isdir(path+'_{:02d}'.format(i)):
                    path += '_{:02d}'.format(i)
                    break

        path += '/'
        os.mkdir(path)
        os.mkdir(path+'energy')
        os.mkdir(path+'nnp')
        os.mkdir(path+'ckpt')
        return path

def print_input(Inp):
    with open(Inp.train['output_path']+'train.in','w') as f:
        f.write('TRAININGSET    {}\n'.format(Inp.train['TrnSet_fn']))
        f.write('TESTPERCENT    {}\n\n'.format(Inp.train['TestRate']*100))
        f.write('MAXENERGY      {}\n\n'.format(Inp.train['E_cut']))
        f.write('METHOD\n')
        f.write('{}\n\n'.format(Inp.train['method']['algorithm']))
        f.write('EPOCH.MAX      {}\n'.format(Inp.train['epoch_max']))
        f.write('EPOCH.CKPT     {}\n\n'.format(Inp.train['epoch_ckpt']))
        f.write('BATCH.SIZE     {}\n\n'.format(Inp.train['batch_size']))
        f.write('FLOAT.PREC     {}\n'.format(Inp.train['float_prec']))
        f.write('SHUFFLE.DATA   {}\n\n'.format(Inp.train['shuffle']))
        f.write('OUTPUT.PATH    {}\n\n'.format(Inp.train['output_path']))
        f.write('NETWORKS\n')
        for Atom, fn, Net in zip(Inp.network['Atoms'], Inp.network['fns'], Inp.network['Nets']):
            f.write('  {:4s} {:15s} {:4d}'.format(Atom, fn.split('/')[-1], Net['nLay']-2))
            for nNode, ActFunc in zip(Net['nNodes'][1:-1], Net['ActFunc'][:-1]):
                f.write('   {}:{}'.format(nNode, conv_to_str(ActFunc)))
            f.write('\n')
        
        

