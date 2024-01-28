"""
aenet_tool
------------------------------------------------------------------
This module interfaces AENET with Python3.

    .read_aenet_genout(filename)
    .write_aenet_genout(data)
    .get_stp_scaler(Stps)
    .normalize(data, E_max=1e10)
    .build_nnp(inp=None)
    .save_networks(NNP, iIter=None)
    .load_networks(fns)

"""
import numpy as np
import time
from scipy.io import FortranFile
import pyutil.util as pu

def _convert_str_to_fbytes(x, Len=None):
    """
_convert_str_to_fbytes(x, Len)
------------------------------------------------------------------
take string or list (including numpy.ndarray and numpy.str) 'x' 
and return combined numpy.fbytes_ data.
if 'Len' is given, the length of resulting data adjusted to 'Len'

    x               Len     Return
    ---------------------------------------
    'a bc'          None    b'a bc'
    'a bc'          10      b'a bc      '
    ['a',' b','c']  None    b'a bc'
    ['a',' b','c']  10      b'a bc      '

    """
    if 'str' in str(type(x)):
        string = x
    elif ('list' in str(type(x))) or ('array' in str(type(x))):
        string = ''
        for item in x:
            string += item
    else:
        pu.call_error('input error - aenet_tool._convert_str_to_fbytes',
            'not a list or string - '+str(type(x)))
    if pu.assert_num(Len):
        string = string.ljust(int(Len))
    string = np.array([string],dtype=np.bytes_)
    return string

def _read_fbytes_header(fo):
    """
_read_fbytes_header(fo)
------------------------------------------------------------------
read header data of AENET GENERATION OUTPUTFILE written by fortran 
binary format by using scipy.io.FortranFile().
'fo' is a file object of scipy.io.FortranFile()
    """
    Hdr = {}
    Hdr['nType']   = fo.read_ints('i4')[0]
    Hdr['nStru']   = fo.read_ints('i4')[0]
    Hdr['Atoms']   = fo.read_record('a2').astype(np.str)
    Hdr['E_Atom']  = fo.read_reals('f8')
    Hdr['isNorm']  = fo.read_record('?')
    Hdr['E_scale'] = fo.read_reals('f8')[0]
    Hdr['E_shift'] = fo.read_reals('f8')[0]
    return Hdr

def _write_fbytes_header(fo, Hdr):
    """
_write_fbytes_header(fo, Hdr)
------------------------------------------------------------------
write header data of AENET GENERATION OUTPUTFILE saved as 'Hdr'
in fortran binary format by using scipy.io.FortranFile().
'fo' is a file object of scipy.io.FortranFile()
    """
    fo.write_record(np.int32(Hdr['nType']))
    fo.write_record(np.int32(Hdr['nStru']))
    fo.write_record(_convert_str_to_fbytes(Hdr['Atoms']))
    fo.write_record(np.float64(Hdr['E_Atom']))
    fo.write_record(Hdr['isNorm'])
    fo.write_record(np.float64(Hdr['E_scale']))
    fo.write_record(np.float64(Hdr['E_shift']))

def _read_fbytes_setup(fo):
    """
_read_fbytes_setup(fo)
------------------------------------------------------------------
read SETUP data of AENET written by fortran binary format by using
scipy.io.FortranFile().
'fo' is a file object of scipy.io.FortranFile()
    """
    Stp = {}
    Stp['Des']       = fo.read_record('a1024').astype(np.str)[0].replace('$','').strip()
    Stp['Atom']      = fo.read_record('a2').astype(np.str)[0]
    Stp['nEnv']      = fo.read_ints('i4')[0]
    Stp['Atoms']     = fo.read_record('a2').astype(np.str)
    Stp['r_min']     = fo.read_reals('f8')[0]
    Stp['r_max']     = fo.read_reals('f8')[0]
    Stp['sfType']    = fo.read_record('a100').astype(np.str)[0].strip()
    Stp['nsf']       = fo.read_ints('i4')[0]
    Stp['nsfParam']  = fo.read_ints('i4')[0]
    Stp['sf']        = fo.read_ints('i4')
    Stp['sfParam']   = fo.read_reals('f8').reshape((Stp['nsf'],4))
    Stp['sfEnv']     = fo.read_ints('i4').reshape((Stp['nsf'],2))
    Stp['nEval']     = fo.read_ints('i4')[0]
    Stp['sfval_min'] = fo.read_reals('f8')
    Stp['sfval_max'] = fo.read_reals('f8')
    Stp['sfval_avg'] = fo.read_reals('f8')
    Stp['sfval_cov'] = fo.read_reals('f8')
    return Stp

def _write_fbytes_setup(fo, Stp):
    """
_write_fbytes_setup(fo, Stp)
------------------------------------------------------------------
write SETUP data of AENET saved as 'Stp' in fortran binary format 
by using scipy.io.FortranFile().
'fo' is a file object of scipy.io.FortranFile()
    """
    fo.write_record(_convert_str_to_fbytes(Stp['Des'],Len=1024))
    fo.write_record(_convert_str_to_fbytes(Stp['Atom'],Len=2))
    fo.write_record(np.int32(Stp['nEnv']))
    fo.write_record(_convert_str_to_fbytes(Stp['Atoms']))
    fo.write_record(np.float64(Stp['r_min']))
    fo.write_record(np.float64(Stp['r_max']))
    fo.write_record(_convert_str_to_fbytes(Stp['sfType'],Len=100))
    fo.write_record(np.int32(Stp['nsf']))
    fo.write_record(np.int32(Stp['nsfParam']))
    fo.write_record(np.int32(Stp['sf']))
    fo.write_record(np.float64(Stp['sfParam']))
    fo.write_record(np.int32(Stp['sfEnv']))
    fo.write_record(np.int32(Stp['nEval']))
    fo.write_record(np.float64(Stp['sfval_min']))
    fo.write_record(np.float64(Stp['sfval_max']))
    fo.write_record(np.float64(Stp['sfval_avg']))
    fo.write_record(np.float64(Stp['sfval_cov']))

def _read_fbytes_trnset(fo):
    """
_read_fbytes_trnset(fo)
------------------------------------------------------------------
read TRAIN SET data (TrnSet) of AENET written by fortran binary 
format by using scipy.io.FortranFile().
'fo' is a file object of scipy.io.FortranFile()
    """
    TrnSet = {}
    TrnSet['fn']      = fo.read_record('a1024').astype(np.str)[0].strip()
    TrnSet['isNorm']  = fo.read_record('?')
    TrnSet['E_scale'] = fo.read_reals('f8')[0]
    TrnSet['E_shift'] = fo.read_reals('f8')[0]
    TrnSet['nType']   = fo.read_ints('i4')[0]
    TrnSet['Atoms']   = fo.read_record('a2').astype(np.str)
    TrnSet['E_Atom']  = fo.read_reals('f8')
    TrnSet['nAtom']   = fo.read_ints('i4')[0]
    TrnSet['nStru']   = fo.read_ints('i4')[0]
    TrnSet['E_min'], TrnSet['E_max'], TrnSet['E_avg'] = fo.read_reals('f8')
    return TrnSet

def _write_fbytes_trnset(fo, TrnSet):
    """
_write_fbytes_trnset(fo, TrnSet)
------------------------------------------------------------------
write TRAIN SET data (TrnSet) of AENET saved as 'TrnSet' in 
fortran binary format by using scipy.io.FortranFile().
'fo' is a file object of scipy.io.FortranFile()
    """
    fo.write_record(_convert_str_to_fbytes(TrnSet['fn'],Len=1024))
    fo.write_record(TrnSet['isNorm']) 
    fo.write_record(np.float64(TrnSet['E_scale']))
    fo.write_record(np.float64(TrnSet['E_shift']))
    fo.write_record(np.int32(TrnSet['nType']))
    fo.write_record(_convert_str_to_fbytes(TrnSet['Atoms']))
    fo.write_record(np.float64(TrnSet['E_Atom']))
    fo.write_record(np.int32(TrnSet['nAtom']))
    fo.write_record(np.int32(TrnSet['nStru']))
    fo.write_record(np.float64([TrnSet['E_min'], TrnSet['E_max'], TrnSet['E_avg']]))

def _read_fbytes_network(fo):
    """
_read_fbytes_network(fo)
------------------------------------------------------------------
read NETWORK data of AENET written by fortran binary format by 
using scipy.io.FortranFile().
'fo' is a file object of scipy.io.FortranFile()
    """
    Net = {}
    Net['nLay']       = fo.read_ints('i4')[0]
    Net['nNodes_max'] = fo.read_ints('i4')[0]
    Net['nW']         = fo.read_ints('i4')[0]
    Net['nNode']      = fo.read_ints('i4')[0]
    Net['nNodes']     = fo.read_ints('i4')
    Net['ActFunc']    = fo.read_ints('i4')
    Net['iW_max']     = fo.read_ints('i4')
    Net['iNode_max']  = fo.read_ints('i4')
    W                 = fo.read_reals('f8')

    Net['W'] = []
    iW_1     = 0
    nNode1   = Net['nNodes'][0]
    for iLay in range(Net['nLay']-1):
        iW_2   = Net['iW_max'][iLay+1]
        nNode2 = Net['nNodes'][iLay+1]
        Net['W'].append(W[iW_1:iW_2].reshape((nNode1+1,nNode2)))
        iW_1   = iW_2
        nNode1 = nNode2
    return Net

def _write_fbytes_network(fo, Net):
    """
_write_fbytes_network(fo, Net)
------------------------------------------------------------------
write NETWORK data of AENET saved as 'Net' in fortran binary 
format by using scipy.io.FortranFile().
'fo' is a file object of scipy.io.FortranFile()

    """
    fo.write_record(np.int32(Net['nLay']))
    fo.write_record(np.int32(Net['nNodes_max']))
    fo.write_record(np.int32(Net['nW']))
    fo.write_record(np.int32(Net['nNode']))
    fo.write_record(np.int32(Net['nNodes']))
    fo.write_record(np.int32(Net['ActFunc']))
    fo.write_record(np.int32(Net['iW_max']))
    fo.write_record(np.int32(Net['iNode_max']))
    W = np.hstack([np.hstack(Net_W) for Net_W in Net['W']])
    fo.write_record(np.float64(W))

class read_aenet_genout:
    """
read_aenet_genout(filename)
------------------------------------------------------------------
read AENET GENERATION OUTPUTFILE written by fortran binary that is 
open with scipy.io.FortranFile(). The resulting object includes
'.Hdr', '.Strs', and '.Stps'.
    
.Hdr  : dict, collection of header data
    fn            : filename of AENET GENERATION OUTPUT file
    nType         : number of atom type
    nStru         : number of structure
    Atoms         : name of atom species
    E_Atom(nType) : atomic energies of atom species
    isNorm(4)     : True if data is normalized
    E_scale       : scaleing factor of energy
    E_shift       : shift of energy
    nAtom         : number of total atom
    E_avg         : average of cohesive energy
    E_min         : minimum of cohesive energy
    E_max         : maximum of cohesive energy
    isHaveStp(4)  : True if file includes SETUP data

.Strs : dict, collection of structure data
    fn               : filename of structure
    nType            : number of atom type
    nAtom            : number of total atom
    E                : total energy
    Type(nAtom)      : type of each atom (dim : nAtom)
    coord(nAtom,3)   : coordinate of each atom (dim : nAtom)
    force(nAtom,3)   : force of each atom (dim : nAtom)
    sfval(nAtom,nsf) : symmetry function values of each atom

.Stps(nType) : dict, collection of SETUP data
    Type                  : type of SETUP data
    Des                   : description 
    Atom                  : name of central atom
    nEnv                  : number of environment
    Atoms                 : name of environment atoms
    r_min                 : minimum distance btw atoms
    r_max                 : maximum distance btw atoms
    sfType                : type of symmetry function, Behler only supported
    nsf                   : number of symmetry functions
    nsfParam              : maximum number of symmetry function parameter
    sf(nsf)               : index of symmetry function 
    sfParam(nsf,nsfParam) : parameter of symmetry function 
    sfEnv(nsf,2)          : type of environment atom species 
    nEval                 : number of performed evaluation 
    sfval_min(nsf)        : mininum of symmetry function value 
    sfval_max(nsf)        : maximum of symmetry function value 
    sfval_avg(nsf)        : average of symmetry function value 
    sfval_cov(nsf)        : covariance of symmetry function value

    """
    def __init__(self, f_name):
        t1 = time.time()
        self.Stps = {}
        with FortranFile(f_name,'r') as f:
            self.Hdr       = _read_fbytes_header(f)
            self.Hdr['fn'] = f_name
            self.Strs = {'Group':np.zeros((self.Hdr['nStru']),dtype=np.int32),
                         'nType':np.zeros((self.Hdr['nStru']),dtype=np.int32),
                         'nAtom':np.zeros((self.Hdr['nStru']),dtype=np.int32),
                         'E':np.zeros((self.Hdr['nStru']),dtype=np.float64),
                         'nSpec':np.zeros((self.Hdr['nStru'],self.Hdr['nType']),dtype=np.int32),
                         'fn':np.ndarray((self.Hdr['nStru']),dtype=np.object),
                         'Type':np.ndarray((self.Hdr['nStru']),dtype=np.object),
                         'coord':np.ndarray((self.Hdr['nStru']),dtype=np.object),
                         'force':np.ndarray((self.Hdr['nStru']),dtype=np.object),
                         'sfval':np.ndarray((self.Hdr['nStru']),dtype=np.object)}
            self._read_fbytes_str(f, 0)
            group_ids = [self.Strs['nSpec'][0].tobytes()]

            for iStru in range(1,self.Hdr['nStru']):
                self._read_fbytes_str(f, iStru)
                group_id  = self.Strs['nSpec'][iStru].tobytes()
                group_ids = self._check_group(group_ids, group_id, iStru)
            self.Strs['fn']   = np.stack(self.Strs['fn'])
            self.Hdr['nAtom'] = f.read_ints('i4')[0]
            E = f.read_reals('f8')
            self.Hdr['E_avg'] = E[0]
            self.Hdr['E_min'] = E[1]
            self.Hdr['E_max'] = E[2]
            self.Hdr['isHaveStp'] = f.read_record('?')
            if self.Hdr['isHaveStp'][0]:
                for iType in range(self.Hdr['nType']):
                    f.read_ints('i4')
                    self.Stps[iType] = _read_fbytes_setup(f)
        t2 = time.time()
        print('aenet generator outfile loaded : %.2f sec'%(t2-t1))

    def _read_fbytes_str(self, f, iStru):
        len_file     = f.read_ints('i4')
        fn           = f.read_record('a%i'%(len_file)).astype(np.str)[0]
        nAtom, nType = f.read_ints('i4')

        self.Strs['nType'][iStru] = nType
        self.Strs['nAtom'][iStru] = nAtom
        self.Strs['E'][iStru]     = f.read_reals('f8')[0]
        nSpec = np.zeros((self.Hdr['nType']),dtype=np.int32)
        Type  = np.zeros((nAtom),dtype=np.int32)
        coord = np.zeros((nAtom,3),dtype=np.float64)
        force = np.zeros((nAtom,3),dtype=np.float64)
        sfval = []
        for iAtom in range(nAtom):
            Type[iAtom]  = f.read_ints('i4')[0] - 1
            coord[iAtom] = f.read_reals('f8')
            force[iAtom] = f.read_reals('f8')
            nsf          = f.read_ints('i4')
            sfval.append(f.read_reals('f8'))
            nSpec[Type[iAtom]] += 1

        self.Strs['fn'][iStru]    = fn
        self.Strs['nSpec'][iStru] = nSpec
        self.Strs['Type'][iStru]  = Type
        self.Strs['coord'][iStru] = coord
        self.Strs['force'][iStru] = force
        self.Strs['sfval'][iStru] = np.array(sfval)

    def _check_group(self, group_ids, group_id, iStru):
        for iGroup, group_key in enumerate(group_ids):
            if group_key == group_id:
                self.Strs['Group'][iStru] = iGroup
                return group_ids
        group_ids.append(group_id)
        self.Strs['Group'][iStru] = len(group_ids) - 1
        return group_ids

def write_aenet_genout(Dat_in):
    """
write_aenet_genout(Dat_in)
------------------------------------------------------------------
write AENET GENERATION OUTPUT file in fortran binary file using
class object 'Dat_in' from aenet_tool.read_aenet_genout().
The 'Dat_in' should include '.fn', '.Hdr', '.Strs', and '.Stps'
    """
    t1 = time.time()
    print('aenet writer',Dat_in.Hdr['fn'])
    with FortranFile(Dat_in.Hdr['fn'],'w') as f:
        nsfs = np.array([Dat_in.Stps[iType]['nsf'] 
            for iType in range(Dat_in.Hdr['nType'])], dtype=np.int32)
        _write_fbytes_header(f, Dat_in.Hdr)
        for fn, nType, nAtom, E, nSpec, Type, coord, force, sfval in zip(
                Dat_in.Strs['fn'], Dat_in.Strs['nType'], Dat_in.Strs['nAtom'], 
                Dat_in.Strs['E'], Dat_in.Strs['nSpec'], Dat_in.Strs['Type'], 
                Dat_in.Strs['coord'], Dat_in.Strs['force'], Dat_in.Strs['sfval']):
            f.write_record(np.int32(len(fn)))
            f.write_record(_convert_str_to_fbytes(fn))
            f.write_record(np.int32([nAtom, nType]))
            f.write_record(np.float64([E]))
            for iAtom in range(nAtom):
                iType = Type[iAtom]
                f.write_record(np.int32([iType+1]))
                f.write_record(np.float64(coord[iAtom]))
                f.write_record(np.float64(force[iAtom]))
                f.write_record(np.int32(nsfs[iType]))
                f.write_record(np.float64(sfval[iAtom]))
        f.write_record(np.int32(Dat_in.Hdr['nAtom']))
        f.write_record(np.float64(
            [Dat_in.Hdr['E_avg'], Dat_in.Hdr['E_min'], Dat_in.Hdr['E_max']]))
        if not Dat_in.Hdr['isHaveStp'][0]:
            f.write_record([False for i in range(4)])
            return
        f.write_record([True for i in range(4)])
        for iType in range(Dat_in.Hdr['nType']):
            f.write_record(np.int32([iType+1]))
            _write_fbytes_setup(f, Dat_in.Stps[iType])
    t2 = time.time()
    print('aenet generator outfile is written : %.2f sec'%(t2-t1))

def get_stp_scaler(Stps):
    """
get_stp_scaler(Stps):
------------------------------------------------------------------
calculate scaling factor and shift for all sfvals of input 'Stps'
and return 'sf_scale', 'sf_shift'

sf_scale(nType,nsf) : scaling factor
sf_shift(nType,nsf) : shift value

    """
    sf_scale = {}
    sf_shift = {}
    sf_mask  = {}
    for iType, Stp in Stps.items():
        mask   = np.ones((Stp['nsf']),dtype=np.bool)
        sf_sum = np.add(Stp['sfval_cov'], 
            -np.multiply(Stp['sfval_avg'], Stp['sfval_avg']))
        if min(sf_sum) <= 1e-10:
            pu.call_error('scaling factor warning - aenet_tool.nomalization',
                'Zero scaling factor of symmetry function - Atom : %s '%(
                Stp['Atom']),isexit=False)
            print('   - Excluded Symmetry Function Number : ',end='')
            for i,val in enumerate(sf_sum):
                if val <= 1e-10:
                    mask[i]   = False
                    sf_sum[i] = np.inf
                    print('%5i'%(i),end='')
            print()
        sf_shift[iType] = -Stp['sfval_avg']
        sf_scale[iType] = 1/np.sqrt(sf_sum)
        sf_mask[iType]  = mask
    return sf_scale, sf_shift, sf_mask

def normalize(Dat_in, E_cut):
    """
normalize(Dat_in, E_cut=1e10)
------------------------------------------------------------------
normalize energy of each structure (Dat_in.Strs['E'][0:nStru] and
symm. func. val.s (Dat_in.Strs['sfval'][0:nStru][0:nAtom][0:nsf]).
'E_cut' is cutoff energy for energy normalization.
Energy and symm. func. val.s of each structure nomalized to the 
range from -1 to 1, and the results are return as new dataset 
'Dat_out'.
    """
    import copy
    if '.scaled' in Dat_in.Hdr['fn'] or Dat_in.Hdr['isNorm'][0]:
        return Dat_in

    Dat_out            = copy.deepcopy(Dat_in)
    Dat_out.Hdr['fn'] += '.scaled'
    Strs               = copy.deepcopy(Dat_in.Strs)
    mask_str           = np.ones((Dat_in.Hdr['nStru']),dtype=np.bool)

    if Dat_out.Hdr['E_max'] < E_cut: 
        E_cut = Dat_out.Hdr['E_max']
    E_scale  = 2.0/(E_cut - Dat_out.Hdr['E_min'])
    E_shift  = 0.5*(E_cut + Dat_out.Hdr['E_min'])
    sf_scale, sf_shift, sf_mask = get_stp_scaler(Dat_in.Stps)
    for iStru in range(Dat_in.Hdr['nStru']):
        E      = Strs['E'][iStru]
        nAtom  = Strs['nAtom'][iStru]
        E_norm = E_scale*(E - E_shift*nAtom)
        if np.float64(E_norm/nAtom) > 1e0:
            mask_str[iStru] = False
            continue
        Strs['E'][iStru]     = E_norm
        Strs['force'][iStru] = np.multiply(E_scale, Strs['force'][iStru])
        for iAtom in range(nAtom):
            Type = Strs['Type'][iStru][iAtom]
            Strs['sfval'][iStru][iAtom] = np.multiply(sf_scale[Type][sf_mask[Type]],
                    np.add(Strs['sfval'][iStru][iAtom][sf_mask[Type]], 
                    sf_shift[Type][sf_mask[Type]]))
    
    for key,item in Strs.items():
        Dat_out.Strs[key] = item[mask_str]
    
    for iType in range(Dat_out.Hdr['nType']):
        for key in ['sf','sfParam','sfEnv','sfval_min','sfval_max','sfval_avg','sfval_cov']:
            Dat_out.Stps[iType][key] = Dat_out.Stps[iType][key][sf_mask[iType]]
        Dat_out.Stps[iType]['nsf']   = np.sum(sf_mask[iType])

    nExStr  = np.sum(~mask_str)
    nExAtom = np.sum(Strs['nAtom'][~mask_str])
    E = Strs['E'][mask_str]/Strs['nAtom'][mask_str]
    if nExStr != 0:
        print('  {} high energy Structures are excluded'.format(nExStr))
    Dat_out.Hdr['E_scale']   = np.float64(E_scale)
    Dat_out.Hdr['E_shift']   = np.float64(E_shift)
    Dat_out.Hdr['nAtom']    -= np.int32(nExAtom)
    Dat_out.Hdr['nStru']    -= np.int32(nExStr)
    Dat_out.Hdr['E_min']     = np.float64(np.min(E))
    Dat_out.Hdr['E_max']     = np.float64(np.max(E))
    Dat_out.Hdr['E_avg']     = np.float64(np.average(E))
    Dat_out.Hdr['isNorm'][:] = True

    return Dat_out

class build_nnp:
    """
build_nnp(network=None, setup=None, trainset=None)
------------------------------------------------------------------
build neural network potential (NNP) dataset that consists of 
'.fns', '.Nets', '.Stps', and '.TrnSet'
    
default: build empty NNP
    network  = None
    setup    = None
    trainset = None

network  : dict, collection of '.fns' and '.Nets'.
           corresponds to 'inp.network'
setup    : dict, collection of SETUP data.
           corresponds to 'inp.network'
trainset : dict, part of 'data.Hdr'

.fns(nType)  : list, collection of filenams. 

.Nets(nType) : dict, collection of network data.
    nLay            : number of layers
    nNodes_max      : maximum number of layer
    nW              : number of weights
    nNode           : total number of nodes
    nNodes(nLay)    : number of nodes of each layer
    ActFunc(nLay)   : type of activation function of each layer
    iW_max(nLay)    : index of last weight of each layer
    iNode_max(nLay) : index of last node of each layer
    W(nW)           : all weights

.Stps(nType) : dict, collection of SETUP data.
    Type                  : type of SETUP data
    Des                   : description 
    Atom                  : name of central atom
    nEnv                  : number of environment
    Atoms                 : name of environment atoms
    r_min                 : minimum distance btw atoms
    r_max                 : maximum distance btw atoms
    sfType                : type of symmetry function, Behler only supported
    nsf                   : number of symmetry functions
    nsfParam              : maximum number of symmetry function parameter
    sf(nsf)               : index of symmetry function 
    sfParam(nsf,nsfParam) : parameter of symmetry function 
    sfEnv(nsf,2)          : type of environment atom species 
    nEval                 : number of performed evaluation 
    sfval_min(nsf)        : mininum of symmetry function value 
    sfval_max(nsf)        : maximum of symmetry function value 
    sfval_avg(nsf)        : average of symmetry function value 
    sfval_cov(nsf)        : covariance of symmetry function value

.TrnSet      : dict, collection of TrnSet data.
    fn            : filename of training set data
    nType         : number of atom type
    nStru         : number of structure
    Atoms         : name of atom species
    E_Atom(nType) : atomic energies of atom species
    isNorm        : True if data is normalized
    E_scale       : scaleing factor of energy
    E_shift       : shift of energy
    nAtom         : number of total atom
    E_avg         : average of cohesive energy
    E_min         : minimum of cohesive energy
    E_max         : maximum of cohesive energy

    """
    def __init__(self, network=None, setup=None, trainset=None):
        import copy
        self.fns    = []
        self.Nets   = {}
        self.Stps   = {}
        self.TrnSet = {}
        if network is not None:
            self.fns    = network['fns']
            self.Nets   = network['Nets']
        if setup is not None:
            self.Stps   = setup
        if trainset is not None:
            key_ts = ['fn','nType','nStru','Atoms','E_Atom',
                'isNorm','E_scale','E_shift','nAtom','E_avg','E_min','E_max']
            for key in key_ts:
                self.TrnSet[key] = trainset[key]

def save_networks(NNP, iEpoch=None):
    """
save_networks(NNP, iIter)
------------------------------------------------------------------
write AENET NEURAL NETWORK POTENTIAL files named by '.fns' in 
fortran binary format for each atom species. 
NNP is aenet_tool.build_nnp() object that consists of '.fns', 
'.Nets', '.Stps', and 'TrnSet'.
if iIter is given, additional AENET NEURAL NETWORK POTENTIAL files
named by '.fns'+'-%5.5i'%(iIter) will be written.
    """
    for Type, fn in enumerate(NNP.fns):
        with FortranFile(fn,'w') as f:
            _write_fbytes_network(f, NNP.Nets[Type])
            _write_fbytes_setup(f, NNP.Stps[Type])
            _write_fbytes_trnset(f, NNP.TrnSet)
    if iEpoch is None:
        return
    if not pu.assert_num(iEpoch):
        pu.call_error('input error - aenet_tool.save_network',
            'given \'iIter\' is not a number - '+str(iEpoch))
    for Type, fn in enumerate(NNP.fns):
        with FortranFile(fn+'-{:06d}'.format(int(iEpoch)),'w') as f:
            _write_fbytes_network(f, NNP.Nets[Type])
            _write_fbytes_setup(f, NNP.Stps[Type])
            _write_fbytes_trnset(f, NNP.TrnSet)

def load_networks(fns):
    """
load_networks(fns)
------------------------------------------------------------------
read AENET NEURAL NETWORK POTENTIAL files written by fortran
binary format by using scipy.io.FortranFile().
the function returns aenet_tool.build_nnp() object that consists of
'.fns', '.Nets', '.Stps', and 'TrnSet'.
    """
    NNP = build_nnp()
    for fn in fns:
        pu.assert_file(fn,'aenet_tool.load_network')
    NNP.fns = fns
    nType   = len(fns)
    TrnSet  = {}
    for Type,fn in enumerate(fns):
        with FortranFile(fn,'r') as f:
            NNP.Nets[Type] = _read_fbytes_network(f)
            NNP.Stps[Type] = _read_fbytes_setup(f)
            TrnSet[Type]   = _read_fbytes_trnset(f)
    for iType1 in range(nType):
        for iType2 in range(nType):
            if iType1 == iType2:
                continue
            for key in TrnSet[iType1]:
                if TrnSet[iType1][key] != TrnSet[iType2][key]:
                    pu.call_error('TrnSet info. missmatch - aenet_tool.load_network',
                        'key - %s \ %i - %s \ %i - %s'%(key, iType1, 
                        str(TrnSet[iType1][key]), iType2, str(TrnSet[iType2][key])))
    NNP.TrnSet = TrnSet[0]
    return NNP

def update_network(NNP, iEpoch, info, values):
    iW = -2
    ib = -1
    for iType, Net in enumerate(NNP.Nets):
        for iLay in range(1,Net['nLay']):
            iW += 2
            ib += 2
            infoNNP_W = '\'Type_{:02d}/layer_{}/W:0\' shape={}'.format(
                    iType, iLay, values[iW].shape)
            infoNNP_b = '\'Type_{:02d}/layer_{}/B:0\' shape={}'.format(
                    iType, iLay, values[ib].shape)
            if infoNNP_W not in str(info[iW]):
                pu.call_error(infoNNP_W, str(info[iW]))
            if infoNNP_b not in str(info[ib]):
                pu.call_error(infoNNP_b, str(info[ib]))
            Net['W'][iLay-1][:-1] = values[iW]
            Net['W'][iLay-1][-1]  = values[ib]
        NNP.Nets[iType] = Net
    save_networks(NNP, iEpoch)
    return NNP