"""
lp34vasp
------------------------------------------------------------------
This module interfaces VASP INPUTS/OUTPUTS with Python3 for
easier data analysis

usage : 

    import pyutil.lp34vasp as MODULENAME
    from pyutil.lp34vasp import FUNCTION

functions :
    .get_outcar(VASP_OUTCAR)
    .get_doscar(VASP_DOSCAR, input)
    .get_eigenval(VASP_EIGENVAL)
    .coo()
    .get_poscar(VASP_POSCAR)
    .write_poscar(filename, dataset)
    .get_xsf(XCRYSDEN_XSF)
    .write_xsf(filename, dataset)

"""

import os
import numpy as np
from .util import *

class get_outcar:
    """ 
get_outcar(f_name)
------------------------------------------------------------------
read VASP OUTCAR and return dataset as python object

outputs :
    .f_name  : filename of OUTCAR
    .input   : INCAR variables
    .loop    : total elapsed loop
    .is_band : 'True' if band structure calculation
    .is_stat : 'True' if calculation is static self-consistant
    .is_end  : 'True' if calculation finished

internal functions :
    .get_energy()    : VASP total energies
    .get_lattice()   : real and reciprocal lattice vectors
    .get_magnetize() : magnetization vectors for each atom
            (LORBIT >= 10 required)
    .get_total_magnetize() : total magnetization
    .get_force()     : force vectors for each atom
 
    """
    def __init__(self, f_name):
        self.f_name  = f_name
        assert_file(self.f_name,'get_outcar')
        self.input   = self._get_inputs()
        self.energy  = self._get_energy()
        self.is_band = False
        self.is_stat = False
        if get_line('interpolating k-points', self.f_name, iLine=1) is not None: 
            self.is_band = True
        if not self.is_band and self.input['nsw'] == 0: 
            self.is_stat = True
        self.loop   = 0 
        self.is_end = self._get_ended()

    def _get_inputs(self):
        input = {}
        input['nsw']     = int(chk_len(get_line('  NSW', self.f_name, 1), 2))
        input['nedos']   = int(chk_len(get_line('  NEDOS', self.f_name, 1), 5))
        input['nions']   = int(chk_len(get_line('  NIONS', self.f_name, 1), 11))
        input['ispin']   = int(chk_len(get_line('  ISPIN', self.f_name, 1), 2))
        input['lorbit']  = int(chk_len(get_line('  LORBIT', self.f_name, 1), 2))
        input['lwave']   = chk_len(get_line('  LWAVE', self.f_name, 1), 2)
        input['lcharg']  = chk_len(get_line('  LCHARG', self.f_name, 1), 2)
        input['lvtot']   = chk_len(get_line('  LVTOT', self.f_name, 1), 2)
        input['lsorbit'] = chk_len(get_line('  LSORBIT', self.f_name, 1), 2)
        input['ncollin'] = chk_len(get_line('  LNONCOLLINEAR', self.f_name, 1), 2)
        input['nkpts']   = int(chk_len(get_line('  NKPTS', self.f_name, 1), 3))
        input['nbands']  = int(chk_len(get_line('  NBANDS=', self.f_name, 1), 14))
        input['nelect']  = float(chk_len(get_line('  NELECT =', self.f_name, 1), 2))
        input['atoms']   = self._get_atoms()
        input['nspec']   = [int(item) for item in get_line('ions per type =', self.f_name, 1)[4:]]
        return input

    def _get_energy(self):
        enes = {}
        enes['E_fermi']  = np.float64(chk_len(get_line('E-fermi', self.f_name), 2))
        enes['E_free']   = np.float64(chk_len(get_line('free  energy', self.f_name), 4))
        line = get_line('energy  without', self.f_name)
        enes['E_w/oS']   = np.float64(chk_len(line, 3))
        enes['E_sigma0'] = np.float64(chk_len(line, 6))
        return enes

    def _get_ended(self):
        isend = False
        self.loop = int(os.popen('grep \'EDIFF is reached\' %s | wc -l'%(self.f_name)).read())
        if self.is_stat:
            if self.loop == 1: 
                isend = True
        else:
            line = get_line('reached required accuracy', self.f_name, 1)
            if line is not None: 
                isend = True
        return isend

    def _get_atoms(self):
        lines = get_line('VRHFIN', self.f_name, iLine='all')
        Atoms = []
        for line in lines:
            Atoms.append(line[1].replace(':','').replace('=',''))
        return Atoms

    def get_lattice(self):
        """
.get_lattice()
------------------------------------------------------------------
return (2,3,3) optimized real and reciprocal lattice vectors

    get_lattice()[0] : (3,3) real lattice vector
    get_lattice()[1] : (3,3) reciprocal lattice vector

        """
        lines = get_lines('direct lattice vectors',self.f_name,3)
        if lines is None: 
            call_error('reading error - get_lattice',
                'lattice information is missing')
        lvec = np.zeros([2,3,3],dtype=np.float64)
        for i,line in enumerate(lines):
            lvec[0][i] = np.array(line[0:3],dtype=np.float64)
            lvec[1][i] = np.array(line[3:6],dtype=np.float64)
        return lvec

    def get_magnetize(self):
        """
.get_magnetize()
------------------------------------------------------------------
return (nAtom,3) atomic magnetization if spin-polarized or non-collinear
calculation is performed with LORBIT >= 10

    nAtom : number of total atoms in POSCAR order

    - collinear calculation with ISPIN = 2 -
        get_magnetize()[0:nAtom,2]   : + is spin up and - is spin down

    - non collinear calculation
        get_magnetize()[0:nAtom,0:2] : magnetization vector

        """
        iscollin = True
        mag = np.zeros([self.input['nions'], 3],dtype=np.float64)
        nline = mag.shape[0]+3
        if self.input['ncollin'] == 'T': 
            nline = 3*nline + 12
            iscollin = False
        if not (self.input['ispin'] == 2 or (not iscollin)): 
            call_error('input error - get_magnetize','spin polarization required')
        if self.input['lorbit'] < 10: 
            call_error('input error - get_magnetize','LORBIT >= 10 required')
        lines = get_lines('magnetization (x)',self.f_name,nline,0)
        if lines is None: 
            call_error('reading error - get_magnetize','polarization is missing')
        for iAtom in range(self.input['nions']):
            if iscollin:
                mag[iAtom][2] = np.float64(lines[iAtom+3][4])
            else:
                mag[iAtom] = np.float64([lines[iblock*(self.input['nions']+
                             9)+iAtom+3][4] for iblock in range(3)])
        return mag

    def get_total_magnetize(self):
        """
.get_total_magnetize()
------------------------------------------------------------------
return total magnetization (or vector) if spin-polarized or non-collinear
calculation is performed

    collinear with ISPIN = 2 : total magnetization 
    non-collinear            : total magnetization vector

        """
        iscollin = True
        if self.input['ncollin'] == 'T': iscollin = False
        if not (self.input['ispin'] == 2 or (not iscollin)):
            call_error('input error - get_total_magnetize',
                'spin polarization required')
        mag = np.zeros((3), dtype=np.float64)
        line = get_line(' number of electron ',self.f_name,n=0)
        if line is None: 
            call_error('reading error - get_total_magnetize',
                'magnetization is missing')
        else:
            line = line.split()
        if iscollin:
            mag[2] = np.float64(line[5])
        else:
            mag = np.float64(line[5:8])
        return mag
 
    def get_force(self):
        """
.get_force()
------------------------------------------------------------------
return (n,3) atomic force 

        """
        force = np.zeros([self.input['nions'], 3], dtype=np.float64)
        lines = get_lines('TOTAL-FORCE (eV/Angst)',self.f_name,force.shape[0]+1,0)
        if lines is None: 
            call_error('reading error - get_force','force is missing')
        for iAtom,line in enumerate(lines[1:]):
            force[iAtom] = np.float64(line[3:6])
        return force

class get_doscar:
    """
get_doscar(f_name, input)
------------------------------------------------------------------
read VASP DOSCAR named 'f_name' and return dataset.

input : INCAR tags get from lp.get_outcar().input

outputs :
    .tdos     : total DOS
    .pdos     : lm and site projected DOS (LORBIT != 0 required)
    .nedos    : number of grid
    .nAtom    : number of atoms, required by PDOS
    .ef       : fermi level
    .header_T : header of total DOS
    .header_P : header of PDOS

function :
    .write_dos(filename, dataset, header) :
            write dos datafile named by 'filename' by using DOS 
            data 'dataset' and column names 'header'.

    """
    def __init__(self, f_name, inp):
        self.f_name = f_name
        assert_file(self.f_name,'get_doscar')
        self.ispin    = inp['ispin']
        self.iscollin = True
        if inp['ncollin'] == 'T':
            self.iscollin = False
        self.lorbit   = inp['lorbit']

        with open(self.f_name,'r') as f:
            self.nAtom = int(chk_len(f.readline().split(), 0))
            for i in range(4): 
                f.readline() # read buffer lines
            line = f.readline().split() # data line
            self.nedos = int(chk_len(line, 2)) - 1
            self.ef    = np.float64(chk_len(line, 3))
            if self.nedos+1 != inp['nedos']:
                call_error('input error - get_doscar',
                    'NEDOS does not match - %i %i'%(inp['nedos'],self.nedos))
            if self.nAtom != inp['nions']:
                call_error('input error - get_doscar',
                    'NIONS does not match - %i %i'%(inp['nions'],self.nAtom))
            self.tdos = self._get_tdos(f)
            self.pdos = self._get_pdos(f)

    def _get_tdos(self,fo):
        # TDOS header
        self.header_T = ['E-EF','Tot']
        if self.ispin == 2:
            self.header_T.append(self.header_T[1]+'_D')
            self.header_T[1] += '_U'

        tdos = np.zeros((self.nedos,self.ispin+1),dtype=np.float64)

        fo.readline() # first line, remove
        for i in range(self.nedos):
            line = self.f.readline().split()
            try:
                tdos[i][0]  = np.float64(line[0]) - self.ef
                tdos[i][1:] = np.float64([line[j+1] for j in range(self.ispin)])
            except ValueError:
                call_error('reading error - get_doscar._get_tdos',
                    'TDOS value error line %i'%(i))
            except IndexError:
                call_error('reading error - get_doscar._get_tdos',
                    'TDOS indexing error line %i - %i'%(i,len(line)))
        return tdos

    def _get_pdos(self, fo):
        # PDOS header
        if self.lorbit > 10:
            header_tmp = ['s','px','py','pz','dxy','dyz','dz2','dxz','dx2y2']
        elif self.lorbit != 0:
            header_tmp = ['s','p','d']
        else:
            return None
        self.header_P = ['E-EF']
        if self.iscollin and self.ispin == 1:
            self.header_P += header_tmp
        elif self.ispin == 2:
            for item in header_tmp:
                self.header_P.append(item+'_U')
                self.header_P.append(item+'_D')
        else:
            for item in header_tmp:
                self.header_P.append(item+'-Tot')
                self.header_P.append(item+'-mx')
                self.header_P.append(item+'-my')
                self.header_P.append(item+'-mz')
        self.ncol = len(self.header_P)
        pdos = np.zeros([self.nAtom,self.nedos,self.ncol],dtype=np.float64)
        for iAtom in range(self.nAtom):
            fo.readline() # data line
            fo.readline() # first line, remove
            for idos in range(self.nedos):
                line = fo.readline().split()
                pdos[iAtom,idos]    = np.float64(line)
                pdos[iAtom,idos,0] -= self.ef
        return pdos

    def write_dos(self, f_name, dos, header):
        """
.write_dos(f_name, dos, header)
------------------------------------------------------------------
write dos datafile with filename 'f_name' using data 'dos' and 
data description header 'header'

        """
        with open(f_name,'w') as f:
            f.write('#%9s'%(header[0]))
            for str in header[1:]: 
                f.write(' %10s'%(str))
            f.write('\n')
            for line in dos:
                for dat in line:
                    f.write('%10.4f '%(dat))
                f.write('\n')

class get_eigenval:
    """
get_eigenval(f_name)
------------------------------------------------------------------
read VASP EIGENVAL named by 'f_name' and return dataset.

outputs :
    .f_name : filename of EIGENVAL
    .ispin  : ISPIN
    .nelect : number of electrons
    .nkpts  : number of kpoints
    .nbands : number of bands
    .kpts   : kpoint coordination in xyz (.nkpts, 3)
    .weight : weight of kpoint (.nkpts)
    .eigens : eigenvalues (.nkpts, .nbands, .ispin)

    """
    def __init__(self, f_name, inp):
        self.f_name = f_name
        assert_file(self.f_name,'get_eigenval')
        with open(self.f_name,'r') as f:
            # get information
            line = f.readline()
            try:
                self.ispin = int(line.split()[3])
            except ValueError:
                call_error('reading error - get_eigenval',
                    'ISPIN value - %s'%(line))
            except IndexError:
                call_error('reading error - get_eigenval',
                    'ISPIN index - %i'%(len(line.split())))
            for i in range(4):
                f.readline()
            line_orig = f.readline()
            line = line_orig.split()
            try:
                self.nelect = int(line[0])
                self.nkpts  = int(line[1])
                self.nbands = int(line[2])
            except ValueError:
                call_error('reading error - get_eigenval',
                    'NELECT, NKPTS, NBANDS value - %s'%(line_orig))
            except IndexError:
                call_error('reading error - get_eigenval',
                    'NELECT, NKPTS, NBANDS index - %i'%(len(line)))
            self.eigens = np.zeros([self.nkpts,self.nbands,self.ispin],dtype=np.float64)
            self.kpts   = np.zeros([self.nkpts,3], dtype=np.float64)
            self.weight = np.zeros([self.nkpts], dtype=np.float64)
 
            if inp['nelect'] != self.nelect:
                call_error('input error - get_eigenval',
                    'NELECT does not match - %i %i'%(inp['nelect']), self.nelect)
            if inp['ispin'] != self.ispin:
                call_error('input error - get_eigenval',
                    'ISPIN does not match - %i %i'%(inp['ispin']), self.ispin)
             # get eigenvalues  
            for j in range(self.nkpts):
                f.readline() # empty line
                line_orig = f.readline()
                line = line_orig.split()
                try:
                    self.kpts[j]   = np.float64(line[:3])
                    self.weight[j] = np.float64(line[3])
                except ValueError:
                    call_error('reading error - get_eigenval',
                        'KPOINTS value - %s'%(line_orig))
                except IndexError:
                    call_error('reading error - get_eigenval',
                        'KPOINTS index - %i'%(len(line)))
                for i in range(self.nbands):
                    line_orig = f.readline()
                    line = line_orig.split()
                    try:
                        self.eigens[j][i] = np.float64(line[1:1+self.ispin])
                    except ValueError:
                        call_error('reading error - get_eigenval',
                            'EIGENS value - %s'%(line_orig))
                    except IndexError:
                        call_error('reading error - get_eigenval',
                            'EIGENS index - %i'%(len(line)))

class build_coo:
    """
coo
------------------------------------------------------------------
build empty coordinate dataset for .write_poscar or .write_xsf

    """
    def __init__(self):
        self.f_name  = None
        self.des     = None
        self.lvec    = None
        self.ispbc   = None
        self.nAtom   = None
        self.Atoms   = None
        self.Atoms_l = None
        self.nSpec   = None
        self.iscart  = None
        self.coord   = None
        self.vector  = None

class get_poscar(build_coo):
    """
get_poscar(f_name)
------------------------------------------------------------------
read VASP POCSCAR 'f_name' and return dataset.

outputs :
    .f_name : filename of POSCAR
    .des    : description of POSCAR at the first line
    .lvec   : lattice vector
    .ispbc  : periodic boundary condition, True
    .style  : POSCAR style, 4 or 5
    .nAtom  : total number of atoms
    .Atoms   : atom species (None if .style = 4, (.nSpec))
    .Atoms_l : atom species (None if .style = 4, (.nAtom))
    .nSpec  : number of each atom species
    .isselc : True if selective dynamics used
    .iscart : True if coordinate is cartesian
    .coord  : atomic coordination in xyz (.nAtom, 3)
    .vector : vector in xyz (.nAtom, 3)
    .selc   : if isselc is True, return selective dynamics flags
            (.nAtoms, 3)
    """
    def __init__(self, f_name):
        super().__init__()
        self.f_name = f_name
        assert_file(self.f_name, 'get_poscar')
        with open(self.f_name) as f:
            self.des = f.readline()
            # get lattice constant
            latcont = np.float64(f.readline().split()[0])
            # get lattice vector
            self.lvec = np.float64([f.readline().split()[:3] 
                for j in range(3)])
            self.ispbc = True
            # get atom types and numbers
            line = f.readline().split()
            if assert_num(line[0]):
                self.style = 4
            else:
                self.style = 5
                self.Atoms  = line
                line = f.readline().split()
            self.nSpec = np.int64(line)
            self.nAtom = np.sum(self.nSpec)
            if self.Atoms is not None:
                self.Atoms_l = []
                for i,j in enumerate(self.nSpec):
                    self.Atoms_l += [self.Atoms[i] for k in range(j)]
            self.isselc = False
            self.iscart = False
            # get POSCAR type (selective dynamics // fractional // cartesian)
            line = f.readline()
            if line[0] == 'S' or line[0] == 's':
                self.isselc = True
                line = f.readline()
            if line[0] == 'C' or line[0] == 'c':
                self.iscart = True
            elif line[0] == 'D' or line[0] == 'd':
                self.iscart = False
            else:
                call_error('reading error - get_poscar','coord type - '+line)
            # get coordinate
            self.coord = np.zeros((self.nAtom,3), dtype=np.float64)
            self.vector = np.zeros((self.nAtom,3), dtype=np.float64)
            if self.isselc: 
                self.selc = [['T' for i in range(3)] for j in range(self.nAtom)]
            for iAtom in range(self.nAtom):
                line_orig = f.readline()
                line = line_orig.split()
                try:
                    self.coord[iAtom] = np.float64(line[:3])
                except ValueError:
                    call_error('reading error - get_poscar',
                        'COO value - %s'%(line_orig))
                except IndexError:
                    call_error('reading error - get_poscar',
                        'COO index - %i'%(len(line)))
                if self.isselc:
                    try:
                        self.selc[iAtom] = line[3:6]
                    except IndexError:
                        call_error('reading error - get_poscar',
                            'SELC index - %i'%(len(line)))
            # get vectors if exists
            f.readline()
            for iAtom in range(self.nAtom):
                line = f.readline().split()
                try:
                    self.vector[iAtom] = np.float64(line[:3])
                except (EOFError, ValueError):
                    break
                except IndexError:
                    continue

class write_poscar:
    """
write_poscar(dat, islong=True)
------------------------------------------------------------------
get coordinate dataset 'dat' and write VAPS POSCAR.

dataset inputs :
    .f_name : filename 
    .des    : description of file at first line
    .lvec   : lattice vector
    .nAtom  : total number of atoms
    .Atoms  : atom species
    .nSpec  : number of each atom species
    .isselc : whether write selective dynamics or not
    .iscart : if dat.coord is fractional, False
    .coord  : atomic coordination in xyz (nAtom, 3)
    .vector : vector in xys (nAtom, 3) 
    .selc   : selective dynamics flags

    """
    def __init__(self, dat, islong=True):
        if dat.f_name is None:
            call_error('input error - write_poscar',
                'dataset is empty - f_name')
        if not (('POSCAR' in dat.f_name) or ('CONTCAR' in dat.f_name)):
            dat.f_name = 'CONTCAR'+dat.f_name
        with open(dat.f_name,'w') as f:
            if dat.isselc and (np.shape(dat.selc) != np.shape(dat.coord)):
                call_error('input error - write_poscar',
                    'dimension does not match btw COO and SELC - %s %s'%(
 	                np.shape(dat.coord), np.shape(dat.selc),))
            # write description
            f.write(dat.des)
            # write lattice inform.
            f.write('1.000000000\n')
            for i in dat.lvec: 
                f.write(' ')
                if islong: 
                    f.write(''.join(['%22.16f'%(j) for j in i])+'\n')
                else: 
                    f.write(''.join(['%12.6f'%(j) for j in i])+'\n')
             # write atom info 
            f.write(''.join(['%4s'%(i) for i in dat.Atoms])+'\n')
            f.write(''.join(['%4i'%(i) for i in dat.nSpec])+'\n')
            if hasatter(dat,'isselc'):
                if dat.isselc: 
                    f.write('Selective Dynamics\n')
            else:
                dat.isselc = False
            f.write('Direct\n')
            # write coordinate
            if dat.iscart:
                dat.coord = cart_to_dir(dat.lvec, dat.coord)
            for idx,coo in enumerate(dat.coord):
                if islong: 
                    f.write(''.join(['%20.16f'%(j) for j in coo]))
                else: 
                    f.write(''.join(['%10.6f'%(j) for j in coo]))
                if dat.isselc: 
                    f.write(''.join(['%3s'%(j) for j in dat.selc[idx]]))
                f.write('\n')
            # write vector
            if dat.vector is not None:
                f.write('\n')
                for vec in dat.vector:
                    if islong: 
                        f.write(''.join(['%20.16f'%(j) for j in vec])+'\n')
                    else: 
                        f.write(''.join(['%10.6f'%(j) for j in vec])+'\n')

class get_xsf(build_coo):
    """
get_xsf(f_name, dat)
------------------------------------------------------------------
read XCrysDen XSF 'f_name' and return dataset 'dat'

outputs :
    .f_name  : filename of XSF
    .des     : description of XSF at the first line
    .lvec    : lattice vector
    .iscart  : cartesian coordinate, True
    .ispbc   : periodic boundary condition
    .nAtom   : total number of atoms
    .Atoms   : atom species (len(.nSpec))
    .Atoms_l : atom species (.nAtom)
    .nSpec   : number of each atom species
    .coord   : atomic coordination in xyz (.nAtom, 3)
    .vector  : vector, drawn in VESTA (.nAtom, 3)

    """
    def __init__(self, f_name):
        super().__init__()
        assert_file(f_name,'get_xsf')
        self.f_name = f_name
        self.des    = ''
        self.iscart = True
        self.ispbc  = True
        self.nAtom  = 0
        self.Atoms  = []
        self.nSpec  = [0]
        with open(self.f_name,'r') as f:
            self.des   = f.readline()
            self.ispbc = False
            for line in f:
                if 'CRYSTAL' in line:
                    self.ispbc = True
                if 'PRIMVEC' in line or 'CONVVEC' in line:
                    self.lvec  = np.float64([f.readline().split() for i in range(3)])
                if 'PRIMCOORD' in line:
                    self.nAtom = int(f.readline().split()[0])
                    coord   = np.zeros((self.nAtom,3),dtype=np.float64)
                    vector  = np.zeros((self.nAtom,3),dtype=np.float64)
                    Atoms_l = ['' for i in range(self.nAtom)]
                    break
            if self.nAtom == 0:
                call_error('reading error - get_xsf','PRIMCOORD not found')
            for i in range(self.nAtom):
                line_orig = f.readline()
                line = line_orig.split()
                try:
                    Atoms_l[i] = line[0]
                    coord[i] = np.array(line[1:4],dtype=np.float64)
                except ValueError:
                    call_error('reading error - get_xsf',
                        'COO value - %s'%(line_orig))
                except IndexError:
                    call_error('reading error - get_xsf',
                        'COO index - %i'%(len(line)))
                try:
                    vector[i] = np.array(line.split()[4:7],dtype=np.float64)
                except:
                    continue
                
        if self.lvec is None:
            call_error('reading error - get_xsf','LATTICE VECTORS not found')
        self.Atoms.append(Atoms_l[0])
        for atom_i in Atoms_l:
            isdup = False
            for atom_j in self.Atoms:
                if atom_i == atom_j: 
                    isdup = True
                    break
            if not isdup:
                self.Atoms.append(atom_i)
                self.nSpec.append(0)
        i = 0
        self.Atoms_l = ['' for j in range(self.nAtom)]
        self.coord   = np.zeros((self.nAtom,3),dtype=np.float64)
        self.vector  = np.zeros((self.nAtom,3),dtype=np.float64)
        for ispec,atom_i in enumerate(self.Atoms):
            for j,atom_j in enumerate(Atoms_l):
                if atom_i != atom_j:
                    continue
                self.nSpec[ispec] += 1
                self.Atoms_l[i] = Atoms_l[j]
                self.coord[i]   = coord[j]
                self.vector[i]  = vector[j]
                i += 1

class write_xsf:
    """
write_xsf(dat, iimg=1)
------------------------------------------------------------------
get coordinate dataset 'dat' and write XCrysDen .xsf file

dataset inputs :
    .f_name  : filename
    .des     : description of XSF at first line
    .lvec    : lattice vector
    .ispbc   : periodic boundary or not
    .nAtom   : total number of Atom
    .Atoms_l : atom species
    .nspec   : number of each atom species
    .iscart  : if dat.coord is fractional, False
    .coord   : atomic coordination in xyz (.nAtom, 3)
    .vector  : vector, drawn in VESTA (.nAtom, 3)

    """
    def __init__(self, dat, iimg=1):
        if not '.xsf' in dat.f_name: 
            dat.f_name += '.xsf'
        with open(dat.f_name,'w') as f:
            f.write('# %s\n'%(dat.des))
            if (dat.ispbc):
                f.write('CRYSTAL\n')
                f.write('PRIMVEC\n')
                for vec in dat.lvec: 
                    f.write(' '+''.join(['%22.16f'%(j) for j in vec])+'\n')
                f.write('PRIMCOORD\n')
                f.write('%i %i\n'%(dat.nAtom, iimg))
            else:
                f.write('ATOMS\n')
            if not dat.iscart: 
                dat.coord = dir_to_cart(dat.lvec, dat.coord)
            if dat.Atoms_l is None:
                dat.Atoms_l = [dat.Atoms[i] for k in range(j) 
                    for i,j in enumerate(dat.nSpec)]
            for i in range(dat.nAtom):
                f.write('%2s'%(dat.Atoms_l[i]))
                f.write(''.join([' %23.16f'%(j) for j in dat.coord[i]]))
                f.write(''.join([' %13.6f'%(j) for j in dat.vector[i]]))
                f.write('\n')

