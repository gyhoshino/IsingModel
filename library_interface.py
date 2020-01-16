from lattice_pylib import IsingLattice as LatticeLib_Python
import numpy as np
import sys
from ctypes import *
from os.path import isfile
import time
'''File that returns an object IsingLattice which either
interfaces with the C++ library (if available) or a pure python 
version (if not).'''

# Interface to use the C++ library
def IsingLattice(inp, set_seed = 0):
    if set_seed == 0:
        seed = int(time.time())
    elif set_seed < 0:
        seed = abs(set_seed)
    elif set_seed > 0:
        seed = int(time.time())+set_seed
        
    '''return an IsingLattice object, use the cpp object if desired and possible'''
    if not inp['use_cpp']:
        return LatticeLib_Python(inp['N'], inp['r_flip'])
    elif not isfile('./lattice_lib.so'):
        print('Warning: shared library "lattice_lib.so" is not present.\n'
              'Therefore the c++ library cannot be used. Reverting to\n'
              'the python library')
        return LatticeLib_Python(inp['N'], inp['r_flip'])

    # else return the c++ library version
    libc = CDLL("./lattice_lib.so")
    libc.newMatrix.restype  = c_void_p
    libc.newMatrix.argtypes = [c_int, c_int, c_int]

    libc.step.argtypes        = [ c_void_p, c_float, c_float ]
    libc.wolff.argtypes       = [ c_void_p, c_float ]
    libc.nsteps.argtypes      = [ c_void_p, c_int, c_float, c_float ]

    libc.get_site_flip_prob.argtypes   = [ c_void_p, c_int, c_int, c_float, c_float ]
    libc.get_site_flip_prob.restype    = c_float

    libc.get_E.restype = c_float
    libc.get_M.restype = c_float
    libc.n_spincorr.restype = c_int
    libc.auto_correlation.restype = c_float

    libc.spincorr.restype = c_float
    libc.auto_correlation.restype = c_float

    class LatticeLib_cpp(object):
        def __init__(self, N, flip_prop, seed):
            n_flip = N*N*flip_prop
            self.N = N
            if n_flip != int(n_flip):
                n_flip = int(n_flip) + 1
            self.matrix_ptr = c_void_p(libc.newMatrix(int(N),int(n_flip), seed))
            self.auto_correlation = []
        def get_site_flip_prob(self, i, j, T, B):
            return libc.get_site_flip_prob( self.matrix_ptr, i, j, T, B)
        def free_memory(self):
            # print ('deleting ising_lattice_lib object')
            libc.delMatrix( self.matrix_ptr ) 
        def step(self, T, B):
            return libc.step(self.matrix_ptr, c_float(T), c_float(B))
        def wolff(self, T):
            return libc.wolff(self.matrix_ptr, c_float(T))
        def nsteps(self, T, B, n):
            return libc.nsteps(self.matrix_ptr,n,T,B)
        def get_E(self):
            return float( libc.get_E(self.matrix_ptr) )
        def get_M(self):
            return float( libc.get_M(self.matrix_ptr) )

        def get_SC_v0(self):
            # calculate the spin correlation
            # by the legacy definition
            libc.calc_auto_correlation(self.matrix_ptr)
            self.auto_correlation = []
            for i in range(1,int(self.N/2)):
                val = libc.auto_correlation(self.matrix_ptr, i)
                self.auto_correlation.append(val)
            return self.auto_correlation

        def get_SC_v1(self):
            # calculate correlation by
            # <a*b>-<a><b>
            libc.calc_spincorr(self.matrix_ptr)
            self.spincorr = []
            for i in range(1,int(self.N/2)):
                val = libc.spincorr(self.matrix_ptr, i)
                self.spincorr.append(val)
            return self.spincorr

        def set_Nflip(self,npick):
            return  libc.set_Npick(self.matrix_ptr, c_int(npick))
        def set_flip_prop(self,flip_prop):
            n_flip = self.N*flip_prop
            if n_flip != int(n_flip):
                n_flip = int(n_flip) + 1
            else:
                n_flip = int(n_flip)
            return (libc.set_Npick(self.matrix_ptr, c_int(n_flip)))
        def get_Nspin(self):
            return libc.get_Nspin(self.matrix_ptr)
        def get_Nalign(self):
            return libc.get_Nbond(self.matrix_ptr)
        def get_spin(self, i, j):
            return libc.get_spin(self.matrix_ptr, i, j)
        def set_spin(self, i, j, k):
            return libc.set_spin(self.matrix_ptr, i, j, k)
        def print_spins(self):
            return libc.print_spins(self.matrix_ptr)
        def print_aligned(self):
            return libc.print_bonds(self.matrix_ptr)
        def randomize_spins(self):
            return libc.initialize_spins(self.matrix_ptr)
        def get_numpy_spin_matrix(self):
            matrix = np.ones((self.N, self.N))
            for i in range(self.N):
                for j in range(self.N):
                    matrix[i,j] = libc.get_spin(self.matrix_ptr,i,j)
            return matrix
    return LatticeLib_cpp(inp['N'],inp['r_flip'],seed)

'''Some simple tests to see that it worked.'''

if __name__ == '__main__':
    print('first')
    from sys import argv
    n = 10
    inp = {'N':50,'r_flip':0.10,'use_cpp':True}
    if len(argv) == 1:
        pass
    elif argv[1] == "1":
        lattice = IsingLattice(inp)
        lattice.print_spins()
        print("----")
        for i in range(n):
            for j in range(n):
                continue
                lattice.set_spin(i,j,1)
        lattice.print_spins()
        print("Nspin : %i"%lattice.get_Nspin())
        print(" This is mag: %f"%lattice.get_M())
        print("auto_corr:")
        x = lattice.calc_auto_correlation()
        for val in x:
            print(val)
        lattice.free_memory()
        
    elif argv[1] == "0":
        inp['N'] = 10
        lattice = IsingLattice(inp)
        print(argv[0])
        lattice.print_spins()
        print('---')
        lattice.set_spin(0,0,1)
        lattice.print_spins()
        print('---')
        lattice.set_spin(0,0,-1)
        lattice.print_spins()
        print('---')

        flip = True
        for i in range(n):
            for j in range(n):
                flip = not flip
                if flip:
                    lattice.set_spin(i,j,1)
                else:
                    lattice.set_spin(i,j,-1)
        lattice.print_spins()
        print('---')

        lattice.free_memory()
    else:
        inp['N']=10
        inp['r_flip'] = 0.01
        lattice = IsingLattice(inp)
        lattice.print_spins()
        lattice.print_aligned()
        lattice.nsteps(4.,2.1,10000)

        # lattice.randomize_spins()
        print('E: ',lattice.get_E())
        print('M: ',lattice.get_M())
        print('autocorrelation ', lattice.calc_auto_correlation())
        print("set_Nflip:      ", lattice.set_Nflip(20))
        print("get_Nspin:      " , lattice.get_Nspin())
        print("get_Nalign:     " , lattice.get_Nalign())
        # print("step:           " , lattice.step(2.9,1.0))
        print("get_Nspin:      " , lattice.get_Nspin())
        print("get_Nalign:     " , lattice.get_Nalign())
        print('E: '              , lattice.get_E())
        print("N               " , lattice.N)
        print('aligned')
        # lattice.print_aligned()
        print('spins')
        # lattice.print_spins()
        lattice.free_memory()
    print('exit')
    sys.exit()
