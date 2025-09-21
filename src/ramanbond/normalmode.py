import os
import numpy as np
import math
from numpy.typing import NDArray
from typing import Tuple
from numpy._typing import NDArray
from decimal import Decimal

from .util import find_pm_results, angle_to_rgb
from .pol import polarizability_bond, AMS
from .conversion import ANGSTROM2BOHR

class raman_bond(polarizability_bond):

    def __init__(self,outputfile,parameters:Tuple[float,float]=(1.0,3.0)):
        self.chemObj = AMS(name=outputfile)
        self.parameters = parameters
        self.chemObj._collect(abort =False)
        self.atoms = self.fragment_atoms
        self.coordinates = self.fragment_coordinates
        self.__generate_dis_matrix()

    def collect_raman_derivatives(self, dir=None, sR=0.01):
        if dir is None:
            dir = os.curdir
        pnames, mnames = find_pm_results(dir)
        self.chemObj.collect_raman_derivatives(dir = dir)
        raman_atoms = np.zeros((self.nmodes, self.natoms,3, 3), dtype=np.complex128)
        raman_bonds = np.zeros((self.nmodes, self.natoms,self.natoms, 3,3), dtype=np.complex128)
        sQ = self.step_size(sR)
        for i in range(self.nmodes):
            p = polarizability_bond(pnames[i],parameters=self.parameters)
            m = polarizability_bond(mnames[i],parameters=self.parameters)
            # Atomic polarizability tensor (natoms, 3, 3)
            temp_atoms = p.atomic_polarizability - m.atomic_polarizability
            temp_atoms[:]  = temp_atoms[:]/(2 * sQ[i])
            raman_atoms[i] = temp_atoms
            # Inter atomic polarizability tensor (natoms, natoms 3, 3)
            if "VELOCITY" not in p.calctype:
                temp_bonds = (p.bond_polarizability - m.bond_polarizability)
                temp_bonds = temp_bonds/(2 * sQ[i])
                raman_bonds[i] = temp_bonds

        self.raman_atoms = raman_atoms
        self.raman_bonds = raman_bonds

        # def __generate_cos_matrix(self,dis_matrix) -> NDArray[np.float64]:
        #     cos_matrix = np.zeros((3,self.natoms, self.natoms))
        #     for i in range(self.natoms):
        #         for j in range(i):
        #             for ii in range(3):
        #                 cos_matrix[ii,i,j] = abs(np.dot((self.coordinates[i] -
        #                     self.coordinates[j]), np.identity(3)[ii]) /
        #                     dis_matrix[i,j])
        #                 cos_matrix[ii,j,i] = cos_matrix[ii,i,j]
        #     return cos_matrix

        # def __solve_lagrange(self,dis_matrix, cos_matrix,parameters:Tuple[float,float]=(1.0,3.0)):
        #     lambda_matrix = np.zeros((self.natoms, self.natoms),dtype = object)
        #     for i in range(self.natoms):
        #         for j in range(i):
        #             # NOTE: I don't know why we are using Decimal here.
        #             lambda_matrix[i,j] = Decimal(inter_atomic_lambda(parameters,self.atoms[i],
        #                                                              self.atoms[j], 
        #                                                              dis_matrix[i,j],
        #                                                              cos_matrix[i,j]))
        #             lambda_matrix[j,i] = lambda_matrix[i,j]
        #         lambda_matrix[i,i] = Decimal(-1*sum(lambda_matrix[i,:])) # Add an arbitrary constant C = 1.0
        #     lambda_matrix += Decimal('1.0')
        #
        #     lambda_solution = np.linalg.solve(lambda_matrix, self.hirshfeld_induced_charges.T)
        #
        #     charge_flow = np.zeros((self.natoms, self.natoms), dtype = object)
        #
        #     for i in range(self.natoms):
        #         for j in range(i):
        #             charge_flow[i, j] = calc_charge_flow(parameters, self.atoms[i], self.atoms[j],
        #                                                  dis_matrix[i,j], cos_matrix[i,j],
                    #                                      lambda_solution[i], lambda_solution[j])
                    # charge_flow[j,i] = charge_flow[i,j]
                    #
            # return
except ImportError:
    import numpy as np
    class normalmode(object):
        '''
        A mini-chemPackage. Targeted to collect information from AMS2021 above results.
        Information to collect:
        - Geometry
        - Normal modes
        '''
        def __init__(self, freqout):
            self.filename = 'freqout'
            f = open(freqout)
            fl = f.readlines()
            self.natoms = 0
            self.atoms = None
            self.coordinates = None
            self.collect_geometry(fl)
            f.close()
        # Only to be overide across different engine I guess
        def collect_geometry(self, fl):
            # self.atomnote = 
            # self.coord = 
            s=e=0
            for i in range(len(fl)):
                if "Geometry" in fl[i]:
                    ii=0
                    i+=5
                    s=i
                    while True:
                        ii+=1
                        if len(fl[i].split()) == 0:
                            e = i
                            break
                        else:
                            i+=1

                break
            self.natoms = e-s
            self.coordinates = np.empty((self.natoms,3), dtype=float)
            self.atoms = np.empty((self.natoms),dtype='<U1')
            for i, line in enumerate(fl[s:e]):
                self.coordinates[i] = line.strip("\n").split()[1:]
                self.atoms[i] = line.strip("\n").split()[0]
            return

        # This is different for ADF and AMS
        def vib_freq(self,line):
            ln = line.split()
            if not ln:
                return False
            elif line == " Index  Atom      ---- Displacements (x/y/z) ----":
                return True
            else:
                return False

        # TODO: Implementation
        def collect_normalmode(self, fl):
            return

        def calc_polbond(self, p_files, m_files, component = "all"):
# TODO: Check if there is a one-to-one correspondence between plus and minus
            from .pol import pol

            return

