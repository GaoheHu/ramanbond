import numpy as np
from numpy._typing import NDArray
from typing import Tuple
from decimal import Decimal

from .util import inter_atomic_lambda, calc_charge_flow
from .conversion import ANGSTROM2BOHR
# This is basically a minimal version of chemPacakge collection
# Works for ADF and AMS only
try:
    from chemPackage.ams import AMS
    class polarizability(AMS):
        def __init__(self, pol_out):
            AMS.__init__(self, name=pol_out)
            self._collect(abort=True)

        def atomic_polarizability(self):
        # Local part of Hirshfeld partitioned polarizability
            return self.hirshfeld_induced_dipoles_loc

        def inter_atomic_polarizability(self,parameters:Tuple[float,float]=(1.0,3.0)):
            charge_flow = self.__generate_charge_flow(parameters)
            # print("charge flow")
            # print(charge_flow)
            inter_atomic_polarizability =np.zeros((self.natoms, self.natoms, 3, 3),dtype=float)
            for i in range(self.natoms):
                for j in range(i):
                    for idir in range(3):
                        for jdir in range(3):
                            inter_atomic_polarizability[i,j,idir,jdir] = charge_flow[idir,i,j] *(
                                ANGSTROM2BOHR(
                                    self.coordinates[i,jdir] - self.coordinates[j, jdir]
                                )
                            )
                    inter_atomic_polarizability[j,i,:,:] = inter_atomic_polarizability[i,j,:,:]

            return inter_atomic_polarizability

        def __generate_charge_flow(self,parameters:Tuple[float, float]=(1.0,3.0)):
            dis_matrix= ANGSTROM2BOHR(self.__generate_dis_matrix())
            cos_matrix = self.__generate_cos_matrix()
            print(cos_matrix)
            charge_flow = np.zeros((3,self.natoms,self.natoms),dtype=float)
            for idir in range(3):
                charge_flow[idir] = self.__solve_lagrange(dis_matrix,cos_matrix[idir],idir,parameters)

            return charge_flow

        def __generate_dis_matrix(self) -> NDArray[np.float64]:
            # Generate atom distance matrix R_AB
            dis_matrix = np.zeros((self.natoms, self.natoms),dtype=float)
            for i in range(self.natoms):
                for j in range(i):
                    dis_matrix[i][j] = np.linalg.norm(self.coordinates[i]-self.coordinates[j])
                    dis_matrix[j][i] = dis_matrix[i][j]
            return dis_matrix

        def __generate_cos_matrix(self) -> NDArray[np.float64]:
            cos_matrix = np.zeros((3,self.natoms, self.natoms))
            for i in range(self.natoms):
                for j in range(i):
                    for ii in range(3):
                        cos_matrix[ii,i,j] = abs(np.dot((self.coordinates[i] -
                            self.coordinates[j]), np.identity(3)[ii]) /
                            np.linalg.norm(self.coordinates[i]-self.coordinates[j]))
                        cos_matrix[ii,j,i] = cos_matrix[ii,i,j]
            return cos_matrix

        def __solve_lagrange(self,dis_matrix,cos_matrix,idir,parameters:Tuple[float,float]=(1.0,3.0)):
            lambda_matrix = np.zeros((self.natoms, self.natoms),dtype = float)
            for i in range(self.natoms):
                for j in range(i):
                    # NOTE: I don't know why we are using Decimal here.
                    lambda_matrix[i,j] = Decimal(inter_atomic_lambda(parameters,self.atoms[i],
                                                                     self.atoms[j], 
                                                                     dis_matrix[i,j],
                                                                     cos_matrix[i,j]))
                    lambda_matrix[j,i] = lambda_matrix[i,j]
                # lambda_matrix[i,i] = Decimal(-1*sum(lambda_matrix[i,:])) # Add an arbitrary constant C = 1.0
            for i in range(self.natoms):
                lambda_matrix[i,i] = Decimal(-1*sum(lambda_matrix[i]))
            lambda_matrix += 1.0

            lambda_solution = np.linalg.solve(lambda_matrix, self.hirshfeld_induced_charges.T[idir])

            charge_flow = np.zeros((self.natoms, self.natoms), dtype = object)

            for i in range(self.natoms):
                for j in range(i):
                    charge_flow[i, j] = calc_charge_flow(parameters, self.atoms[i], self.atoms[j],
                                                         dis_matrix[i,j], cos_matrix[i,j],
                                                         lambda_solution[i], lambda_solution[j])
                    charge_flow[j,i] = -charge_flow[i,j]
            return charge_flow


except ImportError:
    # This "mini chem" is not written in the same form as in chemPackage for simplicity
    class polarizability():
        def __init__(self, filename) -> None:
            self.nfreq = 0
            f = open(filename)
            pass

        def _collect_polarizability(self,lines):
            self.polarizability = np.array((self.nfreq, 3, 3), dtype = float)
