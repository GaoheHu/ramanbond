import numpy as np
from numpy._typing import NDArray
from typing import Tuple
from decimal import Decimal
import math, cmath
import textwrap

from .util import inter_atomic_lambda, calc_charge_flow
from .conversion import ANGSTROM2BOHR
from .util import angle_to_rgb
# This is basically a minimal version of chemPackage collection
# Works for ADF and AMS only
try:
    from chemPackage.ams import AMS

except ImportError:
    # This "mini chem" is not written in the same form as in chemPackage for simplicity
    class AMS():
        def __init__(self, name):
            print("Not yet implemented")
        # def _collect(abort=False):
        #     return
    
class polarizability_bond():
    def __init__(self, outputfile, parameters:Tuple[float,float]=(1.0,3.0)):
        self.chemObj=AMS(name=outputfile)
        self.parameters = parameters
        self.chemObj._collect(abort=False)
        self.atomic_polarizability = self._calc_atomic_polarizability()

        # Overwrite initial geometry with fragment geometry
        self.atoms = self.fragment_atoms
        self.coordinates = self.fragment_coordinates
        if "VELOCITY" not in self.calctype:
            self.bond_polarizability = self._calc_bond_polarizability()

    def __getattr__(self, name: str, /):
        if hasattr(self.chemObj, name):
            return getattr(self.chemObj, name)
        raise AttributeError(f"{self.__class__.__name__} has no attribute '{name}'")
        
    def plot_bond(self, numbond=70, SCALE=1, component="all"):
        # if mode not in self.v_frequencies:
        #     print("The requested normal mode is not calculated")

        import sys
        dirs = {"xx":0, "yy":1, "zz":2}
        if component in dirs:
            raman_atom = self.atomic_polarizability[:, dirs[component], dirs[component]]
            raman_bond = self.bond_polarizability[:,:,dirs[component], dirs[component]]
        else:
            print("Considering the isotropic component")
            raman_atom = np.trace(self.atomic_polarizability, axis1=1, axis2=2)
            raman_bond = np.trace(self.bond_polarizability, axis1=2, axis2=3)
        # Adjust the sign
        if np.sum(raman_bond).real+np.sum(raman_atom).real <= 0:
            raman_atom = - raman_atom
            raman_bond = - raman_bond 

        # Limit number of bonds shown
        bond_magnitude = np.sort(raman_bond,axis = None)
        if numbond < len(bond_magnitude):
            threshold = bond_magnitude[numbond-1]
        else:
             threshold = bond_magnitude[-1]


        self.coordinates = self.fragment_coordinates
        self.atoms = self.fragment_atoms
        self.writeCoords()
        with open(self.filename[:-4]+".pml", "w") as f:
            sys.stdout = f
            print(textwrap.dedent("""
            load {}.xyz
            preset.ball_and_stick, all 
            set sphere_scale, 0.00001, all
            set stick_radius, 0.0000000001, all
            """.format(self.filename[:-4] )))

            counter = 0
            for i in range(self.natoms):
                # print("sef spere_scale, 0.0000001, {}".format(self.atoms[i]))
                r_atom = ((abs(raman_atom[i]/SCALE)*3 / 
                    (4*math.pi*(ANGSTROM2BOHR(1))**3))**(1/3))
                print("set sphere_scale, {}, {}/{}".format(r_atom, i+1, self.atoms[i]))
                
                R, G, B = angle_to_rgb(cmath.polar(raman_atom[i])[1])
                print("set_color atomcolor{}, [{},{},{}]".format(i+1, R, G, B))
                print("color atomcolor{}, {}/{}".format(i+1, i+1, self.atoms[i]))
                for j in range(i):
                    if cmath.polar(raman_bond[i,j])[0] > threshold:
                        counter+=1
                        r_bond = math.sqrt((abs(raman_bond[i,
                            j]/SCALE)/(self.dis_matrix[i,j]*math.pi*(ANGSTROM2BOHR**3.0))))
                        print("bond {}/{}, {}/{}".format(i+1,self.atoms[i],  j+1, self.atoms[j]))
                        print("select bond{}, {}/{} {}/{}".format(counter, i+1, self.atoms[i], j+1, self.atoms[j]))
                        print("set_bond stick_radius, {}, bond{}".format(r_bond, counter))
                        R,G,B = angle_to_rgb(cmath.polar(raman_bond[i,j])[1])
                        print("set_color bondcolor{}, [{}, {},{}]".format(counter, R,G,B))
                        print("set_bond stick_color, bondcolor{}, bond{}".format(counter, counter))
                print("show sticks, all")
        return

    def _calc_atomic_polarizability(self):
        # Local part of Hirshfeld partitioned polarizability
        atomic_polarizability = \
        np.array([self.hirshfeld_induced_dipoles_loc[i,:,:] for i in
            range(self.natoms)])

        return atomic_polarizability

    def _calc_bond_polarizability(self):
        charge_flow = self.__generate_charge_flow()
        bond_polarizability =np.zeros((self.natoms, self.natoms, 3,3),dtype=np.complex128)
        for i in range(self.natoms):
            for j in range(i):
                for idir in range(3):
                    for jdir in range(3):
                        bond_polarizability[i,j,idir,jdir] = charge_flow[idir,i,j] *(
                                ANGSTROM2BOHR(
                                    self.coordinates[i,jdir] - self.coordinates[j, jdir]
                                )
                        )
                        bond_polarizability[j,i:,:] = bond_polarizability[i,j,:,:]
        return bond_polarizability

    def __generate_charge_flow(self):
        dis_matrix= ANGSTROM2BOHR(self.__generate_dis_matrix())
        cos_matrix = self.__generate_cos_matrix()
        charge_flow = np.zeros((3,self.natoms,self.natoms),dtype=np.complex128)
        for idir in range(3):
            charge_flow[idir] = self.__solve_lagrange(dis_matrix,cos_matrix[idir],idir)

        return charge_flow

    def __generate_dis_matrix(self) -> NDArray[np.float64]:
        # Generate atom distance matrix R_AB
        dis_matrix = np.zeros((self.natoms, self.natoms),dtype=float)
        for i in range(self.natoms):
            for j in range(i):
                dis_matrix[i][j] = np.linalg.norm(self.coordinates[i]-self.coordinates[j])
                dis_matrix[j][i] = dis_matrix[i][j]
        self.dis_matrix = dis_matrix
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

    def __solve_lagrange(self,dis_matrix,cos_matrix,idir):
        lambda_matrix = np.zeros((self.natoms, self.natoms),dtype = float)
        for i in range(self.natoms):
            for j in range(i):
                lambda_matrix[i,j] = \
                inter_atomic_lambda(self.parameters,self.atoms[i],
                                            self.atoms[j], dis_matrix[i,j],
                                            cos_matrix[i,j])
                lambda_matrix[j,i] = lambda_matrix[i,j]
            # lambda_matrix[i,i] = Decimal(-1*sum(lambda_matrix[i,:])) # Add an arbitrary constant C = 1.0
        for i in range(self.natoms):
            lambda_matrix[i,i] = -1*sum(lambda_matrix[i])
        lambda_matrix += 1.0

        lambda_solution = np.linalg.solve(lambda_matrix, self.hirshfeld_induced_charges.T[idir])

        charge_flow = np.zeros((self.natoms, self.natoms), dtype = object)

        for i in range(self.natoms):
            for j in range(i):
                charge_flow[i, j] = calc_charge_flow(self.parameters, self.atoms[i], self.atoms[j],
                                                     dis_matrix[i,j], cos_matrix[i,j],
                                                     lambda_solution[i], lambda_solution[j])
                charge_flow[j,i] = -charge_flow[i,j]
        return charge_flow
