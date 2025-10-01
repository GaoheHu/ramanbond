import os, sys
import numpy as np
import math, cmath
from typing import Tuple
import textwrap

from .util import find_pm_results, angle_to_rgb
from .pol import polarizability_bond, AMS
from .conversion import ANGSTROM2BOHR

class raman_bond(polarizability_bond):

    def __init__(self,outputfile,parameters:Tuple[float,float]=(1.0,3.0)):
        self.chemObj = AMS(name=outputfile)
        self.parameters = parameters
        self.chemObj._collect(abort =False)
        # self.atoms = self.fragment_atoms
        # self.coordinates = self.fragment_coordinates
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

        # If you don't have fragment atoms in frequency calculation
        if self.fragment_atoms is None:
            self.fragment_atoms = p.fragment_atoms
            self.fragment_coordinates = 0.5 * (p.fragment_coordinates + m.fragment_coordinates)

        self.raman_atoms = raman_atoms
        self.raman_bonds = raman_bonds

    def plot_bond(self, numbond=70, SCALE=1, component="all"):
        # if mode not in self.v_frequencies:
        #     print("The requested normal mode is not calculated")

        dirs = {"xx":0, "yy":1, "zz":2}
        if component in dirs:
            raman_atom = self.raman_atoms[:,:, dirs[component], dirs[component]]
            raman_bond = self.raman_bonds[:,:,:,dirs[component], dirs[component]]
        else:
            print("Considering the isotropic component")
            raman_atom = np.trace(self.raman_atoms, axis1=2, axis2=3)
            raman_bond = np.trace(self.raman_bonds, axis1=3, axis2=4)
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
        for n, mode in enumerate(self.v_frequencies):
            with open(self.filename[:-4]+".mode{:.2f}".format(mode)+".pml", "w") as f:
                sys.stdout = f
                print(textwrap.dedent("""
                load {}.xyz
                preset.ball_and_stick, all 
                set sphere_scale, 0.00001, all
                set_bond stick_radius, 0.0000000001, all
                """.format(self.filename[:-4] )))

                counter = 0
                for i in range(self.natoms):
                    r_atom = ((abs(raman_atom[n,i]/SCALE)*3 / 
                        (4*math.pi*(ANGSTROM2BOHR(1))**3))**(1/3))
                    print("set sphere_scale, {}, {}/{}".format(r_atom, i+1, self.atoms[i]))
                    
                    R, G, B = angle_to_rgb(cmath.polar(raman_atom[n,i])[1])
                    print("set_color atomcolor{}, [{},{},{}]".format(i+1, R, G, B))
                    print("color atomcolor{}, {}/{}".format(i+1, i+1, self.atoms[i]))
                    for j in range(i):
                        if cmath.polar(raman_bond[n, i,j])[0] > threshold:
                            counter+=1
                            r_bond = math.sqrt((abs(raman_bond[n,i,
                                j]/SCALE)/(self.dis_matrix[i,j]*math.pi*(ANGSTROM2BOHR**3.0))))
                            print("bond {}/{}, {}/{}".format(i+1,self.atoms[i],  j+1, self.atoms[j]))
                            print("select bond{}, {}/{} {}/{}".format(counter, i+1, self.atoms[i], j+1, self.atoms[j]))
                            print("set_bond stick_radius, {}, bond{}".format(r_bond, counter))
                            R,G,B = angle_to_rgb(cmath.polar(raman_bond[n,i,j])[1])
                            print("set_color bondcolor{}, [{}, {},{}]".format(counter, R,G,B))
                            print("set_bond stick_color, bondcolor{}, bond{}".format(counter, counter))
                print("show sticks, all")
        return
    def __generate_dis_matrix(self):
        # Generate atom distance matrix R_AB
        dis_matrix = np.zeros((self.natoms, self.natoms),dtype=float)
        for i in range(self.natoms):
            for j in range(i):
                dis_matrix[i][j] = np.linalg.norm(self.coordinates[i]-self.coordinates[j])
                dis_matrix[j][i] = dis_matrix[i][j]
        self.dis_matrix = dis_matrix
        return 
