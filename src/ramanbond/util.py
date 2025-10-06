import math
from typing import Tuple 
from decimal import getcontext,Decimal
from .conversion import ANGSTROM2BOHR
import math

# Covalent atomic radii in Angstrom
RBS_DICT_ANGSTROM: dict[str, float] = {
'H': 0.31,'He': 0.28,'Li': 1.28,'Be': 0.96,'B': 0.84,'C': 0.70,'N': 0.71,'O':
    0.66,'F': 0.57,'Ne': 0.58,'Na': 1.66, 'Mg': 1.41,'Al': 1.21,'Si': 1.11,'P':
    1.07,'S': 1.05,'Cl': 1.02,'Ar': 1.06, 'Cu':1.38, 'Zn': 1.31, 'Ag':
    1.45,'Au': 1.36
}

RBS_DICT_BOHR: dict[str, float] = {}

for key in RBS_DICT_ANGSTROM:
    RBS_DICT_BOHR[key] = ANGSTROM2BOHR(RBS_DICT_ANGSTROM[key])


def penalty_fucntion(parameters: Tuple[float,float], Za, Zb, distance,
                     anglecos):
    # TODO: Add in reference
    # The penalty function for LoProp.
    sum_cov_distance = RBS_DICT_BOHR[Za] + RBS_DICT_BOHR[Zb]
    try:
        if distance < 1.1*sum_cov_distance:
            return (math.exp(parameters[0]*(distance/sum_cov_distance)**2) +
                    math.exp(parameters[1]*(1-anglecos)))
        else:
            return (math.exp(2*parameters[0]*(distance/sum_cov_distance)**2) +
                    math.exp(parameters[1]*(1-anglecos)))
    except OverflowError:
        return float("inf")

def inter_atomic_lambda(parameters:Tuple[float,float], atom1, atom2, distance, angle_cos):
    # Return lambda between atom a and atom b
    # To constrcut a matrix lambda
    try:
        if distance == 0:
            return 0.0
        else:
            before_L = 1/(2*penalty_fucntion(parameters, atom1, atom2, distance, angle_cos))
    except OverflowError:
        before_L = 0.0

    return before_L

def calc_charge_flow(parameters:Tuple[float,float], atom1, atom2, distance,
                     angle_cos, lambda_i, lambda_j):
    try:
        charge_flow = -1 *( (lambda_i - lambda_j) /
                        (2*penalty_fucntion(parameters, atom1, atom2, distance,
                                        angle_cos)))
    except OverflowError:
        charge_flow = 0.0

    return charge_flow

def find_pm_results(dir=None):

    '''
    Collects all the polarizability files under given directory.
    By default the directory is cwd.
    Greatly inspired by chemPackage
    '''
    from glob import glob
    from natsort import natsort_key
    import os

    if dir is None:
        dir = os.curdir

    pfiles = glob(os.path.join(dir, 'mode*-p.out'))
    mfiles = glob(os.path.join(dir, 'mode*-m.out'))

    assert len(mfiles) == len(pfiles), ("Numbers of plus and minus files are note the same")

    mfiles.sort(key=natsort_key)
    pfiles.sort(key=natsort_key)

    return pfiles, mfiles

def calc_3component(normalmode,polarizability_derivative_atoms,
                    polarizability_derivative_bonds, mol_list =
                    ["C","H","O","N"],clus_list=["Ag","Au"]):

    total_contributions = []
    for mode in range(normalmode.nmodes):
        contributions= {"MOLECULE": 0.0, "CLUSTER":0.0, "INTER":0.0}

        atom_group = {}
        for atom in normalmode.atoms:
            if atom in mol_list:
                atom_group[atom] = "MOLECULE"
            elif atom in clus_list:
                atom_group[atom] = "CLUSTER"
            else:
                print("Atom {} not in either of molecule list or cluster list".format(atom))
                return
        for i, atom in enumerate(normalmode.atoms):
            group = atom_group[atom]
            contributions[group]+= polarizability_derivative_atoms[mode,i]

        for i in range(normalmode.natoms):
            for j in range(i):
                group_i = atom_group[normalmode.atoms[i]]
                group_j = atom_group[normalmode.atoms[j]]

                if group_i == group_j:
                    contributions[group_i] += polarizability_derivative_bonds[mode,i,j]
                else:
                    contributions["INTER"] += polarizability_derivative_bonds[mode,i,j]

        total_contributions.append(contributions)

    return total_contributions

def angle_to_rgb(angle):
    # Normalize angle into [-π, π]
    angle = ((angle + math.pi) % (2 * math.pi)) - math.pi

    # Map angle to range [-3, 3]
    p = (3.0 * angle) / math.pi

    if -3 <= p < -2:      # Red → Yellow
        R, G, B = 1, p + 3, 0
    elif -2 <= p < -1:    # Yellow → Green
        R, G, B = -1 - p, 1, 0
    elif -1 <= p < 0:     # Green → Cyan
        R, G, B = 0, 1, p + 1
    elif 0 <= p < 1:      # Cyan → Blue
        R, G, B = 0, 1 - p, 1
    elif 1 <= p < 2:      # Blue → Magenta
        R, G, B = p - 1, 0, 1
    elif 2 <= p <= 3:     # Magenta → Red
        R, G, B = 1, 0, 3 - p
    else:
        raise ValueError(f"Unexpected p value: {p}")

    return float(R), float(G), float(B)
