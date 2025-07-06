try:
    from chemPackage import constants
    ANGSTROM2BOHR = constants.ANGSTROM2BOHR
except ImportError:
    def ANGSTROM2BOHR(num):
        return
