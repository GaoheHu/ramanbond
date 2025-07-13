
__version__ = '2.0'
__all__ = [
    'pol',
    'normalmode',
    'conversion',
    'collect'
]
try:
    from chemPackage import collect as chemCollect
    def collect(filename):
        d = chemCollect(filename)
        if "FREQUENCIES" in d.calctype:
            from .normalmode import normalmode
            return normalmode(filename)
        elif "POLARIZABILITY" in d.calctype:
            from .pol import polarizability
            return polarizability(filename)
        else:
            print("Calculation not supported in Ramanbond")

except ImportError:
    def collect(filename):
        print("Not implemented yet. Please install chemPackage.")
