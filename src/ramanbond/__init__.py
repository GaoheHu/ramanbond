from .pol import pol

__version__ = '2.0'
__all__ = [
    'pol',
    'conversion',
    'collect'
]
try:
    from chemPackage import collect
except ImportError:
    def collect(filename):
        d = pol(filename)
        return d
