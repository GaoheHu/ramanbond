# ramanbond
The python package processes outputs of ADF or BAND by SCM and generate Raman bonds, which provides an intuitive interpretation of Raman spectra based on bonding networks of studied systems. More details about the model (Raman bond model) are available from J. Chem. Phys., 152, 024126 and J. Chem. Phys., 153, 224704. The model can be extended to other computational chemistry software as long as atomic charges and atomic dipoles are given.

# Dependencies
The package depends on python packages numpy, math, cmath, decimal, and pickle. The plotting of Raman bonds depends on Pymol. 
At this point, it depends on [chemPackage](https://github.com/jensengrouppsu/chemPackage) for collecting properties from the outputs.
Refer to chemPackage for how to install it.

# Usage

The current version of ramanbond mimic the design of chemPackage. 
All features are built around the "collect" function.
A typical workflow of ramanbond is:
```
from ramanbond import collect

d=collect("output")
d.collect_raman_derivatives() # If you are collecting Raman results
d.plot_bond()
```

# additional Notes

We plan to implement a new collection back end that is independent of chemPackage
