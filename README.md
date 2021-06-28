# ramanbond
The python package processes outputs of ADF or BAND by SCM and generate Raman bonds, which provides an intuitive interpretation of Raman spectra based on bonding networks of studied systems. More details about the model (Raman bond model) are available from J. Chem. Phys., 152, 024126 and J. Chem. Phys., 153, 224704. The model can be extended to other computational chemistry softwares as long as atomic charges and atomic dipoles are given.

The core functions for Raman bond model are in \_\_init\_\_.py. The class for localized systems calculated by ADF is given in polbond.py. The class for 2d periodic systems calculated by BAND is in given in staticramanbond_periodic2d.py (currently only static Raman is considered). The scripts normalmode.py analyze the normal modes of localized systems and periodic systems calculated by ADF, BAND, DFTB and help generate input files needed for numerically differentiating polarizabilities. The scripts geometry_periodic.py helps to generate the coordinates of atoms of periodic systems in customizable supercells. The script plotatombond.py, maps Raman bonds to bonding networks using Pymol.

# Dependencies
The package depends on python packages numpy, math, cmath, decimal, and pickle. The plotting of Raman bonds depends on Pymol. 

# additional Notes
The part of the package for periodic systems is still under development.

