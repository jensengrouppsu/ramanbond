# ramanbond
The python package processes outputs of ADF or BAND by SCM and generate Raman bonds, which provides an intuitive interpretation of Raman spectra based on bonding networks of studied systems. More details about the model (Raman bond model) are available from J. Chem. Phys., 152, 024126 and J. Chem. Phys., 153, 224704. The model can be extended to other computational chemistry softwares as long as atomic charges and atomic dipoles are given.

The core functions for Raman bond model are in \_\_init\_\_.py. The class for localized systems calculated by ADF is given in polbond.py. The class for 2d periodic systems calculated by BAND is in given in staticramanbond\_periodic2d.py (currently only static Raman is considered). The scripts normalmode.py analyze the normal modes of localized systems and periodic systems calculated by ADF, BAND, DFTB and help generate input files needed for numerically differentiating polarizabilities. The script plotatombond.py, maps Raman bonds to bonding networks using Pymol.

# Dependencies
The package depends on python packages numpy, math, cmath, decimal, and pickle. The plotting of Raman bonds depends on Pymol. 

# Usage
template ADF/AMS run files are in template folder
example files showing how to obtain Raman bond model variables are in example folder
To use example files, copy them outside of the ramanbond/
when excuating these python and pymol files, they will search for the ramanbond/ in the directory

Example 1. python extract_qinter.py
Print out the inter-fragment charge flow, |qinter|, which was defined in J. Chem. Phys., 14, 2022

Example 2. python plot_3compo.py
Calc and plot the pmol, pinter, pclu, which is the projected contributions of the molecular, of the inter-fragment and of the cluster contribution, respectively, to the Raman intensity. These were defined in Equation(2) in J. Chem. Phys. 153, 224704 (2020)


Example 3a. python calc_ramanbond.py
In Raman bond model, polarizablity derivaties can be written as changes in the atomic induced charge densities R^{atom} and modulation of charge flows between atoms R^{bond}. Both are calculated and stored as pickle file xxx.p to be read in the next step.
R^{atom} and R^{bond} were defined in Equation (10) in J. Chem. Phys. 152, 024126 (2020).

Example 3b. pymol plot_ramanbond.pml
R^{atom} and R^{bond} are plotted/rendered as spheres and cylinders in pymol.
*static Raman bond images can be found in J. Chem. Phys. 152, 024126 (2020)
*frequency-dependent Raman bond can be found in J. Chem. Phys. 153, 224704 (2020)

# additional Notes
The part of the package for periodic systems is still under development.

