from ramanbond import *
from ramanbond.local import *

## collect infor from AOResponse output files and calculate pol bonds
mode=2010.69
pfile = f'mode{mode:.2f}-p.out'
mfile = pfile[:-6] + '-m.out'
p = polbond(collectaoresponse(pfile)[0])
m = polbond(collectaoresponse(mfile)[0])

## calculate the step size used in the numerical differentiation of pol wrt. the vibration
from ramanbond.normalmode import *
freqout = 'freq_Ag8_co2_Dft.out'
vib = normalmode_mbh(freqout)
#vib = normalmode(freqout)
stepsize = calc_stepsize(vib, mode)

## dump ramanatom and ramanbond to a pickle file
raman = ramanbond(p, m, stepsize, 1, 3, 'zz') 
printcoord('Ag8_co2.xyz', raman.atomnote, raman.coord)
dumpdata('Ag8_co2_mode.p', raman.atomnote, raman.coord, raman.ramanatom, raman.ramanbond)

