from ramanbond import *
from ramanbond.local import *
from ramanbond.normalmode import *
import pandas as pd
freqout = 'freq_Ag8_co2_DFT.out'

#normalmode_mbh for mobile block hessian
vib = normalmode_mbh(freqout)
#vib = normalmode(freqout)

def collect(mode):
  ## collect infor from AOResponse output files and calculate pol bonds
  pfile = 'mode{:.2f}-p.out'.format(float(mode))
  mfile = pfile[:-6] + '-m.out'
  p = polbond(collectaoresponse(pfile)[0])
  m = polbond(collectaoresponse(mfile)[0])
  ## calculate the step size used in the numerical differentiation of pol wrt. the vibration
  stepsize = calc_stepsize(vib, mode)
  ## dump ramanatom and ramanbond to a pickle file
  raman = ramanbond(p, m, stepsize, 1, 3, 'zz')

  ## get GS hirshfled charge
  charge=read_gs_hirshfeld_charge(freqout)
  gs_qinter = calc_qinter(raman.atomnote, ['C', 'H', 'O', 'N'], ['Ag', 'Au'],charge)

  ## get -p,-m file charge for inter-frag charge flow: Qinter
  pQ=p.calc_Q3compo(['C','H','O','N'],['Ag','Au'])
  mQ=m.calc_Q3compo(['C','H','O','N'],['Ag','Au'])

  return gs_qinter, pQ, mQ


def call_ramanbond():
  modelist = [145.64, 2010.69]
  gs_qinter, pQ, mQ= collect(modelist[0])

  qinter=(np.abs(pQ)+np.abs(mQ))/2
  gs_qinter=np.abs(gs_qinter)
  return gs_qinter, qinter

def print_results():
    # Call the call_ramanbond function to get the data
    gs_qinter, qinter = call_ramanbond()
    
    # Print the results in the specified format
    print(f"GS_qinter: {gs_qinter:.5f}, qinter: {qinter:.5f}")

# Call the function to execute and print the results
print_results()
