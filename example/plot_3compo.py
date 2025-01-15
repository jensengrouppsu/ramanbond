from ramanbond import *
from ramanbond.local import *
from ramanbond.normalmode import *

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
  ## calc the (mol, inter, metal) contributions to the Pol der
  mol, inter, metal = calc_3compo(raman.ramanatom, raman.ramanbond, raman.atomnote, ['C', 'H', 'O', 'N'], ['Ag', 'Au'])
  total = mol + inter + metal
  mol_proj = np.dot(np.array([mol.real, mol.imag]), np.array([total.real, total.imag])) / np.linalg.norm(np.array([total.real, total.imag]))
  inter_proj = np.dot(np.array([inter.real, inter.imag]), np.array([total.real, total.imag])) / np.linalg.norm(np.array([total.real, total.imag]))
  metal_proj = np.dot(np.array([metal.real, metal.imag]), np.array([total.real, total.imag])) / np.linalg.norm(np.array([total.real, total.imag]))
  return mol_proj, inter_proj, metal_proj


## plot (mol, inter, metal) and total versus the vibrational frequencies
modelist = [145.64, 2010.69]
mollist = np.zeros(len(modelist))
interlist = np.zeros(len(modelist))
metallist = np.zeros(len(modelist))

for i in range(len(modelist)):
  mollist[i], interlist[i], metallist[i] = collect(modelist[i])

x = np.linspace( min(modelist) - 50, max(modelist) + 50, 4000)
y_mol = sum_lorentzian(x, modelist, mollist, 20)
y_inter = sum_lorentzian(x, modelist, interlist, 20)
y_metal = sum_lorentzian(x, modelist, metallist, 20)
y_total = y_mol + y_inter + y_metal

import matplotlib.pyplot as plt
fig, axs = plt.subplots(4, 1, sharex=True, sharey=False)
axs[0].plot(x, y_mol, 'C0', label='mol')
axs[0].legend()
axs[1].plot(x, y_inter, 'C1', label='inter')
axs[1].legend()
axs[2].plot(x, y_metal, 'C2', label='metal')
axs[2].legend()
axs[3].plot(x, y_total, 'C3', label='total')
axs[3].legend()
plt.show()
