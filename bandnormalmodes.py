import numpy as np
import math
from decimal import *
from numpy.linalg import norm
from cmath import polar

A2B = 1.8897261328856432 ## constant converting angstrom to bohr

atomnote2mass = {
'H': 1.00797,
'C': 12.011,
'N': 14.0067,
'O': 15.9994,
'Mg': 24.305,
'S': 32.06,
'Ag': 107.868,
'Au': 196.9665}


class bandnormalmode():
  def __init__(self, filename):
    f = open(filename)
    f1 = f.readlines()
    f.close()

    self.vec = np.zeros((2, 3)) ## 2d periodic systems
    for i in range(len(f1)):
      if 'Lattice vectors' in f1[i]:
        for j in range(self.vec.shape[0]):
          self.vec[j][0] = float(f1[i+1+j].strip('\n').split()[-3])
          self.vec[j][1] = float(f1[i+1+j].strip('\n').split()[-2])
          self.vec[j][2] = float(f1[i+1+j].strip('\n').split()[-1])

      if 'Total nr. of atoms:' in f1[i]:
        self.num = int(f1[i].strip('\n').split()[-1])
        break

    self.coord = np.zeros((self.num, 3))
    self.atomnote = []
    self.freq = []
    for i in range(len(f1)):
      if 'Geometry' in f1[i]:
        for j in range(self.num):
          self.coord[j][0] = float(f1[i+4+j].strip('\n').split()[-3])
          self.coord[j][1] = float(f1[i+4+j].strip('\n').split()[-2])
          self.coord[j][2] = float(f1[i+4+j].strip('\n').split()[-1])
          self.atomnote.append(f1[i+4+j].strip('\n').split()[-4])

      if 'Normal Mode Frequencies' in f1[i]:
        ii = 3
        while True:
          if len(f1[i+ii]) < 3:
            break
          self.freq.append(float(f1[i+ii].strip('\n').split()[-1]))
          ii += 1

    self.normalmodes = np.zeros((len(self.freq), self.num, 3))
    ii = 0
    for i in range(len(f1)):
      if 'Displacements (x/y/z)' in f1[i]: ## not mass weighted, in Bohr
        for j in range(self.num):
          self.normalmodes[ii][j][0] = float(f1[i+1+j].strip('\n').split()[-3])
          self.normalmodes[ii][j][1] = float(f1[i+1+j].strip('\n').split()[-2])
          self.normalmodes[ii][j][2] = float(f1[i+1+j].strip('\n').split()[-1])
        ii += 1


    self.masses = np.zeros(self.num)
    for i in range(len(self.atomnote)):
      self.masses[i] = atomnote2mass[self.atomnote[i]]

    self.rt_masses = np.diag(np.sqrt(self.masses))
    self.mw_normalmodes = np.zeros(self.normalmodes.shape)

    for i in range(self.normalmodes.shape[0]):
      self.mw_normalmodes[i] = np.matmul(self.rt_masses, self.normalmodes[i])

  def calc_stepsize(self, mode):
    for i in range(self.normalmodes.shape[0]):
      if abs(float(mode) - self.freq[i]) < 0.02:
        self.stepsize = (norm(self.mw_normalmodes[i]) / norm(self.normalmodes[i]) * 0.01 ) ##cart_stepsize = 0.01 usually 
      else:
        pass



  def create_modefile(self, mode, template, pfield=0.001, mfield=-0.001):
    for i in range(len(self.freq)):
      if abs(float(mode) - self.freq[i]) < 0.02:
        self.displacedcoord_plus = self.coord + (self.normalmodes[i] / ( A2B * norm(self.normalmodes[i] ))) * 0.01
## convert from Bohr to Angstrom !
        self.displacedcoord_minus = self.coord - (self.normalmodes[i] / ( A2B * norm(self.normalmodes[i]))) * 0.01

    f = open(template)
    f1 = f.readlines()
    f.close()
    
    g = open('mode'+'{0:.3f}'.format(mode)+'-vibpfieldp.run', 'w')
    for i in range(len(f1)):
      g.write(f1[i])
      if 'atoms' in f1[i].lower():
        for j in range(self.num):
          g.write('{0:<3}    {1: .8f}    {2: .8f}    {3: .8f}\n'.format(self.atomnote[j], self.displacedcoord_plus[j][0], self.displacedcoord_plus[j][1], self.displacedcoord_plus[j][2]))
      if 'efield' in f1[i].lower():
         g.write('unit a.u.\n')
         g.write('ez {0: .4f}\n'.format(pfield))
      if 'lattice' in f1[i].lower():
        for j in range(self.vec.shape[0]):
          g.write('{0: .8f}    {1: .8f}    {2: .8f}\n'.format(self.vec[j][0], self.vec[j][1], self.vec[j][2]))

    g.close()

    g = open('mode'+'{0:.3f}'.format(mode)+'-vibpfieldm.run', 'w')
    for i in range(len(f1)):
      g.write(f1[i])
      if 'atoms' in f1[i].lower():
        for j in range(self.num):
          g.write('{0:<3}    {1: .8f}    {2: .8f}    {3: .8f}\n'.format(self.atomnote[j], self.displacedcoord_plus[j][0], self.displacedcoord_plus[j][1], self.displacedcoord_plus[j][2]))
      if 'efield' in f1[i].lower():
         g.write('unit a.u.\n')
         g.write('ez {0: .4f}\n'.format(mfield))

      if 'lattice' in f1[i].lower():
        for j in range(self.vec.shape[0]):
          g.write('{0: .8f}    {1: .8f}    {2: .8f}\n'.format(self.vec[j][0], self.vec[j][1], self.vec[j][2]))
    g.close()

    
    g = open('mode'+'{0:.3f}'.format(mode)+'-vibmfieldp.run', 'w')
    for i in range(len(f1)):
      g.write(f1[i])
      if 'atoms' in f1[i].lower():
        for j in range(self.num):
          g.write('{0:<3}    {1: .8f}    {2: .8f}    {3: .8f}\n'.format(self.atomnote[j], self.displacedcoord_minus[j][0], self.displacedcoord_minus[j][1], self.displacedcoord_minus[j][2]))
      if 'efield' in f1[i].lower():
         g.write('unit a.u.\n')
         g.write('ez {0: .4f}\n'.format(pfield))

      if 'lattice' in f1[i].lower():
        for j in range(self.vec.shape[0]):
          g.write('{0: .8f}    {1: .8f}    {2: .8f}\n'.format(self.vec[j][0], self.vec[j][1], self.vec[j][2]))
    g.close()

    g = open('mode'+'{0:.3f}'.format(mode)+'-vibmfieldm.run', 'w')
    for i in range(len(f1)):
      g.write(f1[i])
      if 'atoms' in f1[i].lower():
        for j in range(self.num):
          g.write('{0:<3}    {1: .8f}    {2: .8f}    {3: .8f}\n'.format(self.atomnote[j], self.displacedcoord_minus[j][0], self.displacedcoord_minus[j][1], self.displacedcoord_minus[j][2]))
      if 'efield' in f1[i].lower():
         g.write('unit a.u.\n')
         g.write('ez {0: .4f}\n'.format(mfield))

      if 'lattice' in f1[i].lower():
        for j in range(self.vec.shape[0]):
          g.write('{0: .8f}    {1: .8f}    {2: .8f}\n'.format(self.vec[j][0], self.vec[j][1], self.vec[j][2]))
    g.close()


  def modes_to_pymol(self, freqmin=400., freqmax=2000.,
    vectorwidth=0.1, scalevector=1.0, component='all', transparency=1.0, unitconvert=1.0):

    # Allow for degenerate modes
    deg_list = ('', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i')
    degeneracy = 0
    previous_mode = ''

    # Cycle over each normal mode
    for i in range(len(self.freq)):

      # Make sure that this mode is within range
      if self.freq[i] < freqmin: continue
      if self.freq[i] > freqmax: continue

      # Check for degeneracy
      strmode = '{0:.3f}'.format(self.freq[i])
      if previous_mode == strmode:
        degeneracy += 1
        strmode = strmode+'_'+deg_list[degeneracy]
      else:
        degeneracy = 0

      # Open filename for this mode
      fw = open('mode{0}.pymol'.format(strmode), 'w')
      previous_mode = strmode

      pymolmode = self.normalmodes[i] / ( A2B * norm(self.normalmodes[i] ))

      # Cycle over each atom
      for j in range(self.coord.shape[0]):
        if component == 'all':
          print('cgo_modevec {0}{2: 11.7f},{3: 11.7f},{4: 11.7f}'
                  '{1}, {0}{5: 11.7f},{6: 11.7f},{7: 11.7f}{1}, radius={8:4.2f}, transparency={9:2.2f}'.format(
                  '[', ']', unitconvert*self.coord[j][0], unitconvert*self.coord[j][1], unitconvert*
                  self.coord[j][2], unitconvert*pymolmode[j][0]*scalevector,
                  unitconvert*pymolmode[j][1]*scalevector, unitconvert*pymolmode[j][2]*scalevector,
                  vectorwidth, transparency), file=fw)

        if component == 'x':
          print('cgo_modevec {0}{2: 11.7f},{3: 11.7f},{4: 11.7f}'
                  '{1}, {0}{5: 11.7f},{6: 11.7f},{7: 11.7f}{1}, radius={8:4.2f}, transparency={9:2.2f}'.format(
                  '[', ']', unitconvert*self.coord[j][0], unitconvert*self.coord[j][1], unitconvert*
                  self.coord[j][2], unitconvert*pymolmode[j][0]*scalevector,
                  0.0, 0.0,
                  vectorwidth, transparency), file=fw)
        if component == 'y':
          print('cgo_modevec {0}{2: 11.7f},{3: 11.7f},{4: 11.7f}'
                  '{1}, {0}{5: 11.7f},{6: 11.7f},{7: 11.7f}{1}, radius={8:4.2f}, transparency={9:2.2f}'.format(
                  '[', ']', unitconvert*self.coord[j][0], unitconvert*self.coord[j][1], unitconvert*
                  self.coord[j][2], 0.0, unitconvert*pymolmode[j][1]*scalevector,
                  0.0,
                  vectorwidth, transparency), file=fw)
        if component == 'z':
          print('cgo_modevec {0}{2: 11.7f},{3: 11.7f},{4: 11.7f}'
                  '{1}, {0}{5: 11.7f},{6: 11.7f},{7: 11.7f}{1}, radius={8:4.2f}, transparency={9:2.2f}'.format(
                  '[', ']', unitconvert*self.coord[j][0], unitconvert*self.coord[j][1], unitconvert*
                  self.coord[j][2], 0.0, 0.0, unitconvert*pymolmode[j][2]*scalevector,
                  vectorwidth, transparency), file=fw)

      fw.close()

