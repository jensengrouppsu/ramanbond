import numpy as np
import math
from decimal import *
from numpy.linalg import norm
from cmath import polar

A2B = 1.8897261328856432 ## constant converting angstrom to bohr

def mode_to_pymol(obj, mode, vectorwidth=0.1, scalevector=1.0, component='all', transparency=1.0, unitconvert=1.0):

  # Open filename for this mode
  fw = open('mode{0: .3f}.pymol'.format(mode), 'w')
  for i in range(len(obj.freq)):
    if abs(obj.freq - mode) < 0.01:
      pymolmode = obj.normalmode[i] / ( A2B * norm(obj.normalmode[i] ))

  # Cycle over each atom
  for j in range(obj.coord.shape[0]):
    if component == 'all':
      print('cgo_modevec {0}{2: 11.7f},{3: 11.7f},{4: 11.7f}'
              '{1}, {0}{5: 11.7f},{6: 11.7f},{7: 11.7f}{1}, radius={8:4.2f}, transparency={9:2.2f}'.format(
              '[', ']', unitconvert*obj.coord[j][0], unitconvert*obj.coord[j][1], unitconvert*
              obj.coord[j][2], unitconvert*pymolmode[j][0]*scalevector,
              unitconvert*pymolmode[j][1]*scalevector, unitconvert*pymolmode[j][2]*scalevector,
              vectorwidth, transparency), file=fw)

    if component == 'x':
      print('cgo_modevec {0}{2: 11.7f},{3: 11.7f},{4: 11.7f}'
              '{1}, {0}{5: 11.7f},{6: 11.7f},{7: 11.7f}{1}, radius={8:4.2f}, transparency={9:2.2f}'.format(
              '[', ']', unitconvert*obj.coord[j][0], unitconvert*obj.coord[j][1], unitconvert*
              obj.coord[j][2], unitconvert*pymolmode[j][0]*scalevector,
              0.0, 0.0,
              vectorwidth, transparency), file=fw)
    if component == 'y':
      print('cgo_modevec {0}{2: 11.7f},{3: 11.7f},{4: 11.7f}'
              '{1}, {0}{5: 11.7f},{6: 11.7f},{7: 11.7f}{1}, radius={8:4.2f}, transparency={9:2.2f}'.format(
              '[', ']', unitconvert*obj.coord[j][0], unitconvert*obj.coord[j][1], unitconvert*
              obj.coord[j][2], 0.0, unitconvert*pymolmode[j][1]*scalevector,
              0.0,
              vectorwidth, transparency), file=fw)
    if component == 'z':
      print('cgo_modevec {0}{2: 11.7f},{3: 11.7f},{4: 11.7f}'
              '{1}, {0}{5: 11.7f},{6: 11.7f},{7: 11.7f}{1}, radius={8:4.2f}, transparency={9:2.2f}'.format(
              '[', ']', unitconvert*obj.coord[j][0], unitconvert*obj.coord[j][1], unitconvert*
              obj.coord[j][2], 0.0, 0.0, unitconvert*pymolmode[j][2]*scalevector,
              vectorwidth, transparency), file=fw)

  fw.close()


def calc_stepsize(obj, mode):
  for i in range(len(obj.freq)):
    if abs(mode - obj.freq[i]) < 0.01:
      return (norm(obj.mw_normalmode[i]) / norm(obj.normalmode[i]) * 0.01 ) 
'''cart_stepsize = 0.01 usually,  normalmodes due to digital precision are not strictly normalzied '''


def calc_displacedcoord(obj, mode):
  for i in range(len(obj.freq)):
    if abs(float(mode) - obj.freq[i]) < 0.01:
      displacedcoord_plus = obj.coord + (obj.normalmode[i] / ( A2B * norm(obj.normalmode[i] ))) * 0.01  ## convert from Bohr to Angstrom !
      displacedcoord_minus = obj.coord - (obj.normalmode[i] / ( A2B * norm(obj.normalmode[i]))) * 0.01 

  return displacedcoord_plus, displacedcoord_minus


def create_modefile(obj, mode, template):
  dispalcedcoord_plus, displacedcoord_minus = calc_displacedcoord(obj, mode)
  f = open(template)
  f1 = f.readlines()
  f.close()
  
  g = open('mode'+'{0:.3f}'.format(mode)+'-p.run', 'w')
  for i in range(len(f1)):
    if 'atoms' in f1[i].strip('\n').lower():
      g.write(f1[i])
      for j in range(obj.num):
        g.write('{0:<3}    {1: .8f}    {2: .8f}    {3: .8f}\n'.format(obj.atomnote[j], displacedcoord_plus[j][0], displacedcoord_plus[j][1], displacedcoord_plus[j][2]))
    else:
      g.write(f1[i])
  g.close()

  g = open('mode'+'{0:.3f}'.format(mode)+'-m.run', 'w')
  for i in range(len(f1)):
    if 'atoms' in f1[i].strip('\n').lower():
      g.write(f1[i])
      for j in range(obj.num):
        g.write('{0:<3}    {1: .8f}    {2: .8f}    {3: .8f}\n'.format(obj.atomnote[j], displacedcoord_minus[j][0], displacedcoord_minus[j][1], displacedcoord_minus[j][2]))
    else:
      g.write(f1[i])
  g.close()


def create_periodic_numdif_modefile(obj, mode, template, pfield=0.001, mfield=-0.001):
  displacedcoord_plus, displacedcoord_minus = calc_displacedcoord(obj, mode)

  f = open(template)
  f1 = f.readlines()
  f.close()
  
  g = open('mode'+'{0:.3f}'.format(mode)+'-vibpfieldp.run', 'w')
  for i in range(len(f1)):
    g.write(f1[i])
    if 'atoms' in f1[i].lower():
      for j in range(obj.num):
        g.write('{0:<3}    {1: .8f}    {2: .8f}    {3: .8f}\n'.format(obj.atomnote[j], displacedcoord_plus[j][0], displacedcoord_plus[j][1], displacedcoord_plus[j][2]))
    if 'efield' in f1[i].lower():
       g.write('unit a.u.\n')
       g.write('ez {0: .4f}\n'.format(pfield))
    if 'lattice' in f1[i].lower():
      for j in range(obj.vec.shape[0]):
        g.write('{0: .8f}    {1: .8f}    {2: .8f}\n'.format(obj.vec[j][0], obj.vec[j][1], obj.vec[j][2]))

  g.close()

  g = open('mode'+'{0:.3f}'.format(mode)+'-vibpfieldm.run', 'w')
  for i in range(len(f1)):
    g.write(f1[i])
    if 'atoms' in f1[i].lower():
      for j in range(obj.num):
        g.write('{0:<3}    {1: .8f}    {2: .8f}    {3: .8f}\n'.format(obj.atomnote[j], displacedcoord_plus[j][0], displacedcoord_plus[j][1], displacedcoord_plus[j][2]))
    if 'efield' in f1[i].lower():
       g.write('unit a.u.\n')
       g.write('ez {0: .4f}\n'.format(mfield))

    if 'lattice' in f1[i].lower():
      for j in range(obj.vec.shape[0]):
        g.write('{0: .8f}    {1: .8f}    {2: .8f}\n'.format(obj.vec[j][0], obj.vec[j][1], obj.vec[j][2]))
  g.close()

  
  g = open('mode'+'{0:.3f}'.format(mode)+'-vibmfieldp.run', 'w')
  for i in range(len(f1)):
    g.write(f1[i])
    if 'atoms' in f1[i].lower():
      for j in range(obj.num):
        g.write('{0:<3}    {1: .8f}    {2: .8f}    {3: .8f}\n'.format(obj.atomnote[j], displacedcoord_minus[j][0], displacedcoord_minus[j][1], displacedcoord_minus[j][2]))
    if 'efield' in f1[i].lower():
       g.write('unit a.u.\n')
       g.write('ez {0: .4f}\n'.format(pfield))

    if 'lattice' in f1[i].lower():
      for j in range(obj.vec.shape[0]):
        g.write('{0: .8f}    {1: .8f}    {2: .8f}\n'.format(obj.vec[j][0], obj.vec[j][1], obj.vec[j][2]))
  g.close()

  g = open('mode'+'{0:.3f}'.format(mode)+'-vibmfieldm.run', 'w')
  for i in range(len(f1)):
    g.write(f1[i])
    if 'atoms' in f1[i].lower():
      for j in range(obj.num):
        g.write('{0:<3}    {1: .8f}    {2: .8f}    {3: .8f}\n'.format(obj.atomnote[j], displacedcoord_minus[j][0], displacedcoord_minus[j][1], displacedcoord_minus[j][2]))
    if 'efield' in f1[i].lower():
       g.write('unit a.u.\n')
       g.write('ez {0: .4f}\n'.format(mfield))

    if 'lattice' in f1[i].lower():
      for j in range(obj.vec.shape[0]):
        g.write('{0: .8f}    {1: .8f}    {2: .8f}\n'.format(obj.vec[j][0], obj.vec[j][1], obj.vec[j][2]))
  g.close()


################
################

class normalmode(object):

  def __init__(self, freqout):

    f = open(freqout) ## freq output
    f1 = f.readlines()
    f.close()

    ## collect atomic masses ##
    mass = [] ## a list for atomic masses
    for i in range(len(f1)):
      if 'Atomic Masses' in f1[i]:
        ii = 2
        while True:
          if f1[i+ii] == '\n':
            break
          else:
            mass.append(float(f1[i+ii].strip('\n').split()[-1]))
            ii = ii + 1
      else:
        pass
    mass = np.array(mass)
    self.num = len(mass)
    self.rt_mass = np.sqrt(mass)
    self.rt_mass = np.diag(self.rt_mass)

    ## collect coordiantes ##
    self.coord = np.zeros((self.num, 3))
    self.atomnote = []
    for i in range(len(f1)):
      if ' G E O M E T R Y' in f1[i] and '***' in f1[i]:
        for j in range(self.num):
          self.coord[j][0] = float(f1[i+8+j].strip('\n').split()[2])
          self.coord[j][1] = float(f1[i+8+j].strip('\n').split()[3])
          self.coord[j][2] = float(f1[i+8+j].strip('\n').split()[4])
          self.atomnote.append(f1[i+8+j].strip('\n').split()[1])

#    for i in range(len(f1)):
#      if 'FRAGMENTS' in f1[i]:
#    	  for j in range(self.num):
#    		  self.coord[j][0] = float(f1[i+3+j].strip('\n').split()[4])
#    		  self.coord[j][1] = float(f1[i+3+j].strip('\n').split()[5])
#    		  self.coord[j][2] = float(f1[i+3+j].strip('\n').split()[6])
#    		  self.atomnote.append(f1[i+3+j].strip('\n').split()[1])
#      else:
#        pass
#
    ## collect frequencies and  cartesian coordinate changes (not mass-weighted) ##
    self.freq = [] ## a list for frequencies
    for i in range(len(f1)):
      if 'Vibrations and Normal Modes' in f1[i]:
        ii = 7
        while True:
          if f1[ii+i] == '\n':
            break
          else:
            for iii in range(len(f1[i+ii].strip('\n').split())):
              self.freq.append(float(f1[i+ii].strip('\n').split()[iii]))
            ii = ii + 1 + self.num + 3
      else:
        pass
    self.freq = np.array(self.freq)
    self.num_freq = len(self.freq) ## the number of frequencies or normal modes
    self.normalmode = np.zeros((self.num_freq, self.num, 3))
    for i in range(len(f1)):
      if 'Vibrations and Normal Modes' in f1[i]:
        ii = 9
        j = 0
        while True:
          if 'List of All Frequencies' in f1[i+ii]:
            break
          else:
            if (len(f1[i+ii].strip('\n').split()) - 1) == 9:
              for jj in range(3):
                for iii in range(self.num):
                  self.normalmode[j+jj][iii][0] = float(f1[i+ii+iii].strip('\n').split()[3*jj+1])
                for iii in range(self.num):
                  self.normalmode[j+jj][iii][1] = float(f1[i+ii+iii].strip('\n').split()[3*jj+2])
                for iii in range(self.num):
                  self.normalmode[j+jj][iii][2] = float(f1[i+ii+iii].strip('\n').split()[3*jj+3])
              j = j + 3
            elif (len(f1[i+ii].strip('\n').split()) - 1) == 6:
              for jj in range(2):
                for iii in range(self.num):
                  self.normalmode[j+jj][iii][0] = float(f1[i+ii+iii].strip('\n').split()[3*jj+1])
                for iii in range(self.num):
                  self.normalmode[j+jj][iii][1] = float(f1[i+ii+iii].strip('\n').split()[3*jj+2])
                for iii in range(self.num):
                  self.normalmode[j+jj][iii][2] = float(f1[i+ii+iii].strip('\n').split()[3*jj+3])
              j = j + 2
  
            elif (len(f1[i+ii].strip('\n').split()) - 1) == 3:
              for jj in range(1):
                for iii in range(self.num):
                  self.normalmode[j+jj][iii][0] = float(f1[i+ii+iii].strip('\n').split()[3*jj+1])
                for iii in range(self.num):
                  self.normalmode[j+jj][iii][1] = float(f1[i+ii+iii].strip('\n').split()[3*jj+2])
                for iii in range(self.num):
                  self.normalmode[j+jj][iii][2] = float(f1[i+ii+iii].strip('\n').split()[3*jj+3])
              j = j + 1
            else:
              print('err about locating normal modes block')
          ii = ii + self.num + 4
      else:
        pass

    self.mw_normalmode = np.zeros((self.num_freq, self.num, 3)) ##  mass-weighted normal modes
    for i in range(self.num_freq):
      self.mw_normalmode[i] = np.matmul(self.rt_mass, self.normalmode[i])


atomnote2mass = {
'H': 1.00797,
'C': 12.011,
'N': 14.0067,
'O': 15.9994,
'Mg': 24.305,
'S': 32.06,
'Ag': 107.868,
'Au': 196.9665}


class normalmode_periodic():
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

    self.coord = []
    self.atomnote = []
    self.freq = []
    for i in range(1, len(f1)-1):
      if 'Geometry' in f1[i] and '---' in f1[i-1] and '---' in f1[i+1]:
        ii = 4
        while True:
          self.coord.append(list(map(float, f1[i+ii].strip('\n').split()[2:])))
          self.atomnote.append(f1[i+ii].strip('\n').split()[1])
          ii += 1
          if len(f1[i+ii]) < 3:
            break
      if 'Index   Frequency (cm-1)' in f1[i]:
        ii = 1
        while True:
          if len(f1[i+ii]) < 2:
            break
          self.freq.append(float(f1[i+ii].strip('\n').split()[-1]))
          ii += 1
    self.coord = np.array(self.coord)
    self.num = self.coord.shape[0]
    self.normalmode = np.zeros((len(self.freq), self.num, 3))
    ii = 0
    for i in range(len(f1)):
      if 'Displacements (x/y/z)' in f1[i]: ## not mass weighted, in Bohr
        for j in range(self.num):
          self.normalmode[ii][j][0] = float(f1[i+1+j].strip('\n').split()[-3])
          self.normalmode[ii][j][1] = float(f1[i+1+j].strip('\n').split()[-2])
          self.normalmode[ii][j][2] = float(f1[i+1+j].strip('\n').split()[-1])
        ii += 1


    self.mass = np.zeros(self.num)
    for i in range(len(self.atomnote)):
      self.mass[i] = atomnote2mass[self.atomnote[i]]

    self.rt_mass = np.diag(np.sqrt(self.mass))
    self.mw_normalmode = np.zeros(self.normalmode.shape)

    for i in range(self.normalmode.shape[0]):
      self.mw_normalmode[i] = np.matmul(self.rt_mass, self.normalmode[i])




