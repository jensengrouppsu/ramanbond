import numpy as np
import math
from decimal import *
from numpy.linalg import norm
from cmath import polar

A2B = 1.8897261328856432 ## constant converting angstrom to bohr

class adfnormalmodes(object):

  def __init__(self, freqout):

    f = open(freqout) ## freq output
    f1 = f.readlines()
    f.close()

    ## collect atomic masses ##
    masses = [] ## a list for atomic masses
    for i in range(len(f1)):
      if 'Atomic Masses' in f1[i]:
        ii = 2
        while True:
          if f1[i+ii] == '\n':
            break
          else:
            masses.append(float(f1[i+ii].strip('\n').split()[-1]))
            ii = ii + 1
      else:
        pass
    masses = np.array(masses)
    self.num = len(masses)
    self.rt_masses = np.sqrt(masses)
    self.rt_masses = np.diag(self.rt_masses)

    ## collect coordiantes ##
    self.coord = np.zeros((self.num, 3))
    self.atomnote = []
    for i in range(len(f1)):
      if ' G E O M E T R Y  ***  3D  Molecule  ***' in f1[i]:
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
    self.normalmodes = np.zeros((self.num_freq, self.num, 3))
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
                  self.normalmodes[j+jj][iii][0] = float(f1[i+ii+iii].strip('\n').split()[3*jj+1])
                for iii in range(self.num):
                  self.normalmodes[j+jj][iii][1] = float(f1[i+ii+iii].strip('\n').split()[3*jj+2])
                for iii in range(self.num):
                  self.normalmodes[j+jj][iii][2] = float(f1[i+ii+iii].strip('\n').split()[3*jj+3])
              j = j + 3
            elif (len(f1[i+ii].strip('\n').split()) - 1) == 6:
              for jj in range(2):
                for iii in range(self.num):
                  self.normalmodes[j+jj][iii][0] = float(f1[i+ii+iii].strip('\n').split()[3*jj+1])
                for iii in range(self.num):
                  self.normalmodes[j+jj][iii][1] = float(f1[i+ii+iii].strip('\n').split()[3*jj+2])
                for iii in range(self.num):
                  self.normalmodes[j+jj][iii][2] = float(f1[i+ii+iii].strip('\n').split()[3*jj+3])
              j = j + 2
  
            elif (len(f1[i+ii].strip('\n').split()) - 1) == 3:
              for jj in range(1):
                for iii in range(self.num):
                  self.normalmodes[j+jj][iii][0] = float(f1[i+ii+iii].strip('\n').split()[3*jj+1])
                for iii in range(self.num):
                  self.normalmodes[j+jj][iii][1] = float(f1[i+ii+iii].strip('\n').split()[3*jj+2])
                for iii in range(self.num):
                  self.normalmodes[j+jj][iii][2] = float(f1[i+ii+iii].strip('\n').split()[3*jj+3])
              j = j + 1
            else:
              print('err about locating normal modes block')
          ii = ii + self.num + 4
      else:
        pass

    self.mw_normalmodes = np.zeros((self.num_freq, self.num, 3)) ##  mass-weighted normal modes
    for i in range(self.num_freq):
      self.mw_normalmodes[i] = np.matmul(self.rt_masses, self.normalmodes[i])

  def calc_stepsize(self, mode):
    for i in range(self.num_freq):
      if abs(float(mode) - self.freq[i]) < 0.02:
        self.stepsize = (norm(self.mw_normalmodes[i]) / norm(self.normalmodes[i]) * 0.01 ) ##cart_stepsize = 0.01 usually 
## normalmodes due to digital precision are not strictly normalzied 
      else:
        pass

  def calc_displacedcoord(self, mode, num_step):
    ## num_step = 1 usually
    for i in range(self.num_freq):
      if abs(float(mode) - self.freq[i]) < 0.02:
        self.displacedcoord_plus = self.coord + (self.normalmodes[i] / ( A2B * norm(self.normalmodes[i] ))) * 0.01 * num_step ## convert from Bohr to Angstrom !
        self.displacedcoord_minus = self.coord - (self.normalmodes[i] / ( A2B * norm(self.normalmodes[i]))) * 0.01 * num_step

  def create_modefile(self, mode, num_step, template):
    self.calc_displacedcoord(mode, num_step)
    f = open(template)
    f1 = f.readlines()
    f.close()
    
    g = open('mode'+'{0:.3f}'.format(mode)+'-p.run', 'w')
    for i in range(len(f1)):
      if 'atoms' in f1[i].strip('\n').lower():
        g.write(f1[i])
        for j in range(self.num):
          g.write('{0:<3}    {1: .8f}    {2: .8f}    {3: .8f}\n'.format(self.atomnote[j], self.displacedcoord_plus[j][0], self.displacedcoord_plus[j][1], self.displacedcoord_plus[j][2]))
      else:
        g.write(f1[i])
    g.close()

    g = open('mode'+'{0:.3f}'.format(mode)+'-m.run', 'w')
    for i in range(len(f1)):
      if 'atoms' in f1[i].strip('\n').lower():
        g.write(f1[i])
        for j in range(self.num):
          g.write('{0:<3}    {1: .8f}    {2: .8f}    {3: .8f}\n'.format(self.atomnote[j], self.displacedcoord_minus[j][0], self.displacedcoord_minus[j][1], self.displacedcoord_minus[j][2]))
      else:
        g.write(f1[i])
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

