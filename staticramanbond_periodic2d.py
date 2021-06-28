import pickle
import numpy as np
import math
from decimal import *
from numpy.linalg import norm
from cmath import polar
import sys
from ramanbond import *

class polbond_periodic2d:
  def __init__(self, p, m, p1, p2):
    self.indcharge = (p.charge - m.charge) / (0.001 * 2)
    self.pol = (p.dipole - m.dipole) / (0.001 * 2)
    self.vec1 = p.vec1
    self.vec2 = p.vec2
    self.coord = p.coord
    self.atomnote = p.atomnote

    self.anglecos = calc_cosmatrix(self.coord)[2]

    self.polatom = {}
    for i in range(p.num):
      self.polatom.update({ i : (p.atomdip[i][0] - m.atomdip[i][0]) / (0.001 * 2) })


    ## store infor temporarily needed to calculate L matrix in before_L
    before_L = np.zeros((p.num, p.num))
    getcontext().prec = 28	
    for i in range(p.num):
      for j in range(p.num):
        before_L[i][j] = ( (1 / (2.0 * pf(p1, p2, p.atomnote[i], p.atomnote[j], calc_dis_periodic2d(p.coord[i], p.coord[j], p.vec1, p.vec2)[0], self.anglecos[i][j] ))))
  
    L = np.zeros((p.num, p.num))
    for i in range(p.num):
      for j in range(p.num):
        if j == i:
          L[i][j] = Decimal(-1*sum(before_L[i])+ 1.0) # add an arbitrary const C = 1.0
        elif j != i:
          L[i][j] = Decimal(before_L[i][j] + 1.0) # add an arbitrary const C 
    ############################################################################################
    ## Matrix equation: L Lambda = (induced charge in one direction)                           #
    ## for induced charge in x direction, we will have a set of lambdas as solution in lam_x   #
    ## similary, we have lam_y, lam_z                                                          #
    ## charge transfer Qij will be calculated based on lam_x, lam_y and lam_z                  #
    ############################################################################################
  
    lam = np.linalg.solve(L, self.indcharge)
  
    self.Q  =  np.zeros((p.num, p.num))

    for i in range(p.num):
      self.Q[i][i] = (-1.0 / (2.0 * pf(p1, p2, p.atomnote[i], p.atomnote[i], calc_dis_periodic2d(p.coord[i], p.coord[i], p.vec1, p.vec2)[0], self.anglecos[i][j] ))) * (lam[i])
  
    for i in range(p.num):
      for j in range(p.num):
        if j != i:
          self.Q[i][j] = (-1.0 / (2.0 * pf(p1, p2, p.atomnote[i], p.atomnote[j], calc_dis_periodic2d(p.coord[i], p.coord[j], p.vec1, p.vec2)[0], self.anglecos[i][j] ))) * (lam[i] - lam[j])


      self.polbond = {}
      for i in range(p.num):
        for j in range(i+1):
          ii = calc_dis_periodic2d(p.coord[i], p.coord[j], p.vec1, p.vec2)[1]
          self.polbond.update( {(i, j) : self.Q[i][j] * (p.coord[i][2] - (p.coord[j] + ii[0]*p.vec1 + ii[1]*p.vec2)[2])  } )
  


class collect_bandspt_periodic2d:
  def __init__(self, filename, closeshell):

    f = open(filename)
    f1 = f.readlines()
    f.close()
    for i in range(len(f1)):
      if 'Total nr. of atoms' in f1[i]:
        self.num = int(f1[i].strip('\n').split()[-1])
  
  
    self.atomdip = np.zeros((self.num, 1))

    self.charge  = np.zeros((self.num, 1))
    for i in range(len(f1)):

      if 'Deformation charges with respect to neutral atoms' in f1[i]:
        for ii in range(self.num):
          self.charge[ii][0] = float(f1[i+5+ii].strip('\n').split()[2])

  
      if 'Hirshfeld atomic dipole z' in f1[i]:
        if closeshell == True:
          for ii in range(self.num):
            self.atomdip[ii][0] = (float(f1[i+1+ii].strip('\n').split()[-1]))
        elif closeshell == False:
          for ii in range(self.num):
            self.atomdip[ii][0] = (float(f1[i+1+ii].strip('\n').split()[-2])) + (float(f1[i+1+ii].strip('\n').split()[-1]))
        else:
          print('error', closeshell)

      if 'D I P O L E' in f1[i]:
        self.dipole = float(f1[i+6].strip('\n').split()[-2])


    self.coord =  np.zeros((self.num, 3))
    self.atomnote = []
    for i in range(len(f1)):
      if 'Geometry' in f1[i]:
        for ii in range(self.num):
          self.coord[ii][0] = float(f1[i+4+ii].strip('\n').split()[2])
          self.coord[ii][1] = float(f1[i+4+ii].strip('\n').split()[3])
          self.coord[ii][2] = float(f1[i+4+ii].strip('\n').split()[4])
          self.atomnote.append(     f1[i+4+ii].strip('\n').split()[1])

        self.vec1 = np.array([float(f1[i+6+self.num].strip('\n').split()[1]),
                              float(f1[i+6+self.num].strip('\n').split()[2]),
                              float(f1[i+6+self.num].strip('\n').split()[3])])
        self.vec2 = np.array([float(f1[i+7+self.num].strip('\n').split()[1]),
                              float(f1[i+7+self.num].strip('\n').split()[2]),
                              float(f1[i+7+self.num].strip('\n').split()[3])])
        break

    self.coord = A2B * self.coord
    self.vec1 =  A2B * self.vec1
    self.vec2 =  A2B * self.vec2




class ramanbond_periodic2d:
  def __init__(self, vibp, vibm, stepsize):
    self.polder = (vibp.pol - vibm.pol) / (stepsize * 2)
    self.coord = 0.5 * (vibp.coord + vibm.coord) / A2B
    self.atomnote = vibp.atomnote
    self.num = len(self.atomnote)
    self.vec1 = vibp.vec1 / A2B
    self.vec2 = vibp.vec2 / A2B

    self.ramanbond = twopoint_numdif(vibp.polbond, vibm.polbond, stepsize)
    self.ramanatom = twopoint_numdif(vibp.polatom, vibm.polatom, stepsize)

    self.superatomnote = self.atomnote * 9

    self.supercoord = np.concatenate((self.coord + -1*self.vec1 + -1*self.vec2,
                                      self.coord +  0*self.vec1 + -1*self.vec2,
                                      self.coord +  1*self.vec1 + -1*self.vec2,
                                      self.coord + -1*self.vec1 +  0*self.vec2, 
                                      self.coord +  0*self.vec1 +  0*self.vec2,
                                      self.coord +  1*self.vec1 +  0*self.vec2,
                                      self.coord + -1*self.vec1 +  1*self.vec2,
                                      self.coord +  0*self.vec1 +  1*self.vec2,
                                      self.coord +  1*self.vec1 +  1*self.vec2), axis=0)

    self.superramanatom = {}
    for key in self.ramanatom:
      self.superramanatom.update({key + 4 * self.num : self.ramanatom[key]})

    supercellhash = {(-1, -1) : 0,
                     (0, -1)  : 1,
                     (1, -1)  : 2,
                     (-1, 0)  : 3,
                     (0, 0)   : 4,
                     (1, 0)   : 5,
                     (-1, 1)  : 6,
                     (0, 1)   : 7,
                     (1, 1)   : 8,}

    self.superramanbond = {}
    for key in self.ramanbond:
      ii = calc_dis_periodic2d(self.coord[key[0]], self.coord[key[1]], self.vec1, self.vec2)[1]
      self.superramanbond.update({(key[0]+4*self.num, key[1] + (supercellhash[ii])*self.num )  : self.ramanbond[key]})

    printcoord('geo_periodic2d.xyz', self.superatomnote, self.supercoord)
    dumpdata('ramanbond_periodic2d.p', self.superatomnote, self.supercoord, self.superramanatom, self.superramanbond)



