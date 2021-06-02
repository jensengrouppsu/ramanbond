import numpy as np
from math import sqrt
from math import pi
from decimal import *
from numpy.linalg import norm
from cmath import polar
import _pickle as cPickle

A2B = 1.8897261328856432 ## constant converting angstrom to bohr

rbs = {
'H': 0.31,'He': 0.28,'Li': 1.28,'Be': 0.96,'B': 0.84,'C': 0.70,'N': 0.71,'O': 0.66,'F': 0.57,'Ne': 0.58,'Na': 1.66,
'Mg': 1.41,'Al': 1.21,'Si': 1.11,'P': 1.07,'S': 1.05,'Cl': 1.02,'Ar': 1.06, 'Cu':1.38, 'Zn': 1.31, 'Ag': 1.45,'Au': 1.36
} ## covalent atomic radii, carbon's radius is tricky: sp3 0.77; sp2 0.73; sp 0.69, here use 0.70
for key in rbs:
  rbs[key] = rbs[key] * A2B
## scale the covalent atomic radii from angstrom to bohr


def twopoint_numdif(pp, mm, step):
  diff = {}
  pluskey= []
  minuskey = []
  for key in pp:
    pluskey.append(key)
  for key in mm:
    minuskey.append(key)
  for i in range(len(pluskey)):
    if pluskey[i] == minuskey[i]:
      diff.update({pluskey[i] : ( (pp.get(pluskey[i]) - mm.get(minuskey[i])) / (2.0*step) ) })
    else:
      print('atompol/bondpol in plus and minus files do not match', pluskey[i], minuskey[i])
      break
  return diff



def calc_3compo(atom, bond, atomnote, mollist, clulist):
   mol = 0.0
   inter = 0.0
   clu = 0.0
   for k in range(len(atomnote)):
     if atomnote[k] in clulist:
       clu += atom[k]
     elif atomnote[k] in mollist:
       mol += atom[k]

   for k in range(len(atomnote)):
     for l in range(k):
       if atomnote[k] in clulist and atomnote[l] in clulist:
         clu += bond[(k, l)]
       elif atomnote[k] in clulist and atomnote[l] in mollist:
         inter += bond[(k, l)]
       elif atomnote[l] in clulist and atomnote[k] in mollist:
         inter += bond[(k, l)]
       elif atomnote[k] in mollist and atomnote[l] in mollist:
         mol += bond[(k, l)]

   return mol, inter, clu




def calc_dis(coord1, coord2):
  return sqrt((coord1[0]-coord2[0])**2.0 +  (coord1[1]-coord2[1])**2.0 + (coord1[2]-coord2[2])**2.0 )


def calc_dismatrix(coord):
  num = coord.shape[0]
  dis = np.zeros((num, num))
  for i in range(num):
    for j in range(num):
      dis[i][j] = calc_dis(coord[i], coord[j])

  return dis


def calc_cosmatrix(coord):
  num = coord.shape[0]
  anglecos = np.zeros((3, num, num))

  for ii in range(3):
    for i in range(num):
      for j in range(num):
        if i == j:
          pass
        else:
          anglecos[ii][i][j] = abs(np.dot((coord[i] - coord[j]), np.identity(3)[ii]) / np.linalg.norm(coord[i]-coord[j]))

  return anglecos



def pf(p1, p2, Za, Zb, distance, anglecos): 
  Ra = rbs[Za]
  Rb = rbs[Zb]
  try:
    return math.exp(p1*(distance/(Ra+Rb))**2 + p2 * (1-anglecos))
  except OverflowError:
    return float('inf')

def solve_lagrange(p1, p2, charge, atomnote, dis, anglecos):
  num = len(charge)
  ## store infor temporarily needed to calculate L matrix in before_L
  before_L = np.zeros((num, num))
  getcontext().prec = 28
  for i in range(num):
    for j in range(num):
      if j != i:
        try:
          before_L[i][j] = 1/(2*pf(p1, p2, atomnote[i], atomnote[j], dis[i][j], anglecos[i][j]))
        except OverflowError as err:
          before_L[i][j] = 0.0
      elif j == i:
          before_L[i][j] = 0.0
      else:
          print('sth is wrong with before_L')


  L = np.zeros((num, num))
  for i in range(num):
    for j in range(num):
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
  lam = np.linalg.solve(L, charge)
  Q = np.zeros((num, num))
  for i in range(num):
    for j in range(num):
      if j == i:
        pass
      else:
        try:			
          Q[i][j] = -1* (Decimal(lam[i]) - Decimal(lam[j]))/(Decimal(2*pf(p1, p2, atomnote[i], atomnote[j], dis[i][j], anglecos[i][j])))
        except OverflowError:
          Q[i][j] = 0.0
  return Q

def calc_chargeflow(p1, p2, charge, atomnote, dis, anglecos):
  return solve_lagrange(p1, p2, charge.real, atomnote, dis, anglecos) + (1j) * (solve_lagrange(p1, p2, charge.imag, atomnote, dis, anglecos))


def sum_bondpol(num, coord_bohr, Q_x, Q_y, Q_z):
  bondpol = np.zeros((3, 3), dtype=np.complex_)
  for i in range(3):
    for j in range(num):
      for jj in range(j):
        bondpol[0][i] += Q_x[j][jj]*(coord_bohr[j][i] - coord_bohr[jj][i])
  for i in range(3):
    for j in range(num):
      for jj in range(j):
        bondpol[1][i] += Q_y[j][jj]*(coord_bohr[j][i] - coord_bohr[jj][i])
  for i in range(3):
    for j in range(num):
      for jj in range(j):
        bondpol[2][i] += Q_z[j][jj]*(coord_bohr[j][i] - coord_bohr[jj][i])
  return bondpol


def collectaoresponse(filename):
  f = open(filename)
  f1  =f.readlines()
  f.close()
  ii = 0
  aoresponse_pointer = {}
  for i in range(len(f1)):
    if 'No. of frequency' in f1[i]:
      ii += 1
      aoresponse_pointer.update({ii : i})

  if ii == 1:
    return [f1]
  elif ii > 1:
    returncontent = []
    for i in range(1, ii):
      returncontent.append(f1[:aoresponse_pointer[1]] + f1[aoresponse_pointer[i]:aoresponse_pointer[i+1]])
    returncontent.append(f1[:aoresponse_pointer[1]] + f1[aoresponse_pointer[ii]:])

    ## if the calculation is finished normally
    if 'END' in f1[-1] and 'NORMAL TERMINATION' in f1[-2]:
      return returncontent
    ## if the calculation is terminated prematurely
    else:
      return returncontent[:-1]

  else:
    print('wrong file collected: '+filename)


def orient_ave(R):
  alpha = (1.0 / 3.0) * (R[0][0] + R[1][1] + R[2][2])
  gamma = sqrt( 0.5 * (polar(R[0][0]-R[1][1])[0]**2 + polar(R[1][1]-R[2][2])[0]**2.0 + polar(R[2][2]-R[0][0])[0]**2 ) + 3*(polar(R[0][1])[0]**2 + polar(R[0][2])[0]**2 + polar(R[1][2])[0]**2 ) )

  return sqrt((45*polar(alpha)[0]**2 + 7*gamma**2) / 45.0)


def lorentzian(x, peakposition, scalefactor, fwhm):
  gamma = fwhm / 2
  ## the lorentzian distribution unscaled is normalized and centered at peakposition ##
  ## the width is controlled by gamma ##
  return  scalefactor * ( 1.0 / pi ) * ( gamma / ( ( x - peakposition )**2 + gamma**2 ) ) 

def sum_lorentzian(x, peakpositions, peakheights, fwhm):
  y = np.array([lorentzian(x, peakpositions[i], peakheights[i], fwhm) for i in range(len(peakpositions))])
  return y.sum(axis=0)

def dumpdata(picklefile, atomnote, coord, atom, bond):
  with open(picklefile, "wb") as output_file:
   cPickle.dump((atomnote, coord, atom, bond), output_file) ## atom/bond can be polatom/bond or ramanatom/bond


def printcoord(xyzname, atomnote, coord):
  num = len(atomnote)
  g = open(xyzname, 'w')
  g.write(str(num)+'\n')
  g.write('\n')
  for i in range(num):
    g.write('{0:}        {1: .8f}        {2: .8f}        {3: .8f}\n'.format(atomnote[i], coord[i][0], coord[i][1], coord[i][2]))
  g.close()



def calc_rotmatrix(A, B): 
## calc the rotation matrix which aligns two vectors A and B parallel 
  A = A / np.linalg.norm(A)
  B = B / np.linalg.norm(B)
  
  tmpcos = np.dot(A, B)
  tmpsin = np.linalg.norm(np.cross(A, B) )
  G = np.zeros((3, 3))
  G[0][0] = tmpcos
  G[0][1] = -1 * tmpsin
  G[1][0] = tmpsin
  G[1][1] = tmpcos
  G[2][2] = 1.0
  u = A.copy()
  v = (B - tmpcos*A) / (np.linalg.norm(B-tmpcos*A ))
  w = np.cross(B, A)  
  F_inv = np.zeros((3, 3))
  for i in range(3):
    F_inv[i][0] = u[i]
  
  for i in range(3):
    F_inv[i][1] = v[i]
  
  for i in range(3):
    F_inv[i][2] = w[i]
  
  F = np.linalg.inv(F_inv)
  U = np.matmul(F_inv , np.matmul(G, F))

  return U


def align_plane(coord, p1, p2, p3, direction): 
## p1 p2 p3 are indices of atoms in a plane
## direction can be x [1, 0, 0] for example which is normal to the atomic plane
## origin is the index of the atom which is selected to be at origin
  vec1 = coord[p2 -1] - coord[p1 - 1] 
  vec2 = coord[p3 -1] - coord[p1 - 1] 
  A = np.cross(vec1, vec2)
  B = np.array(direction)
  U = calc_rotmatrix(A, B)
  num = coord.shape[0]
  for i in range(num):
    coord[i] = np.matmul(U, coord[i])


def align_line(coord, v1, v2, direction):

  A1 = coord[v2-1] - coord[v1-1]
  B1 = np.array(direction)

  U1 = calc_rotmatrix(A1, B1)
  num = coord.shape[0]
  for i in range(num):
    coord[i] = np.matmul(U1, coord[i])



def coordfromxyz(self, filename):
  f = open(filename)
  f1 = f.readlines()
  f.close()
  num = int(f1[0].strip('\n'))
  coord = np.zeros((self.num, 3))
  atomnote = []
  for i in range(self.num):
    atomnote.append(f1[i+2].strip('\n').split()[0])
    coord[i][0] = float(f1[i+2].strip('\n').split()[1])
    coord[i][1] = float(f1[i+2].strip('\n').split()[2])
    coord[i][2] = float(f1[i+2].strip('\n').split()[3])
  return atomnote, coord


## for periodic systems ##

def calc_dis_periodic2d(coord1, coord2, vec1, vec2):
  L = []
  ii = []
  for i in range(-1, 2):
    for j in range(-1, 2):
      newcoord = coord2 + i * vec1 + j * vec2
      L.append( math.sqrt((newcoord[0] - coord1[0])**2 + (newcoord[1] - coord1[1])**2+ (newcoord[2] - coord1[2])**2) )
      ii.append((i, j))

  minL = L[0]
  minii = ii[0]
  for i in range(len(L)):
    if 0 < L[i] < minL:
      minL = L[i]
      minii = ii[i]

  return minL, minii


def periodiccoordfromtxt(self, filename):
  f = open(filename)
  f1 = f.readlines()
  f.close()
  coord = []
  atomnote = []
  vec = []
  for i in range(len(f1)):
    if 'vec' in f1[i].lower():
      vec.append(list(map(float, f1[i].strip('\n').split()[1:])))
    else:
      atomnote.append(f1[i].strip('\n').split()[0])
      coord.append(list(map(float, f1[i].strip('\n').split()[1:])))

  coord = np.array(self.coord)
  vec = np.array(self.vec)
  return atomnote, coord, vec

