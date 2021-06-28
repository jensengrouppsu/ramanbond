import numpy as np
from math import cos
from math import sin
from math import pi
from math import sqrt

class periodic_3d():
  def __init__(self, elem, aa, bb, cc, alpha, beta, gamma):
    ## generate the lattice vectors given the crystal parameters, distance in A and angle in degree
    self.elem = elem
    alpha = alpha * pi / 180.0
    beta = beta * pi / 180.0
    gamma = gamma * pi / 180.0
  
    self.origin = np.array([0, 0, 0])
    self.vec1 = np.array([cc, 0, 0])
    self.vec2 = np.array([bb*cos(alpha), bb*sin(alpha), 0 ])
    self.vec3 = np.array([aa*cos(beta), 
                          aa*(cos(gamma)-cos(beta)*cos(alpha)) / sin(alpha), 
                          aa*sqrt(1-cos(beta)**2 - ( (cos(gamma)-cos(beta) * cos(alpha)) / sin(alpha) )**2)])

    self.V = aa*bb*cc*sqrt(1+2*cos(alpha)*cos(beta)*cos(gamma)-cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2)
    self.S = bb*cc*sin(alpha)
    self.h = self.V / self.S


  def supercell(self, irange, jrange, krange):
    ## generate the coord of atoms in the super cell given the lattice vectors
    self.supercellatomnote = []
    self.supercellcoord = []
    for ii in range(len(self.elem)):
      for i in range(irange):
        for j in range(jrange):
          for k in range(krange):
            self.supercellatomnote.append(self.elem[ii])
            self.supercellcoord.append(
                       [(self.origin + i*self.vec1 + j*self.vec2 + k*self.vec3)[0],
                        (self.origin + i*self.vec1 + j*self.vec2 + k*self.vec3)[1],
                        (self.origin + i*self.vec1 + j*self.vec2 + k*self.vec3)[2]])

    self.supercellcoord = np.array(self.supercellcoord)


## define lattice parameters ##
### Ag ##
#aa = 2.942 ## in A
#bb = 2.942 
#cc = 2.942
#alpha = 60 ## in degree
#beta = 60
#gamma = 60
#########

### Cu ##
#aa = 2.561 ## in A
#bb = 2.561 
#cc = 2.561
#alpha = 60 ## in degree
#beta = 60
#gamma = 60
#elem = 'Cu'
#########


class periodic_2d:
  def __init__(self, filename):
  ## read the lattice vectors, atomnotes and coord of atoms in the unit cell of a 2d periodic system
    f = open(filename)
    f1 = f.readlines()
    f.close()
    self.atomnote = []
    self.coord = []
    
    for i in range(len(f1)):
      if 'VEC1' in f1[i]:
        self.vec1 = np.array([float(f1[i].strip('\n').split()[1]), float(f1[i].strip('\n').split()[2]), float(f1[i].strip('\n').split()[3])])
      elif 'VEC2' in f1[i]:
        self.vec2 = np.array([float(f1[i].strip('\n').split()[1]), float(f1[i].strip('\n').split()[2]), float(f1[i].strip('\n').split()[3])])
      else:
        self.atomnote.append(f1[i].strip('\n').split()[0])
        self.coord.append([float(f1[i].strip('\n').split()[1]), float(f1[i].strip('\n').split()[2]), float(f1[i].strip('\n').split()[3])])
    self.coord = np.array(self.coord)


  def supercell(self, irange, jrange):
    ## generate coord for atoms in the supercell given the lattice vectors
    self.supercellatomnote = []
    self.supercellcoord = []
    for i in irange:
     for j in jrange:
       for k in range(len(self.atomnote)):
         self.supercellatomnote.append(self.atomnote[k])
         self.supercellcoord.append(self.coord[k]+i*self.vec1+j*self.vec2)
    self.supercellcoord = np.array(self.supercellcoord)







