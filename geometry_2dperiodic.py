import numpy as np

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







