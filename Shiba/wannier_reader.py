from Shiba.basic_functions import logger,cell
from dataclasses import dataclass, field
import numpy as np
import os
@dataclass
class WannierHamiltonian:
      norb: int = 0
      norbup: int = 0
      norbdn: int = 0
      degs: list = field(default_factory=list)
      nws: int = 0
      calculation_mode: int = 0
      mels: np.ndarray = np.array([])
      melsup: np.ndarray = np.array([])
      melsdn: np.ndarray = np.array([])
      blocks: int = 0

      def mode_detector(self,Parameters):
          if Parameters.up_h == '' and Parameters.down_h == '' and Parameters.noncolin_h != '':
             self.calculation_mode = 1
          elif Parameters.up_h != '' and Parameters.down_h != '' and Parameters.noncolin_h == '':
             self.calculation_mode = 2
          else:
             logger('\nIncompatible number of Hamiltonian files, exiting.' + '\n')

      def read_wannier_info(self,Parameters):
          self.mode_detector(Parameters)
          if(self.calculation_mode == 1):
            hfile = os.path.join(Parameters.input_dir,Parameters.noncolin_h)
            with open(hfile) as f:
                  for row, line in enumerate(f):
                      if(row==0):
                         continue
                      elif(row==1):
                         self.norb = [int(number) for number in line.split()][0]
                      elif(row==2):
                         self.nws = [int(number) for number in line.split()][0]
                      elif(row<=self.nws//15+3):
                         self.degs.extend([float(number) for number in line.split()])
                      else:
                         break
            self.mels = np.genfromtxt(hfile, skip_header=self.nws//15+4)
            logger('\nWannier orbitals:\t' + str(self.norb) + '\n')
            logger('Wigner-Seitz cells:\t' + str(self.nws) + '\n')
            invdegs = np.zeros_like(self.degs)
            for i in range(len(self.degs)):
                invdegs[i] = 1.0/self.degs[i]
            logger('(Reciprocal grid:\t' + str(int(np.sum(invdegs))) + ')' + '\n')
            logger('Depth of WS overlaps:\t' + str(Parameters.nc) + '\n')
            for idx in range(self.nws):
                if(cell(self.mels[idx*self.norb**2,0],self.mels[idx*self.norb**2,1],self.mels[idx*self.norb**2,2])<=Parameters.nc):
                  self.blocks+=1
            logger('Full Hamiltonian size:\t' + str((self.blocks*4*(self.norb//2))) + 'x' + str((self.blocks*4*(self.norb//2))) + '\n\n')
          elif(self.calculation_mode == 2):
            wfileup = os.path.join(Parameters.input_dir,Parameters.up_h)
            wfiledn = os.path.join(Parameters.input_dir,Parameters.down_h)
            with open(wfiledn) as f:
                for row, line in enumerate(f):
                    if(row==0):
                       continue
                    elif(row==1):
                       self.norbdn = [int(number) for number in line.split()][0]
                    elif(row==2):
                       self.nws = [int(number) for number in line.split()][0]
                    elif(row<=self.nws//15+3):
                       self.degs.extend([float(number) for number in line.split()])
                    else:
                       break
            self.melsdn = np.genfromtxt(wfiledn, skip_header=self.nws//15+4)
            invdegs = np.zeros_like(self.degs)
            for i in range(len(self.degs)):
                invdegs[i] = 1.0/self.degs[i]
            with open(wfileup) as f:
                 for row, line in enumerate(f):
                    if(row==0):
                        continue
                    elif(row==1):
                        self.norbup = [int(number) for number in line.split()][0]
                    else:
                        break
            self.melsup = np.genfromtxt(wfileup, skip_header=self.nws//15+4)
            logger('\nWannier orbitals:\t' + str(self.norbdn)+"+"+str(self.norbup)+"="+str(self.norbdn+self.norbup) + '\n')
            logger('Wigner-Seitz cells:\t' + str(self.nws) + '\n')
            logger('(Reciprocal grid:\t' + str(int(np.sum(invdegs))) + ')' + '\n')
            logger('Depth of WS overlaps:\t' + str(Parameters.nc) + '\n')
            self.norb=self.norbdn+self.norbup
            for idx in range(self.nws):
                if(cell(self.melsdn[idx*self.norbdn**2,0],self.melsdn[idx*self.norbdn**2,1],self.melsdn[idx*self.norbdn**2,2])<=Parameters.nc):
                   self.blocks+=1
            logger('Full Hamiltonian size:\t' + str((self.blocks*4*(self.norb//2))) + 'x' + str((self.blocks*4*(self.norb//2))) + '\n\n')