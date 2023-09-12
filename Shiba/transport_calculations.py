from dataclasses import dataclass, field
from Shiba.basic_functions import cell, logger
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy import linalg, special, sparse, constants

@dataclass
class TransportCalculation:
      h: np.ndarray = np.array([])
      Id: np.ndarray = np.array([])
      Idnorb2: np.ndarray = np.array([])
      Idnorbup: np.ndarray = np.array([])
      Idnorbdn: np.ndarray = np.array([])
      heff: list = field(default_factory=list)
      Gamma1: np.ndarray = np.array([])
      Gamma2: np.ndarray = np.array([])
      Gamma1RR: np.ndarray = np.array([])
      Gamma2LL: np.ndarray = np.array([])
      GammaProd: np.ndarray = np.array([])
      voltage: np.ndarray = np.array([])
      current: np.ndarray = np.array([])
      conductance: np.ndarray = np.array([])
      eps: np.ndarray = np.array([])
      w: np.ndarray = np.array([])
      transmission: np.ndarray = np.array([])
      spectral: np.ndarray = np.array([])

      def construct_nambu(self,Parameters,Wannier_h):
          self.h = np.zeros((Wannier_h.blocks*4*(Wannier_h.norb//2),
                  Wannier_h.blocks*4*(Wannier_h.norb//2)), dtype=complex)
          self.Id = np.eye(Wannier_h.blocks*4*(Wannier_h.norb//2), dtype=complex)
          self.Idnorb2 = np.eye(Wannier_h.norb//2, dtype=complex)
          if(Wannier_h.calculation_mode==2):
             self.Idnorbdn = np.eye(Wannier_h.norbdn, dtype=complex)
             self.Idnorbup = np.eye(Wannier_h.norbup, dtype=complex)
      
          # Loop over WS cells
          ibl=-1
          for idx in range(Wannier_h.nws):
              # Take only the WS overlaps with depth up to NC
              if(Wannier_h.calculation_mode==1):
                if(cell(Wannier_h.mels[idx*Wannier_h.norb**2,0],Wannier_h.mels[idx*Wannier_h.norb**2,1],Wannier_h.mels[idx*Wannier_h.norb**2,2])<=Parameters.nc):
                    ibl+=1
                    # For each WS overlap, extract different spin components
                    huu = np.zeros((Wannier_h.norb//2,Wannier_h.norb//2), dtype=complex)
                    hdd = np.zeros((Wannier_h.norb//2,Wannier_h.norb//2), dtype=complex)
                    hud = np.zeros((Wannier_h.norb//2,Wannier_h.norb//2), dtype=complex)
                    hdu = np.zeros((Wannier_h.norb//2,Wannier_h.norb//2), dtype=complex)

                    ie=-1
                    io=-1
                    for i in range(Wannier_h.norb):
                       if(i%2==0):
                          ie+=1
                       else:
                          io+=1
                       je=-1
                       jo=-1
                       for j in range(Wannier_h.norb):
                           if(j%2==0):
                              je+=1
                           else:
                              jo+=1

                           if(i%2==0 and j%2==0):
                              huu[ie,je]=Wannier_h.mels[idx*Wannier_h.norb**2+i*Wannier_h.norb+j,5]+1.0j*Wannier_h.mels[idx*Wannier_h.norb**2+i*Wannier_h.norb+j,6]
                           if(i%2==0 and j%2!=0):
                              hud[ie,jo]=Wannier_h.mels[idx*Wannier_h.norb**2+i*Wannier_h.norb+j,5]+1.0j*Wannier_h.mels[idx*Wannier_h.norb**2+i*Wannier_h.norb+j,6]
                           if(i%2!=0 and j%2==0):
                              hdu[io,je]=Wannier_h.mels[idx*Wannier_h.norb**2+i*Wannier_h.norb+j,5]+1.0j*Wannier_h.mels[idx*Wannier_h.norb**2+i*Wannier_h.norb+j,6]
                           if(i%2!=0 and j%2!=0):  
                              hdd[io,jo]=Wannier_h.mels[idx*Wannier_h.norb**2+i*Wannier_h.norb+j,5]+1.0j*Wannier_h.mels[idx*Wannier_h.norb**2+i*Wannier_h.norb+j,6]

                    # Align Fermi level at zero
                    huu -= Parameters.MU*Idnorb2
                    hdd -= Parameters.MU*Idnorb2

                    # Construct the corresponding blocks of the full Nambu Hamiltonian
                    row = ibl*4*(Wannier_h.norb//2) # Block diagonal super matrix
                    col = ibl*4*(Wannier_h.norb//2)
                    for i in range(0,4*(Wannier_h.norb//2),4):
                        for j in range(0,4*(Wannier_h.norb//2),4):
                            self.h[row+i+0,col+j+0] = huu[i//4,j//4]
                            self.h[row+i+0,col+j+1] = -Parameters.delta*(i==j)
                            self.h[row+i+0,col+j+2] = hud[i//4,j//4]
                            self.h[row+i+0,col+j+3] = 0.0
                            self.h[row+i+1,col+j+0] = -Parameters.delta*(i==j)
                            self.h[row+i+1,col+j+1] = -hdd[i//4,j//4]
                            self.h[row+i+1,col+j+2] = 0.0
                            self.h[row+i+1,col+j+3] = -hdu[i//4,j//4]
                            self.h[row+i+2,col+j+0] = hdu[i//4,j//4]
                            self.h[row+i+2,col+j+1] = 0.0
                            self.h[row+i+2,col+j+2] = hdd[i//4,j//4]
                            self.h[row+i+2,col+j+3] = Parameters.delta*(i==j)
                            self.h[row+i+3,col+j+0] = 0.0
                            self.h[row+i+3,col+j+1] = -hud[i//4,j//4]
                            self.h[row+i+3,col+j+2] = Parameters.delta*(i==j)
                            self.h[row+i+3,col+j+3] = -huu[i//4,j//4]
          
              else:
                  if(cell(Wannier_h.melsdn[idx*Wannier_h.norbdn**2,0],Wannier_h.melsdn[idx*Wannier_h.norbdn**2,1],Wannier_h.melsdn[idx*Wannier_h.norbdn**2,2])<=Parameters.nc):
                     ibl+=1
                     # For each WS overlap, extract different spin components
                     if(Wannier_h.norbup>=Wannier_h.norbdn):
                        huu = np.zeros((Wannier_h.norbup,Wannier_h.norbup), dtype=complex)
                        hdd = np.zeros((Wannier_h.norbup,Wannier_h.norbup), dtype=complex)
                     else:
                        huu = np.zeros((Wannier_h.norbdn,Wannier_h.norbdn), dtype=complex)
                        hdd = np.zeros((Wannier_h.norbdn,Wannier_h.norbdn), dtype=complex)

                     for i in range(Wannier_h.norbdn):
                         for j in range(Wannier_h.norbdn):
                             hdd[i,j]=Wannier_h.melsdn[idx*Wannier_h.norbdn**2+i*Wannier_h.norbdn+j,5]+1.0j*Wannier_h.melsdn[idx*Wannier_h.norbdn**2+i*Wannier_h.norbdn+j,6]

                     for i in range(Wannier_h.norbup):
                         for j in range(Wannier_h.norbup):
                             huu[i,j]=Wannier_h.melsup[idx*Wannier_h.norbup**2+i*Wannier_h.norbup+j,5]+1.0j*Wannier_h.melsup[idx*Wannier_h.norbup**2+i*Wannier_h.norbup+j,6]

                     # Align Fermi level at zero
                     hdd -= Parameters.MU*self.Idnorbdn
                     huu -= Parameters.MU*self.Idnorbup

                     # Construct the corresponding blocks of the full Nambu Hamiltonian
                     row = ibl*4*(Wannier_h.norb//2) # Block diagonal super matrix
                     col = ibl*4*(Wannier_h.norb//2)
                     for i in range(0,4*(Wannier_h.norb//2),4):
                         for j in range(0,4*(Wannier_h.norb//2),4):
                             self.h[row+i+0,col+j+0] = huu[i//4,j//4]
                             self.h[row+i+0,col+j+1] = -Parameters.delta*(i==j)
                             self.h[row+i+0,col+j+2] = Parameters.alpha
                             self.h[row+i+0,col+j+3] = 0.0
                             self.h[row+i+1,col+j+0] = -Parameters.delta*(i==j)
                             self.h[row+i+1,col+j+1] = -hdd[i//4,j//4]
                             self.h[row+i+1,col+j+2] = 0.0
                             self.h[row+i+1,col+j+3] = -Parameters.alpha
                             self.h[row+i+2,col+j+0] = Parameters.alpha
                             self.h[row+i+2,col+j+1] = 0.0
                             self.h[row+i+2,col+j+2] = hdd[i//4,j//4]
                             self.h[row+i+2,col+j+3] = Parameters.delta*(i==j)
                             self.h[row+i+3,col+j+0] = 0.0
                             self.h[row+i+3,col+j+1] = -Parameters.alpha
                             self.h[row+i+3,col+j+2] = Parameters.delta*(i==j)
                             self.h[row+i+3,col+j+3] = -huu[i//4,j//4]

          #        if(wannier_h.calculation_mode==1):
          #           del mels
          #        elif(wannier_h.calculation_mode==2):
          #           del melsup, melsdn
          #              del degs, invdegs

          # The Wannier matrix-element data may not always be symmetric, so
          # the Hamiltonian is symmetrized to assure proper analytic structure
          self.h = 0.5*(self.h+self.h.conj().T)
          logger("The full Hamiltonian is hermitian? " + str(np.all(np.abs(self.h-self.h.conj().T) < 1.0e-12)) + '\n')

      def construct_coupling_matrices(self,Parameters,Wannier_h):
          Id4     = np.eye(4, dtype=complex)
          gam1    = np.zeros((4*(Wannier_h.norb//2),4*(Wannier_h.norb//2)), dtype=complex)
          gam2    = np.zeros((4*(Wannier_h.norb//2),4*(Wannier_h.norb//2)), dtype=complex)

          # Each orbital from NSTM is coupled to the STM orbital (full)
          for i in Parameters.nstm:
              for j in Parameters.nstm:
                  gam1[4*(i-1):4*i,4*(j-1):4*j] = Parameters.gamma*Id4

          # Each orbital from NSUB is coupled to individual substrate orbitals (diagonal)
          for i in Parameters.nsub:
              gam2[4*(i-1):4*i,4*(i-1):4*i] = Parameters.frac*Parameters.gamma*Id4

          # Each WS overlap block has the same coupling matrix structure,
          # so we duplicate the individual blocks as kronecker products
          blockstructure=np.eye(Wannier_h.blocks, dtype=complex) # Block diagonal super matrix
          self.Gamma1 = np.kron(blockstructure,gam1)
          self.Gamma2 = np.kron(blockstructure,gam2)
          #del gam1, gam2, blockstructure

      def compute_spectral_function(self,Parameters):
          self.w = np.linspace(-Parameters.vran, Parameters.vran, Parameters.vpts)
          self.spectral = np.zeros_like(self.w)
          self.transmission = np.zeros_like(self.w)
          for i in range(len(self.w)):
              # Spectral function is calculated as the imaginary part of the retarded Green function
              # Isolated system with artificial broadening (i*eta)
              # spectral[i] = (-1.0/np.pi)*np.imag(np.trace(np.linalg.inv((w[i]+1.0e-12j)*Id - h)))
              # Coupled system with physical broadening (i*Gamma)
              A = np.linalg.inv(self.w[i]*self.Id - self.heff)
              self.spectral[i] = (-1.0/np.pi)*np.imag(np.trace(A))
              # Transmission function
              self.transmission[i] = np.real(np.trace(self.Gamma1 @ A @ self.Gamma2 @ A.conj().T))

          np.savetxt(os.path.join(Parameters.outdir,'specdata.out'), np.array([self.w,self.spectral,self.transmission]).T, fmt=['%.8e', '%.8e', '%.8e'])
          #del h, heff, Gamma1, Gamma2
   
      def solve_hamiltonian(self,Parameters):
          # Effective Hamiltonian (coupled system)
          self.heff = self.h - 0.5j*(self.Gamma1+self.Gamma2)

          # Solve non-hermitian eigensystem (separate left/right eigenvectors)
          logger('\n')
          logger('Diagonalization ...\n')
          self.eps, psil, psir = linalg.eig(self.heff, left=True, right=True)
          overlap_lr = psil.conj().T @ psir
          invover_lr = np.linalg.inv(overlap_lr)
          overlap_rl = psir.conj().T @ psil
          invover_rl = np.linalg.inv(overlap_rl)
          self.Gamma1RR = psir.conj().T @ self.Gamma1 @ psir @ invover_lr
          self. Gamma2LL = psil.conj().T @ self.Gamma2 @ psil @ invover_rl
          if(Parameters.opt==1):
             self.GammaProd = sparse.coo_matrix(self.Gamma1RR*self.Gamma2LL.T)
          #del psil, psir, overlap_lr, invover_lr, overlap_rl, invover_rl

          if(Parameters.spec):
             logger('\n')
             logger('Computing spectral function ...\n')
             self.compute_spectral_function(Parameters)

          logger('\n')
          logger('Current-voltage ...\n')
          self.voltage = np.linspace(-Parameters.vran, Parameters.vran, Parameters.vpts)
          self.current = np.zeros_like(self.voltage)
          self.conductance = np.zeros_like(self.voltage)
          for i in range(len(self.voltage)):
              if(Parameters.opt==0):
                 self.current[i] = self.integralOriginal(self.voltage[i],0.0,self.eps,Parameters.BETA)
              elif(Parameters.opt==1):
                 self.current[i] = self.integralSparse(self.voltage[i],0.0,self.eps,Parameters.BETA)
              elif(Parameters.opt==2):
                 self.current[i] = self.integralContract(self.voltage[i],0.0,self.eps,Parameters.BETA)
              else:
                 logger('\nIncompatible optimization flag, exiting.' + '\n')
                 sys.exit(0) 
          # The above structure can be written perhaps more elegantly as below (for Python 3.10 onward)
          #    match opt:
          #        case 0:
          #            current[i] = integralOriginal(voltage[i],0.0)
          #        case 1:
          #            current[i] = integralSparse(voltage[i],0.0)
          #        case 2:
          #            current[i] = integralContract(voltage[i],0.0)
          #        case _:
          #            logger('\nIncompatible optimization flag, exiting.' + '\n')
          #            sys.exit(0)

          # Conductance is calculated as a derivative of the current w.r.t. voltage
          self.conductance = np.gradient(self.current, self.voltage[1]-self.voltage[0])

          np.savetxt(os.path.join(Parameters.outdir,'ivdata.out'), np.array([self.voltage,self.current,self.conductance]).T, fmt=['%.8e', '%.8e', '%.8e'])

      def integralOriginal(self,V1, V2,BETA):
          # Original result as a nested loop over eigenvalues
          csum = 0.0
          for i in range(len(self.eps)):
              for j in range(len(self.eps)):
                  csum += (self.Gamma1RR[i,j]*self.Gamma2LL[j,i]
                    *(special.digamma(0.5+(BETA/(2.0j*np.pi))*(self.eps[i].conj()-V1))
                     -special.digamma(0.5-(BETA/(2.0j*np.pi))*(self.eps[j]-V1))
                     -special.digamma(0.5+(BETA/(2.0j*np.pi))*(self.eps[i].conj()-V2))
                     +special.digamma(0.5-(BETA/(2.0j*np.pi))*(self.eps[j]-V2))
                     )/(self.eps[i].conj()-self.eps[j])
                    )
          return np.real(csum)

      def integralSparse(self,V1, V2,eps,BETA):
          # Nested loop replaced by a zipped loop over sparse matrix
          csum = 0.0
          for i,j,v in zip(self.GammaProd.row, self.GammaProd.col, self.GammaProd.data):
              if(abs(np.imag(self.eps[i]))>1.0e-15 and abs(np.imag(self.eps[j]))>1.0e-15):
                  csum += (v*
                           (special.digamma(0.5+(BETA/(2.0j*np.pi))*(self.eps[i].conj()-V1))
                           -special.digamma(0.5-(BETA/(2.0j*np.pi))*(self.eps[j]-V1))
                           -special.digamma(0.5+(BETA/(2.0j*np.pi))*(self.eps[i].conj()-V2))
                           +special.digamma(0.5-(BETA/(2.0j*np.pi))*(self.eps[j]-V2))
                           )/(self.eps[i].conj()-self.eps[j])
                          )
          return np.real(csum)

      def integralContract(self,V1, V2,eps,BETA):
          # Explicit loops replaced by numpy broadcasting and array operations
          C = ((special.digamma(0.5+(BETA/(2.0j*np.pi))*(self.eps.conj()-V1))
               -special.digamma(0.5+(BETA/(2.0j*np.pi))*(self.eps.conj()-V2)))[:,None]
              +(special.digamma(0.5-(BETA/(2.0j*np.pi))*(self.eps-V2))
               -special.digamma(0.5-(BETA/(2.0j*np.pi))*(self.eps-V1)))[None,:]
               )/(self.eps.conj()[:,None]-self.eps[None,:])
          return np.real(np.trace((self.Gamma1RR * C) @ self.Gamma2LL))
