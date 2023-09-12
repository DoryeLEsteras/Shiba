'''

This code calculates current-voltage characteristics of a molecular structure 
coupled to wide bandwidth electrodes. The molecular structure is described by 
Wannier orbitals (input data from Wannier90) and Nambu spinor representation. 
The current is obtained from the Landauer-Büttiker formula via an efficient 
calculation of nonequilibrium Green's functions. The code can read spinful 
(one cmd line argument for input file) and spin-resolved (two cmd line arguments 
for the spin-down and spin-up input files) Wannier input data. The code can be 
used to study, for instance, magnetic molecular effects in the vicinity of a 
superconductor, such as Yu-Shiba-Rusinov states as in-gap conductance peaks.

Example runs

CrBr3-NbSe2 (single Wannier input):
python3 wstm.py crbr3.nbse2.bilayer_hr.dat

Fe-NbSe2 (separate Wannier input for spin-down and spin-up):
python3 wstm.py nbse2.down_hr.dat nbse2.up_hr.dat

riku.m.s.tuovinen@jyu.fi
dorye.esteras@uv.es

'''

from dataclasses import dataclass
from Shiba.wannier_reader import WannierHamiltonian
from Shiba.basic_functions import logger,get_time,stop_watch,estimate_memory
from Shiba.parameter_manager import parser, ParameterManager
from Shiba.transport_calculations import TransportCalculation
from Shiba.plotter import PlotManager

import sys
import numpy as np
from scipy import linalg, special, sparse, constants
import matplotlib.pyplot as plt
import time





if __name__ == '__main__':
   start = time.time()

   logfile = 'out.log'  
   open(logfile, 'w').close()

   Parameters = ParameterManager()
   Parameters.extract_input_information()
   WannierParameters = WannierHamiltonian()
   WannierParameters.read_wannier_info(Parameters.input_dir,Parameters.up_h, 
                     Parameters.down_h,Parameters.noncolin_h,Parameters.nc)

   TransportHamiltonian = TransportCalculation()
   TransportHamiltonian.construct_nambu(Parameters,WannierParameters)
   TransportHamiltonian.construct_coupling_matrices(Parameters,WannierParameters)
   TransportHamiltonian.solve_hamiltonian(Parameters)
   Plotter = PlotManager()
   Plotter.manage_plots(Parameters,TransportHamiltonian)
   estimate_memory(Parameters,WannierParameters)
    
####################
# Timer and finish #
####################

logger('\n')
logger('Done! ' + stop_watch(start, time.time()) + '\n')





"""
# Solve non-hermitian eigensystem (separate left/right eigenvectors)
logger('\n')
logger('Diagonalization ...\n')
eps, psil, psir = linalg.eig(heff, left=True, right=True)
overlap_lr = psil.conj().T @ psir
invover_lr = np.linalg.inv(overlap_lr)
overlap_rl = psir.conj().T @ psil
invover_rl = np.linalg.inv(overlap_rl)
Gamma1RR = psir.conj().T @ Gamma1 @ psir @ invover_lr
Gamma2LL = psil.conj().T @ Gamma2 @ psil @ invover_rl
if(opt==1):
    GammaProd = sparse.coo_matrix(Gamma1RR*Gamma2LL.T)

del psil, psir, overlap_lr, invover_lr, overlap_rl, invover_rl

###############
# Calculation #
###############

# Analytically resolved frequency integral
def integralOriginal(V1, V2):
    # Original result as a nested loop over eigenvalues
    csum = 0.0
    for i in range(len(eps)):
        for j in range(len(eps)):
            csum += (Gamma1RR[i,j]*Gamma2LL[j,i]
                    *(special.digamma(0.5+(beta/(2.0j*np.pi))*(eps[i].conj()-V1))
                     -special.digamma(0.5-(beta/(2.0j*np.pi))*(eps[j]-V1))
                     -special.digamma(0.5+(beta/(2.0j*np.pi))*(eps[i].conj()-V2))
                     +special.digamma(0.5-(beta/(2.0j*np.pi))*(eps[j]-V2))
                     )/(eps[i].conj()-eps[j])
                    )
    return np.real(csum)

def integralSparse(V1, V2):
    # Nested loop replaced by a zipped loop over sparse matrix
    csum = 0.0
    for i,j,v in zip(GammaProd.row, GammaProd.col, GammaProd.data):
        if(abs(np.imag(eps[i]))>1.0e-15 and abs(np.imag(eps[j]))>1.0e-15):
            csum += (v*
                     (special.digamma(0.5+(beta/(2.0j*np.pi))*(eps[i].conj()-V1))
                     -special.digamma(0.5-(beta/(2.0j*np.pi))*(eps[j]-V1))
                     -special.digamma(0.5+(beta/(2.0j*np.pi))*(eps[i].conj()-V2))
                     +special.digamma(0.5-(beta/(2.0j*np.pi))*(eps[j]-V2))
                     )/(eps[i].conj()-eps[j])
                    )
    return np.real(csum)

def integralContract(V1, V2):
    # Explicit loops replaced by numpy broadcasting and array operations
    C = ((special.digamma(0.5+(beta/(2.0j*np.pi))*(eps.conj()-V1))
         -special.digamma(0.5+(beta/(2.0j*np.pi))*(eps.conj()-V2)))[:,None]
        +(special.digamma(0.5-(beta/(2.0j*np.pi))*(eps-V2))
         -special.digamma(0.5-(beta/(2.0j*np.pi))*(eps-V1)))[None,:]
         )/(eps.conj()[:,None]-eps[None,:])
    return np.real(np.trace((Gamma1RR * C) @ Gamma2LL))

if(spec):
    logger('\n')
    logger('Spectral function ...\n')
    w = np.linspace(-vran, vran, vpts)
    spectral = np.zeros_like(w)
    transmission = np.zeros_like(w)
    for i in range(len(w)):
        # Spectral function is calculated as the imaginary part of the retarded Green function
        # Isolated system with artificial broadening (i*eta)
    #    spectral[i] = (-1.0/np.pi)*np.imag(np.trace(np.linalg.inv((w[i]+1.0e-12j)*Id - h)))
        # Coupled system with physical broadening (i*Gamma)
        A = np.linalg.inv(w[i]*Id - heff)
        spectral[i] = (-1.0/np.pi)*np.imag(np.trace(A))
        # Transmission function
        transmission[i] = np.real(np.trace(Gamma1 @ A @ Gamma2 @ A.conj().T))

    np.savetxt('specdata.out', np.array([w,spectral,transmission]).T, fmt=['%.8e', '%.8e', '%.8e'])

del h, heff, Gamma1, Gamma2

logger('\n')
logger('Current-voltage ...\n')
voltage = np.linspace(-vran, vran, vpts)
current = np.zeros_like(voltage)
conductance = np.zeros_like(voltage)
for i in range(len(voltage)):
    if(opt==0):
        current[i] = integralOriginal(voltage[i],0.0)
    elif(opt==1):
        current[i] = integralSparse(voltage[i],0.0)
    elif(opt==2):
        current[i] = integralContract(voltage[i],0.0)
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
conductance = np.gradient(current, voltage[1]-voltage[0])

np.savetxt('ivdata.out', np.array([voltage,current,conductance]).T, fmt=['%.8e', '%.8e', '%.8e'])

#################
# Visualization #
#################

if(plot):
    logger('\n')
    logger('Plotting ...\n')

    fig, axs = plt.subplots(spec+2, 1, figsize=(spec+2,2*(spec+2)))

    if(spec):
        axs[0].plot(w/1.0e-3, spectral) # energy in millielectronvolts
        axs[0].set_yscale('log')
        axs[0].set_title(r'Spectral function', loc='left')
        axs[0].set_xlabel(r'$E-E_{\mathrm{F}}$ [meV]')
        axs[0].set_ylabel(r'$A(E)$')

    axs[spec].plot(voltage/1.0e-3, current*eVtoA/1.0e-9) # voltage in millivolts, current in nanoamperes
    axs[spec].set_title(r'Current', loc='left')
    axs[spec].set_xlabel(r'$V$ [mV]')
    axs[spec].set_ylabel(r'$I$ [nA]')

    axs[spec+1].plot(voltage/1.0e-3, conductance) # voltage in millivolts, conductance in 2e^2/h
    if(spec):
        axs[spec+1].plot(w/1.0e-3, transmission, '--')
    axs[spec+1].set_title(r'Conductance', loc='left')
    axs[spec+1].set_xlabel(r'$V$ [mV]')
    axs[spec+1].set_ylabel(r'$dI/dV$ [$2e^2/h$]')

    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.7)
    fig.savefig("plot.pdf", bbox_inches='tight')
"""

