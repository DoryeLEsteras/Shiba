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

30 Aug 2023
riku.m.s.tuovinen@jyu.fi

'''

import sys
import numpy as np
from scipy import linalg, special, sparse, constants
import matplotlib.pyplot as plt
import time
start = time.time()

logfile = 'out.log'

def logger(s):
    print(s, end='')
    with open(logfile, 'a') as f:
        f.write(s)

if(len(sys.argv) == 2):
    wfile = str(sys.argv[1])    # Coupled-spin Wannier data
    calc = 1
elif(len(sys.argv) == 3):
    wfiledn = str(sys.argv[1])  # Spin-down Wannier data
    wfileup = str(sys.argv[2])  # Spin-up Wannier data
    calc = 2
else:
    logger('\nIncompatible number of input files, exiting.' + '\n')
    sys.exit(0)

# System parameters (units are in eV, when applicable)
mu      = -1.1141           # Fermi level
Delta   = 0.0               # Superconducting pair potential (phenomenological parameter)
alpha   = 0.0               # Spin-orbit interaction (phenomenological parameter when calc==2)
gamma   = 0.1e-3            # Coupling strength to the electrodes (estimate from peak width)
frac    = 1.0               # Relative coupling: SUB vs. STM (1.0 means equal coupling)
vran    = 5.0e-1            # Voltage range
vpts    = 1000              # Number of voltage points
NC      = 0                 # Depth of Wigner-Seitz overlaps in the supermatrix
spec    = 0                 # Calculate spectral and transmission functions?
opt     = 2                 # Optimization: 0=Eigenvalue loop, 1=Sparse-zip loop, 2=Broadcasting
plot    = 0                 # Make plots?
NSTM    = np.genfromtxt('stm.in', dtype=None)   # Orbitals coupled to the STM
NSUB    = np.genfromtxt('sub.in', dtype=None)   # Orbitals coupled to the substrate

# Set either temperature or beta (inverse temperature)
temperature = 350.0e-3
beta = 1.0/((constants.k/constants.eV)*temperature)
#beta = 33000.0  # Inverse temperature (33000 eV^(-1) ~> 350 mK)
eVtoA = constants.eV**2/constants.hbar  # Conversion factor from electronvolts to amperes

def cell(x,y,z):
    return (abs(x)+abs(y)+abs(z))

##################
# Initialization #
##################
open(logfile, 'w').close()
# Read Wannier data
if(calc==1):
    degs = []
    with open(wfile) as f:
        for row, line in enumerate(f):
            if(row==0):
                continue
            elif(row==1):
                norb = [int(number) for number in line.split()][0]
            elif(row==2):
                nws = [int(number) for number in line.split()][0]
            elif(row<=nws/15+3):
                degs.extend([float(number) for number in line.split()])
            else:
                break
    mels = np.genfromtxt(wfile, skip_header=nws//15+4)
    logger('\nWannier orbitals:\t' + str(norb) + '\n')
    logger('Wigner-Seitz cells:\t' + str(nws) + '\n')
    invdegs = np.zeros_like(degs)
    for i in range(len(degs)):
        invdegs[i] = 1.0/degs[i]
    logger('(Reciprocal grid:\t' + str(int(np.sum(invdegs))) + ')' + '\n')
    logger('Depth of WS overlaps:\t' + str(NC) + '\n')
    blocks=0
    for idx in range(nws):
        if(cell(mels[idx*norb**2,0],mels[idx*norb**2,1],mels[idx*norb**2,2])<=NC):
            blocks+=1
    logger('Full Hamiltonian size:\t' + str((blocks*4*(norb//2))) + 'x' + str((blocks*4*(norb//2))) + '\n\n')
else:
    degs = []
    with open(wfiledn) as f:
        for row, line in enumerate(f):
            if(row==0):
                continue
            elif(row==1):
                norbdn = [int(number) for number in line.split()][0]
            elif(row==2):
                nws = [int(number) for number in line.split()][0]
            elif(row<=nws/15+3):
                degs.extend([float(number) for number in line.split()])
            else:
                break
    melsdn = np.genfromtxt(wfiledn, skip_header=nws//15+4)
    invdegs = np.zeros_like(degs)
    for i in range(len(degs)):
        invdegs[i] = 1.0/degs[i]
    with open(wfileup) as f:
        for row, line in enumerate(f):
            if(row==0):
                continue
            elif(row==1):
                norbup = [int(number) for number in line.split()][0]
            else:
                break
    melsup = np.genfromtxt(wfileup, skip_header=nws//15+4)
    logger('\nWannier orbitals:\t' + str(norbdn)+"+"+str(norbup)+"="+str(norbdn+norbup) + '\n')
    logger('Wigner-Seitz cells:\t' + str(nws) + '\n')
    logger('(Reciprocal grid:\t' + str(int(np.sum(invdegs))) + ')' + '\n')
    logger('Depth of WS overlaps:\t' + str(NC) + '\n')
    norb=norbdn+norbup
    blocks=0
    for idx in range(nws):
        if(cell(melsdn[idx*norbdn**2,0],melsdn[idx*norbdn**2,1],melsdn[idx*norbdn**2,2])<=NC):
            blocks+=1
    logger('Full Hamiltonian size:\t' + str((blocks*4*(norb//2))) + 'x' + str((blocks*4*(norb//2))) + '\n\n')

###########################
# System characterization #
###########################

h       = np.zeros((blocks*4*(norb//2),blocks*4*(norb//2)), dtype=complex)
Id      = np.eye(blocks*4*(norb//2), dtype=complex)
Idnorb2 = np.eye(norb//2, dtype=complex)
if(calc==2):
    Idnorbdn = np.eye(norbdn, dtype=complex)
    Idnorbup = np.eye(norbup, dtype=complex)

# Loop over WS cells
ibl=-1
for idx in range(nws):
    # Take only the WS overlaps with depth up to NC
    if(calc==1):
        if(cell(mels[idx*norb**2,0],mels[idx*norb**2,1],mels[idx*norb**2,2])<=NC):
            ibl+=1
            # For each WS overlap, extract different spin components (see WF_Order.txt)
            huu = np.zeros((norb//2,norb//2), dtype=complex)
            hdd = np.zeros((norb//2,norb//2), dtype=complex)
            hud = np.zeros((norb//2,norb//2), dtype=complex)
            hdu = np.zeros((norb//2,norb//2), dtype=complex)

            ie=-1
            io=-1
            for i in range(norb):
                if(i%2==0):
                    ie+=1
                else:
                    io+=1
                je=-1
                jo=-1
                for j in range(norb):
                    if(j%2==0):
                        je+=1
                    else:
                        jo+=1

                    if(i%2==0 and j%2==0):
                        huu[ie,je]=mels[idx*norb**2+i*norb+j,5]+1.0j*mels[idx*norb**2+i*norb+j,6]
                    if(i%2==0 and j%2!=0):
                        hud[ie,jo]=mels[idx*norb**2+i*norb+j,5]+1.0j*mels[idx*norb**2+i*norb+j,6]
                    if(i%2!=0 and j%2==0):
                        hdu[io,je]=mels[idx*norb**2+i*norb+j,5]+1.0j*mels[idx*norb**2+i*norb+j,6]
                    if(i%2!=0 and j%2!=0):
                        hdd[io,jo]=mels[idx*norb**2+i*norb+j,5]+1.0j*mels[idx*norb**2+i*norb+j,6]

            # Align Fermi level at zero
            huu -= mu*Idnorb2
            hdd -= mu*Idnorb2

            # Construct the corresponding blocks of the full Nambu Hamiltonian
            row = ibl*4*(norb//2) # Block diagonal super matrix
            col = ibl*4*(norb//2)
            for i in range(0,4*(norb//2),4):
                for j in range(0,4*(norb//2),4):
                    h[row+i+0,col+j+0] = huu[i//4,j//4]
                    h[row+i+0,col+j+1] = -Delta*(i==j)
                    h[row+i+0,col+j+2] = hud[i//4,j//4]
                    h[row+i+0,col+j+3] = 0.0
                    h[row+i+1,col+j+0] = -Delta*(i==j)
                    h[row+i+1,col+j+1] = -hdd[i//4,j//4]
                    h[row+i+1,col+j+2] = 0.0
                    h[row+i+1,col+j+3] = -hdu[i//4,j//4]
                    h[row+i+2,col+j+0] = hdu[i//4,j//4]
                    h[row+i+2,col+j+1] = 0.0
                    h[row+i+2,col+j+2] = hdd[i//4,j//4]
                    h[row+i+2,col+j+3] = Delta*(i==j)
                    h[row+i+3,col+j+0] = 0.0
                    h[row+i+3,col+j+1] = -hud[i//4,j//4]
                    h[row+i+3,col+j+2] = Delta*(i==j)
                    h[row+i+3,col+j+3] = -huu[i//4,j//4]
    else:
        if(cell(melsdn[idx*norbdn**2,0],melsdn[idx*norbdn**2,1],melsdn[idx*norbdn**2,2])<=NC):
            ibl+=1
            # For each WS overlap, extract different spin components (see WF_Order.txt)
            if(norbup>=norbdn):
                huu = np.zeros((norbup,norbup), dtype=complex)
                hdd = np.zeros((norbup,norbup), dtype=complex)
            else:
                huu = np.zeros((norbdn,norbdn), dtype=complex)
                hdd = np.zeros((norbdn,norbdn), dtype=complex)

            for i in range(norbdn):
                for j in range(norbdn):
                    hdd[i,j]=melsdn[idx*norbdn**2+i*norbdn+j,5]+1.0j*melsdn[idx*norbdn**2+i*norbdn+j,6]

            for i in range(norbup):
                for j in range(norbup):
                    huu[i,j]=melsup[idx*norbup**2+i*norbup+j,5]+1.0j*melsup[idx*norbup**2+i*norbup+j,6]

            # Align Fermi level at zero
            hdd -= mu*Idnorbdn
            huu -= mu*Idnorbup

            # Construct the corresponding blocks of the full Nambu Hamiltonian
            row = ibl*4*(norb//2) # Block diagonal super matrix
            col = ibl*4*(norb//2)
            for i in range(0,4*(norb//2),4):
                for j in range(0,4*(norb//2),4):
                    h[row+i+0,col+j+0] = huu[i//4,j//4]
                    h[row+i+0,col+j+1] = -Delta*(i==j)
                    h[row+i+0,col+j+2] = alpha
                    h[row+i+0,col+j+3] = 0.0
                    h[row+i+1,col+j+0] = -Delta*(i==j)
                    h[row+i+1,col+j+1] = -hdd[i//4,j//4]
                    h[row+i+1,col+j+2] = 0.0
                    h[row+i+1,col+j+3] = -alpha
                    h[row+i+2,col+j+0] = alpha
                    h[row+i+2,col+j+1] = 0.0
                    h[row+i+2,col+j+2] = hdd[i//4,j//4]
                    h[row+i+2,col+j+3] = Delta*(i==j)
                    h[row+i+3,col+j+0] = 0.0
                    h[row+i+3,col+j+1] = -alpha
                    h[row+i+3,col+j+2] = Delta*(i==j)
                    h[row+i+3,col+j+3] = -huu[i//4,j//4]

# The Wannier matrix-element data may not always be symmetric, so
# the Hamiltonian is symmetrized to assure proper analytic structure
h = 0.5*(h+h.conj().T)
logger("The full Hamiltonian is hermitian? " + str(np.all(np.abs(h-h.conj().T) < 1.0e-12)) + '\n')

# Construct coupling matrices
Id4     = np.eye(4, dtype=complex)
gam1    = np.zeros((4*(norb//2),4*(norb//2)), dtype=complex)
gam2    = np.zeros((4*(norb//2),4*(norb//2)), dtype=complex)

# Each orbital from NSTM is coupled to the STM orbital (full)
for i in NSTM:
    for j in NSTM:
        gam1[4*(i-1):4*i,4*(j-1):4*j] = gamma*Id4

# Each orbital from NSUB is coupled to individual substrate orbitals (diagonal)
for i in NSUB:
    gam2[4*(i-1):4*i,4*(i-1):4*i] = frac*gamma*Id4

# Each WS overlap block has the same coupling matrix structure,
# so we duplicate the individual blocks as kronecker products
blockstructure=np.eye(blocks, dtype=complex) # Block diagonal super matrix
Gamma1 = np.kron(blockstructure,gam1)
Gamma2 = np.kron(blockstructure,gam2)

if(plot):
    # Sparsity maps of the Hamiltonian matrix
    fig, axs = plt.subplots(1, 2)
    axs[0].spy(np.real(h), precision=1.0e-10)
    axs[0].set_title(r'Re$[h]$', loc='center')
    axs[1].spy(np.imag(h), precision=1.0e-10)
    axs[1].set_title(r'Im$[h]$', loc='center')
    fig.savefig("sparsity_h.pdf", bbox_inches='tight')
    fig, axs = plt.subplots(1, 2)
    axs[0].spy(np.real(Gamma1), precision=1.0e-10)
    axs[0].set_title(r'$\Gamma_{\mathrm{STM}}$', loc='center')
    axs[1].spy(np.real(Gamma2), precision=1.0e-10)
    axs[1].set_title(r'$\Gamma_{\mathrm{SUB}}$', loc='center')
    fig.savefig("sparsity_gamma.pdf", bbox_inches='tight')

# Effective Hamiltonian (coupled system)
heff = h - 0.5j*(Gamma1+Gamma2)

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

####################
# Timer and finish #
####################

def getTime(value):
    timestr = ''
    if(value == 0):
        timestr += '00'
    elif(value < 10):
        timestr += '0%s'%value
    else:
        timestr += str(value)
    return timestr

def StopWatch(start, finish):
    s_in_h = 3600
    s_in_min = 60
    elapsed = int(finish-start)
    mystr = ''
    mystr += 'Time wasted: '
    if(elapsed > 0):
        mystr += getTime(int(elapsed / s_in_h))
        mystr += ':'
        mystr += getTime(int((elapsed % s_in_h) / s_in_min))
        mystr += ':'
        mystr += getTime(int((elapsed % s_in_h) % (s_in_min)))
    else:
        mystr += '00:00:00'
    mystr += '\n'
    return mystr

logger('\n')
logger('Done! ' + StopWatch(start, time.time()) + '\n')
