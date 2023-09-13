from dataclasses import dataclass
from Shiba.basic_functions import logger
import numpy as np
import matplotlib.pyplot as plt
import os

@dataclass
class PlotManager: 
      def manage_plots(self,Parameters,TransportHamiltonian) -> None:  
          if(Parameters.plot):               
             self.plot_sparsity_maps(Parameters,TransportHamiltonian)
             self.plot_conductance_and_spectral_function(Parameters,TransportHamiltonian)

      def plot_sparsity_maps(self,Parameters,TransportHamiltonian) -> None:
          logger('\n')
          logger('Plotting sparsity maps ...\n')
          fig, axs = plt.subplots(1, 2)
          axs[0].spy(np.real(TransportHamiltonian.h), precision=1.0e-10)
          axs[0].set_title(r'Re$[h]$', loc='center')
          axs[1].spy(np.imag(TransportHamiltonian.h), precision=1.0e-10)
          axs[1].set_title(r'Im$[h]$', loc='center')
          fig.savefig(os.path.join(Parameters.outdir,"sparsity_h.pdf"), bbox_inches='tight')
          fig, axs = plt.subplots(1, 2)
          axs[0].spy(np.real(TransportHamiltonian.Gamma1), precision=1.0e-10)
          axs[0].set_title(r'$\Gamma_{\mathrm{STM}}$', loc='center')
          axs[1].spy(np.real(TransportHamiltonian.Gamma2), precision=1.0e-10)
          axs[1].set_title(r'$\Gamma_{\mathrm{SUB}}$', loc='center')
          fig.savefig(os.path.join(Parameters.outdir,"sparsity_gamma.pdf"), bbox_inches='tight')
      
      def plot_conductance_and_spectral_function(self,Parameters,TransportHamiltonian) -> None:
          logger('\n')
          logger('Plotting conductance ...\n')
          fig, axs = plt.subplots(Parameters.spec+2, 1, figsize=(Parameters.spec+2,2*(Parameters.spec+2)))
          if(Parameters.spec):
             axs[0].plot(TransportHamiltonian.w/1.0e-3, TransportHamiltonian.spectral) # energy in millielectronvolts
             axs[0].set_yscale('log')
             axs[0].set_title(r'Spectral function', loc='left')
             axs[0].set_xlabel(r'$E-E_{\mathrm{F}}$ [meV]')
             axs[0].set_ylabel(r'$A(E)$')

          axs[Parameters.spec].plot(TransportHamiltonian.voltage/1.0e-3, TransportHamiltonian.current*Parameters.EVTOA/1.0e-9) # voltage in millivolts, current in nanoamperes
          axs[Parameters.spec].set_title(r'Current', loc='left')
          axs[Parameters.spec].set_xlabel(r'$V$ [mV]')
          axs[Parameters.spec].set_ylabel(r'$I$ [nA]')

          axs[Parameters.spec+1].plot(TransportHamiltonian.voltage/1.0e-3, TransportHamiltonian.conductance) # voltage in milivolts, conductance in 2e^2/h
          if(Parameters.spec):
             axs[Parameters.spec+1].plot(TransportHamiltonian.w/1.0e-3, TransportHamiltonian.transmission, '--')
          axs[Parameters.spec+1].set_title(r'Conductance', loc='left')
          axs[Parameters.spec+1].set_xlabel(r'$V$ [mV]')
          axs[Parameters.spec+1].set_ylabel(r'$dI/dV$ [$2e^2/h$]')

          fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.7)
          fig.savefig(os.path.join(Parameters.outdir,"plot.pdf"), bbox_inches='tight')