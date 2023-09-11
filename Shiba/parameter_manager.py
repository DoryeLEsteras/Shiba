from typing import List
import numpy as np
from argparse import ArgumentParser
from dataclasses import dataclass
from scipy import linalg, special, sparse, constants

def parser():
   parser = ArgumentParser(description="Script to comput Yu-Shiba-Rusinov states")
   parser.add_argument("-input", "--input",
                        type=str,
                        required=True,
                        help="""Input file to execute Shiba""")
   
   parser.add_argument("-hfile", "--hfile",
                        type=str,
                        required=True,
                        nargs='+',
                        help="""File (noncolinear) or files (colinear) including
                         the up and down components of the Hamiltonian""")

   parser.add_argument("-stm", "--stm",
                        type=str,
                        required=True,
                        help="File indicating the orbitals to be coupled to the stm tip")
  
   parser.add_argument("-sub", "--sub",
                        type=str,
                        required=True,
                        help="File indicating the orbitals to be coupled to the substrate")

   parser.add_argument("-outdir", "--outdir",
                        type=str,
                        required=False,
                        default='./',
                        help="Path to store the outputs")
   
   args = parser.parse_args()
   return args.input, args.hfile, args.stm, args.sub, args.outdir

@dataclass
class ParameterManager:
      # Constants
      mu: float = -0.1452  # CONSTANTS should be in capital letters
      temperature: float = 350.0e-3
      beta: np.ndarray = 1.0/((constants.k/constants.eV)*temperature)
      eVtoA: float = constants.eV**2/constants.hbar

      #Input parameters (Default values)
      delta: float = 0.0
      alpha: float = 0.0
      gamma: float = 0.1e-3
      frac: float = 1.0
      vran: float = 5.0e-1
      vpts: int = 1000
      nc: int = 0
      spec: bool = False
      opt: int = 2
      plot: bool = True
      NSTM: np.ndarray = np.array([])
      NSUB: np.ndarray = np.array([])
      # Files
      shiba_input: str = ''
      hfile: str = ''
      stm: str = ''
      sub: str = ''
      outdir: str = ''


      def extract_input_information(self) -> None:
          self.shiba_input, self.hfile, self.stm, self.sub, self.outdir = parser()
          self.NSTM = np.genfromtxt(self.stm,dtype=None) # var Should be lower case
          self.NSUB = np.genfromtxt(self.sub,dtype=None)
          input_file = open(self.shiba_input, 'r')
          for line_number, line in enumerate(input_file): 
            splitted_line = line.replace('\n','').replace(' ','').split('='); 
            if len(splitted_line) > 1:
               if splitted_line[0].lower() == 'delta':
                  self.delta = float(splitted_line[1])
               if splitted_line[0].lower() == 'alpha':
                  self.alpha = float(splitted_line[1])
               if splitted_line[0].lower() == 'gamma':
                  self.gamma = float(splitted_line[1])
               if splitted_line[0].lower() == 'frac':
                  self.frac = float(splitted_line[1])
               if splitted_line[0].lower() == 'vran':
                  self.vran = float(splitted_line[1])
               if splitted_line[0].lower() == 'vpts':
                  self.vpts = int(splitted_line[1])
               if splitted_line[0].lower() == 'nc': 
                  self.nc = int(splitted_line[1])
               if splitted_line[0].lower() == 'spec': 
                  self.spec = bool(splitted_line[1])
               if splitted_line[0].lower() == 'opt': 
                  self.opt = int(splitted_line[1])
               if splitted_line[0].lower() == 'plot':  
                  if splitted_line[1] in '0yY' or splitted_line[1].lower()== 'true': 
                     self.plot = True
                  elif splitted_line[1] in '1nN' or splitted_line[1].lower()== 'false':
                     self.plot = False

if __name__ == '__main__':
   print('hi')