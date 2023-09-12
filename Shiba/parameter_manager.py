import numpy as np
import os
from argparse import ArgumentParser
from dataclasses import dataclass
from scipy import linalg, special, sparse, constants

def parser():
   parser = ArgumentParser(description="Script to comput Yu-Shiba-Rusinov states")
   parser.add_argument("-input", "--input",
                        type=str,
                        required=True,
                        help="""Input file to execute Shiba""")
   
   parser.add_argument("-outdir", "--outdir",
                        type=str,
                        required=False,
                        default='',
                        help="Path to store the outputs")
   
   args = parser.parse_args()
   return args.input, args.outdir

@dataclass
class ParameterManager: 
      #Input parameters (Default values)
      temperature: float = 350.0e-3
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
      nstm: np.ndarray = np.array([])
      nsub: np.ndarray = np.array([])

      # Files
      input_dir: str = ''
      outdir: str = ''
      shiba_input: str = ''
      up_h: str = ''
      down_h: str = ''
      noncolin_h: str = ''
      stm: str = ''
      sub: str = ''

      # Constants
      MU: float = -0.1452  # CONSTANTS should be in capital letters
      EVTOA: float = constants.eV**2/constants.hbar
      BETA: np.ndarray = 1.0/((constants.k/constants.eV)*temperature)
      
      def extract_input_information(self) -> None:
          shiba_input_and_dir, self.outdir = parser()
          self.input_dir = os.path.abspath(os.path.dirname(shiba_input_and_dir))
          self.shiba_input = os.path.basename(shiba_input_and_dir)
          if self.outdir == '':
             self.outdir = self.input_dir
          input_file = open(os.path.join(self.input_dir,self.shiba_input), 'r')
          for line_number, line in enumerate(input_file): 
            splitted_line = line.replace('\n','').replace(' ','').split('='); 
            if len(splitted_line) > 1:
               if splitted_line[0].lower() == 'stm_file':
                  self.stm = splitted_line[1]
               if splitted_line[0].lower() == 'sub_file':
                  self.sub = splitted_line[1]
               if splitted_line[0].lower() == 'up_h':
                  self.up_h = splitted_line[1]
               if splitted_line[0].lower() == 'down_h':
                  self.down_h = splitted_line[1]
               if splitted_line[0].lower() == 'noncolin_h':
                  self.noncolin_h = splitted_line[1]
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
          self.NSTM = np.genfromtxt(os.path.join(self.input_dir,self.stm),dtype=None) # var Should be lower case
          self.NSUB = np.genfromtxt(os.path.join(self.input_dir,self.sub),dtype=None)

if __name__ == '__main__':
   print('hi')