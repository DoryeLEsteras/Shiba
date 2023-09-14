from argparse import ArgumentParser
import os

def parser():
   parser = ArgumentParser(description="Script to generate input template")
   parser.add_argument("-name", "--name",
                        type=str,
                        required=True,
                        help="""Input file to execute Shiba""")
   
   parser.add_argument("-outdir", "--outdir",
                        type=str,
                        required=False,
                        default='',
                        help="Path to store the outputs")
   
   args = parser.parse_args()
   return args.name, args.outdir


if __name__ == '__main__':
   input_name, outdir = parser()

   opened_file = open(os.path.join(outdir,input_name),'w')

   opened_file.write('#################################\n')
   opened_file.write('############# INPUT #############\n')
   opened_file.write('#################################\n')
   opened_file.write('sub_file = sub.in\n')
   opened_file.write('stm_file = stm.in\n')
   opened_file.write('up_h = nbse2.up_hr.dat\n')
   opened_file.write('down_h = nbse2.down_hr.dat\n')
   opened_file.write('# noncolin_h = nbse2.up_hr.dat\n')
   opened_file.write('delta = 0.1\n')
   opened_file.write('mu = 0.0\n')
   opened_file.write('alpha = 0.01\n')
   opened_file.write('gamma = 0.0001\n')
   opened_file.write('frac = 1\n')
   opened_file.write('vpts = 500\n')
   opened_file.write('vran = 0.5\n')
   opened_file.write('nc= 0\n')
   opened_file.write('spec = 0\n')
   opened_file.write('opt = 2\n')
   opened_file.write('plot = 0\n')
   opened_file.write('temperature = 350e-3\n\n')
   opened_file.write('#################################\n')
   opened_file.write('############## END ##############\n')
   opened_file.write('#################################\n\n')
   opened_file.write('# The input is not sensitive to:\n')
   opened_file.write('#capital letters, commas and spaces\n\n')
   opened_file.write('# Every symbol before the flags\n')
   opened_file.write('# will work as a comment\n')
   opened_file.write('# Exs:  \n')
   opened_file.write('       # frac\n')
   opened_file.write('	   #frac\n')
   opened_file.write('	   ! frac\n')
   
   opened_file.close()
