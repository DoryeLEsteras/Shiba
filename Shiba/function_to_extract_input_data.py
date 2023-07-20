from typing import List
import numpy as np

# We need to discuss what Riku needs from the old WF order file. Fermi and hr file will be readed via parser
def extract_input_information(file_name:str) -> List:
    with open(file_name) as input_file:
        input_vector = input_file.readlines()
    for line_number,line in enumerate(input_vector):
        readed_line = line.split()
        if len(readed_line) >= 2:
           if readed_line[0] == 'Grid' and readed_line[1] == 'size': # kp
              nk = int(readed_line[-1])
              kpoints_rel = np.zeros((nk,3))
              kpoints_cart = np.zeros((nk,3))
              for index in range(line_number +5,line_number +5 + nk):
                  kpoints_rel[index - line_number -5] = input_vector[index].split()[2:5]
                  kpoints_cart[index - line_number -5] = input_vector[index].split()[6:9]
        if line.strip().endswith('Wigner-Seitz supercell:'):     # degeneracies
            nwigner = int(input_vector[line_number].split()[0])
            degeneracies = np.zeros((nwigner,4))
            for index in range(line_number + 2,line_number + nwigner +2):
                degeneracies[index - line_number-2] =   \
                            input_vector[index].split()[1:4] + \
                            input_vector[index].split()[5:6]
    return kpoints_rel,kpoints_cart,degeneracies

if __name__ == '__main__':
    kpoints_rel,kpoints_cart,degeneracies = \
                    extract_input_information('1Fe_centreconf.up.wout')
    print('kpoints (relative)')
    print(kpoints_rel)
    print('\n kpoints (cartesian)')
    print(kpoints_cart)
    print('\n degeneracies')
    print(degeneracies)
