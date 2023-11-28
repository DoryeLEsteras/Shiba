
***********
How to use?
***********

Precalculations
===============

To compute Yu-Shiba-Rusinov states using CHECK_PACKAGE_NAME, previous calculations
based in DFT + Wannier90 are necessary. A DFT package can be used to simulate the system of interest, and to compute the groundstate. At this level of theory, superconductivity is not considered, however magnetism and the Hubbard model should be described correctly.

DFT calculations are followed by the computation of Wannier functions. The objectivehere consist in the obtention of a Wannier Hamiltonian that successfully reproduces the bands around the Fermi level. Both, the fit to the DFT bands and the spreads of Wannier functions should be checked to ensure a correct behaviour of the program.

Starting using CHECK_PACKAGE_NAME
=================================
At this point you are ready to start using CHECK_PACKAGE_NAME!
First, let's start creating an input template with the following command:

.. code-block:: python

   create_template.py -input name_for_you_input -outdir .

.. note::
   The flag 'input' is mandatory, it gives a name to your input. The flag 'outdir' is optional, by default will use your actual directory to store the template.

After this, you will have an input template like this: CHECK_TEMPLATE_UPDATE

.. code-block:: fortran

   up_sub_file = sub_up.in        ! CHECK EXPLAIN EVERYTHIN ANDUNITS
   down_sub_file = sub_down.in    !
   up_stm_file = stm_up.in        !
   down_stm_file = stm_down.in    !
   # comb_sub_file = sub_comb.in  !
   # comb_stm_file = stm_comb.in  !
                                  !
   up_h = nbse2.up_hr.dat         !
   down_h = nbse2.down_hr.dat     !
   # comb_h = nbse2.up_hr.dat     !
                                  !
   delta = 0.1                    !
   mu = 0.0                       !
   alpha = 0.01                   !
   gamma = 0.0001                 !
   frac = 1                       !
   vpts = 500                     !
   vran = 0.5                     !
   nc= 0                          !
   spec = 0                       !
   opt = 2                        !
   plot = 0                       ! 
   temperature = 350e-3           !

.. note::
 The input is not sensitive to capital letters, commas, spaces or empty lines.
 A line can be commented writting any character in front of the keywords.
 Examples: 
 #delta
 ! delta
 This is a comment delta

Running calculations
====================

Nice! You configured your input! Please take care of the next considerations:

CHECK_CALCULATE_DELTA
CHECK_CONVERGENCE_OF_PARAMETERS

With this you can run the code with the following command:

.. code-block:: bash 

   CHECK_COMMAND_TO_RUN_CODE
