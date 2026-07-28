#This code is an attempt at a faithful reproduction of the make_symmdef_file.pl created by Frank DiMaio and others. I, Miles Bradford,
#am creating this code with the help of ChatGPT in order to help me understand how the original code works. Upon completion, I will create my own
#script that reproduces the functionality of this script only insofar as it is necessary to be compatible with Rosetta's symmetry tools, while
#functioning quite differently under the hood than the original script, since it is very complicated. My script will do a lot less than the original,
#and will expect the user to know a little more, but it will be well-annotated and will be perfectly clear about what each step is doing.

import numpy as np

def generate_symmdef(pdb_path,
                     output_symm_path,
                     output_pdb_path,
                     interact_dist):
    """This is the main function of make_symmdef_file.py. It takes a PDB and an interaction distance and generates the symmetry file and symmetric PDB
    that correspond to the passed PDB's contents and CRYST1 line."""
    #1. Read PDB
    #2. Read crystal parameters
    #3. Build backbone trace
    #4. Compute centroid
    #5. Load spacegroup operators
    #6. Find symmetry interfaces
    #7. Pair interfaces
    #8. Write .symm
    #9. Write expanded PDB
    backbone_atom_coords:list[np.ndarray] = []
    cartesian_com = np.zeros(3)
    fractional_com = np.zeros(3)
    
    #backbone_atom_coords is just a list of all the backbone atom coordinates.
