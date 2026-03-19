import sys

import pyrosetta
from pyrosetta import *
from pyrosetta import PyMOLMover
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.scoring import get_score_function

pyrosetta.init()
#Now render in PyMOL
pmm = PyMOLMover()
pmm.keep_history(True)
symm_mover = SetupForSymmetryMover(sys.argv[2])
TELSAM = Pose()
pose_from_pdb(TELSAM,sys.argv[1])
symm_mover.apply(TELSAM)
print(f'Symmetric pose: {is_symmetric(TELSAM)}')
sf = get_score_function()
relax = FastRelax()
relax.set_scorefxn(sf)
relax.apply(TELSAM)
pmm.apply(TELSAM)
