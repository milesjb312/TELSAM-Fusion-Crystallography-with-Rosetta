#Relax

import pyrosetta
from pyrosetta import *
from pyrosetta.toolbox import cleanATOM
import pyrosetta.rosetta.core.io
dir(pyrosetta.rosetta.core.io)
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.scoring import get_score_function
#from pyrosetta.rosetta.core.scoring import get_fa_score_function

pyrosetta.init()

cleanATOM("/home/milesjb/CRY2.pdb")
pose = pose_from_pdb("/home/milesjb/CRY2.clean.pdb")
#Refine
sf = get_score_function()
relax = FastRelax()
relax.set_scorefxn(sf)
relax.apply(pose)
dump_pdb(pose,"/home/milesjb/CRY2.relaxed.pdb")