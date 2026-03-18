import sys
import requests
import os

import pyrosetta
from pyrosetta import *
from pyrosetta.toolbox import cleanATOM
from pyrosetta.rosetta.core.pose import append_subpose_to_pose
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.id import AtomID_Map_AtomID as AtomID_Map
from pyrosetta.rosetta.core.pose import initialize_atomid_map
from pyrosetta.rosetta.core.scoring import superimpose_pose

if len(sys.argv)==2:
	#If one argument is passed in, then just grab the pdb and return both its original format and the cleaned version.
	pyrosetta.init()
	url = f"https://files.rcsb.org/download/{sys.argv[1]}.pdb"
	pdb_text = requests.get(url).text
	with open(f"./TELSetta/{sys.argv[1]}.pdb","w") as file:
		file.write(pdb_text)
	cleanATOM(f"./TELSetta/{sys.argv[1]}.pdb")
	temp_pose = pose_from_pdb(f'./TELSetta/{sys.argv[1]}.clean.pdb')
	final_pose = Pose()
	append_subpose_to_pose(final_pose,temp_pose,temp_pose.chain_begin(1),temp_pose.chain_end(1))
	final_pose.dump_pdb(f'./TELSetta/{sys.argv[1]}.clean.pdb')
	#os.remove(f"./TELSetta/{sys.argv[1]}.pdb")
	
if len(sys.argv) == 3:
	#If two arguments are passed in, then align the first protein (engineered TELSAM) to the second protein (TELSAM with the proper space group),
	#then export the pdb of the first protein along with the cryst1 line of the second protein.
	pyrosetta.init()
	#First, import the engineered TELSAM monomer.
	url = f"https://files.rcsb.org/download/{sys.argv[1]}.pdb"
	pdb_text = requests.get(url).text
	with open("ETEL.pdb","w") as file:
		file.write(pdb_text)
	cleanATOM("ETEL.pdb")
	temp_pose = pose_from_pdb('ETEL.clean.pdb')
	E_pose = Pose()
	#Only take the first part of the TELSAM monomer, not whatever client it was attached to.
	append_subpose_to_pose(E_pose,temp_pose,temp_pose.chain_begin(1),temp_pose.chain_begin(1)+96)

	#Now, import the TELSAM monomer with the right symmetry.
	url = f"https://files.rcsb.org/download/{sys.argv[2]}.pdb"
	pdb_text = requests.get(url).text
	with open("STEL.pdb","w") as file:
		file.write(pdb_text)
	cleanATOM("STEL.pdb")
	S_pose = pose_from_pdb('STEL.clean.pdb')

	#Now, it's time to superimpose the engineered TELSAM unit onto the one with the right symmetry.
	S_residues_to_superimpose = range(S_pose.chain_begin(1),S_pose.chain_begin(1)+E_pose.chain_end(1)-E_pose.chain_begin(1))
	E_residues_to_superimpose = range(E_pose.chain_begin(1),E_pose.chain_end(1))
	atom_map = AtomID_Map()
	initialize_atomid_map(atom_map, E_pose, AtomID())
	# Map CA atoms between residues
	for ER, SR in zip(E_residues_to_superimpose,S_residues_to_superimpose):
		E_atom = AtomID(E_pose.residue(ER).atom_index("CA"), ER)
		S_atom = AtomID(S_pose.residue(SR).atom_index("CA"), SR)
		atom_map.set(E_atom,S_atom)

	superimpose_pose(E_pose,S_pose,atom_map)
	E_pose.dump_pdb(f'./TELSetta/{sys.argv[1]}_in_{sys.argv[2]}_space.pdb')
	with open (f'./TELSetta/{sys.argv[1]}_in_{sys.argv[2]}_space.pdb','a') as file:
		with open("STEL.pdb") as s:
			for line in s:
				if "CRYST1" in line:
					file.write(line)

	os.remove("ETEL.pdb")
	os.remove("ETEL.clean.pdb")
	os.remove("STEL.pdb")
	os.remove("STEL.clean.pdb")

