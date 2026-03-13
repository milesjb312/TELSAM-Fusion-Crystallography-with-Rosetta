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
	pyrosetta.init()
	url = f"https://files.rcsb.org/download/{sys.argv[1]}.pdb"
	pdb_text = requests.get(url).text
	with open("./TELSetta/temp_prot.pdb","w") as file:
		file.write(pdb_text)
	cleanATOM("./TELSetta/temp_prot.pdb")
	temp_pose = pose_from_pdb('./TELSetta/temp_prot.clean.pdb')
	final_pose = Pose()
	append_subpose_to_pose(final_pose,temp_pose,temp_pose.chain_begin(1),temp_pose.chain_end(1))
	final_pose.dump_pdb(f'./TELSetta/{sys.argv[1]}.pdb')
	#os.remove("./TELSetta/temp_prot.pdb")
	os.remove("./TELSetta/temp_prot.clean.pdb")

if len(sys.argv) == 3:
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
	S_pose.dump_pdb(f'./TELSetta/{sys.argv[1]}.pdb')
	os.remove("ETEL.pdb")
	os.remove("ETEL.clean.pdb")
	os.remove("STEL.pdb")
	os.remove("STEL.clean.pdb")

