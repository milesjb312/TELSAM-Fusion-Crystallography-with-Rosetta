#import and initialize PyRosetta and connect to the PyMOLMover server previously set up.

import sys
import requests
import os

import pyrosetta
from pyrosetta import *
from pyrosetta import PyMOLMover
from pyrosetta.toolbox import cleanATOM
from pyrosetta.rosetta.core.import_pose import pose_from_pdbstring
from pyrosetta.rosetta.core.pose import append_subpose_to_pose
from pyrosetta.rosetta.core.pose import append_pose_to_pose
from pyrosetta.rosetta.core.pose import initialize_atomid_map
from pyrosetta.rosetta.core.pose.symmetry import sym_dof_names
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.id import AtomID_Map_AtomID as AtomID_Map
from pyrosetta.rosetta.core.scoring import superimpose_pose
from pyrosetta.rosetta.protocols.grafting import delete_region
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.core.scoring.dssp import Dssp

pyrosetta.init()

pmm = PyMOLMover()

#This line makes it so that all poses sent to PyMol are kept as part of a sequence of poses.
pmm.keep_history(True)

TELSAM_from_pdb = Pose()
client_protein_from_pdb = Pose()

client_and_guest_loaded = False

if sys.argv[1]=="1TEL":
	try:
		url = f"https://files.rcsb.org/download/9DOC.pdb"
		pdb_text = requests.get(url).text
		with open("temp.pdb","w") as file:
			file.write(pdb_text)
		cleanATOM("temp.pdb")
		pose_from_pdb(TELSAM_from_pdb,"temp.clean.pdb")
		os.remove("temp.pdb")
		os.remove("temp.clean.pdb")
		print(f'len(sys.argv)={len(sys.argv)}')
		if len(sys.argv)==3:
			try:
				url = f"https://files.rcsb.org/download/{sys.argv[2]}.pdb"
				pdb_text = requests.get(url).text
				with open("temp_client.pdb","w") as file:
					file.write(pdb_text)
				cleanATOM("temp_client.pdb")
				pose_from_pdb(client_protein_from_pdb,"temp_client.clean.pdb")
				os.remove("temp_client.pdb")
				os.remove("temp_client.clean.pdb")
				client_and_guest_loaded = True
			except Exception as e:
				print(e)
	except Exception as e:
		print(e)

TELSAM_mastercopy = Pose()
client_protein_mastercopy = Pose()

if client_and_guest_loaded:
	append_subpose_to_pose(TELSAM_mastercopy,TELSAM_from_pdb,TELSAM_from_pdb.chain_begin(1),TELSAM_from_pdb.chain_begin(1)+81)
	#Extract first 4-aa helical region from target protein to fuse to TELSAM:
	dssp = Dssp(client_protein_from_pdb)
	dssp.insert_ss_into_pose(client_protein_from_pdb)
	ss_string = client_protein_from_pdb.secstruct()
	first_helix = ss_string.find("HHHH")
	append_subpose_to_pose(client_protein_mastercopy,client_protein_from_pdb,client_protein_from_pdb.chain_begin(1)+first_helix,client_protein_from_pdb.chain_end(1))
	start_residue_to_superimpose = 2
	for resi in range(first_helix,first_helix+9):
		start_residue_to_superimpose += 1
		TELSAM = TELSAM_mastercopy.clone()
		client_protein = client_protein_mastercopy.clone()
		TELSAM_residues_to_superimpose = range(TELSAM.total_residue()-start_residue_to_superimpose,TELSAM.total_residue()-start_residue_to_superimpose+3)
		client_residues_to_superimpose = range(1,4)

		atom_map = AtomID_Map()
		initialize_atomid_map(atom_map, client_protein, AtomID())

		# Map CA atoms between residues
		for CR, TR in zip(client_residues_to_superimpose,TELSAM_residues_to_superimpose):
			client_atom = AtomID(client_protein.residue(CR).atom_index("CA"), CR)
			TELSAM_atom = AtomID(TELSAM.residue(TR).atom_index("CA"), TR)

			atom_map.set(client_atom,TELSAM_atom)

		superimpose_pose(client_protein,TELSAM,atom_map)
		delete_region(TELSAM,TELSAM.chain_end(1)-start_residue_to_superimpose+1,TELSAM.chain_end(1))
		delete_region(client_protein,client_protein.chain_begin(1),client_protein.chain_begin(1))
		append_pose_to_pose(TELSAM,client_protein,new_chain=False)
		TELSAM.conformation().declare_chemical_bond(81-start_residue_to_superimpose+1,"C",82-start_residue_to_superimpose+1,"N")		
		symm_mover = SetupForSymmetryMover("./TELSetta/9DOC_2.symm")
		symm_mover.apply(TELSAM)
		print(f'Symmetric pose: {is_symmetric(TELSAM)}')
		sf = get_score_function()
		relax = FastRelax()
		relax.set_scorefxn(sf)
		#relax.apply(TELSAM)

		pmm.apply(TELSAM)


