#import and initialize PyRosetta and connect to the PyMOLMover server previously set up.

import sys
import requests
import os

import pyrosetta
from pyrosetta import *
from pyrosetta.toolbox import cleanATOM
from pyrosetta.rosetta.core.pose import append_subpose_to_pose
from pyrosetta.rosetta.core.pose import append_pose_to_pose
from pyrosetta.rosetta.core.pose import initialize_atomid_map
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.id import AtomID_Map_AtomID as AtomID_Map
from pyrosetta.rosetta.core.scoring import superimpose_pose
from pyrosetta.rosetta.protocols.grafting import delete_region
from pyrosetta.rosetta.core.scoring.dssp import Dssp

base = os.path.expanduser('~/TELSetta')

pyrosetta.init()

if sys.argv[1]=="1TEL":
	try:
		#Get 2QAR from the pdb (it's our most convenient engineered TELSAM because of the long helical linker at the end.)
		url = f"https://files.rcsb.org/download/2QAR.pdb"
		pdb_text = requests.get(url).text
		with open(os.path.join(base,"ETEL.pdb"),"w") as file:
			file.write(pdb_text)
		cleanATOM(os.path.join(base,"ETEL.pdb"))
		os.remove(os.path.join(base,"ETEL.pdb"))

		#Get 9DOC from the pdb (it has the proper space group that TELSAM usually fits into)
		url = f"https://files.rcsb.org/download/9DOC.pdb"
		pdb_text = requests.get(url).text
		with open(os.path.join(base,"STEL.pdb"),"w") as file:
			file.write(pdb_text)
		cleanATOM(os.path.join(base,"STEL.pdb"))
		#Don't remove the STEL.pdb yet because it has crystallographic information we want to get later.
		
		#Move 2QAR into the 9DOC asymmetric unit.
		TELSAM_in_9DOC = Pose()
		temp_pose = pose_from_pdb(os.path.join(base,'ETEL.clean.pdb'))
		os.remove(os.path.join(base,"ETEL.clean.pdb"))
		append_subpose_to_pose(TELSAM_in_9DOC,temp_pose,temp_pose.chain_begin(1),temp_pose.chain_begin(1)+76)
		S_pose = pose_from_pdb(os.path.join(base,'STEL.clean.pdb'))
		S_residues_to_superimpose = range(S_pose.chain_begin(1),S_pose.chain_begin(1)+TELSAM_in_9DOC.chain_end(1)-TELSAM_in_9DOC.chain_begin(1))
		E_residues_to_superimpose = range(TELSAM_in_9DOC.chain_begin(1),TELSAM_in_9DOC.chain_end(1))
		atom_map = AtomID_Map()
		initialize_atomid_map(atom_map, TELSAM_in_9DOC, AtomID())
		# Map CA atoms between residues
		for ER, SR in zip(E_residues_to_superimpose,S_residues_to_superimpose):
			E_atom = AtomID(TELSAM_in_9DOC.residue(ER).atom_index("CA"), ER)
			S_atom = AtomID(S_pose.residue(SR).atom_index("CA"), SR)
			atom_map.set(E_atom,S_atom)
		superimpose_pose(TELSAM_in_9DOC,S_pose,atom_map)

		#If a target protein was included, continue:
		if len(sys.argv)==3:
			try:
				client = Pose()
				url = f"https://files.rcsb.org/download/{sys.argv[2]}.pdb"
				pdb_text = requests.get(url).text
				with open(os.path.join(base,"temp_client.pdb"),"w") as file:
					file.write(pdb_text)
				cleanATOM(os.path.join(base,"temp_client.pdb"))
				temp_pose = Pose()
				pose_from_pdb(temp_pose,os.path.join(base,"temp_client.clean.pdb"))
				os.remove(os.path.join(base,"temp_client.pdb"))
				os.remove(os.path.join(base,"temp_client.clean.pdb"))
				#Extract first 4-aa helical region from target protein to fuse to TELSAM:
				dssp = Dssp(temp_pose)
				dssp.insert_ss_into_pose(temp_pose)
				ss_string = temp_pose.secstruct()
				first_helix = ss_string.find("HHHH")
				append_subpose_to_pose(client,temp_pose,temp_pose.chain_begin(1)+first_helix,temp_pose.chain_end(1))

				start_residue_to_superimpose = 2
				pdbs = []
				for resi in range(first_helix,first_helix+9):
					TELSAM = TELSAM_in_9DOC.clone()
					start_residue_to_superimpose += 1
					TELSAM_residues_to_superimpose = range(TELSAM.chain_end(1)-start_residue_to_superimpose,TELSAM.chain_end(1)-start_residue_to_superimpose+3)
					client_residues_to_superimpose = range(1,4)

					atom_map = AtomID_Map()
					initialize_atomid_map(atom_map, client, AtomID())

					# Map CA atoms between residues
					for CR, TR in zip(client_residues_to_superimpose,TELSAM_residues_to_superimpose):
						client_atom = AtomID(client.residue(CR).atom_index("CA"), CR)
						TELSAM_atom = AtomID(TELSAM.residue(TR).atom_index("CA"), TR)

						atom_map.set(client_atom,TELSAM_atom)

					superimpose_pose(client,TELSAM,atom_map)
					delete_region(TELSAM,TELSAM.chain_end(1)-start_residue_to_superimpose+1,TELSAM.chain_end(1))
					delete_region(client,client.chain_begin(1),client.chain_begin(1))
					append_pose_to_pose(TELSAM,client,new_chain=False)
					TELSAM.conformation().declare_chemical_bond(77-start_residue_to_superimpose,"C",77-start_residue_to_superimpose+1,"N")
					#relax here?
					#sf = get_score_function()
					#relax = FastRelax()
					#relax.set_scorefxn(sf)
					#relax.apply(TELSAM)
					TELSAM.dump_pdb(f'./TELSetta/{sys.argv[1]}_in_{sys.argv[2]}_space_{resi}.pdb')
					with open (os.path.join(base,f'{sys.argv[1]}_in_{sys.argv[2]}_space_{resi}.pdb'),'r') as file:
						pdb_sans_cryst = file.read()
					with open (os.path.join(base,f'{sys.argv[1]}_in_{sys.argv[2]}_space_{resi}.pdb'),'w') as file:
						with open(os.path.join(base,"STEL.pdb")) as s:
							for line in s:
								if "CRYST1" in line:
									file.write(line)
						file.write(pdb_sans_cryst)
					pdbs.append(os.path.join(base,f'{sys.argv[1]}_in_{sys.argv[2]}_space_{resi}.pdb'))

				with open(os.path.join(base,"pdb_list.txt"),"w") as file:
					for pdb in pdbs:
						file.write(f'{pdb}\n')

				os.remove(os.path.join(base,"STEL.pdb"))
				os.remove(os.path.join(base,"STEL.clean.pdb"))
		
			except Exception as e:
				print(e,file=sys.stderr)
				sys.exit()

	except Exception as e:
		print(e,file=sys.stderr)
		sys.exit()