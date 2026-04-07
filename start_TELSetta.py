

#import and initialize PyRosetta and connect to the PyMOLMover server previously set up.

import sys
import requests
import os
import re
import time

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

from pyrosetta import PyMOLMover
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.core.pose import addVirtualResAsRoot

base = os.path.expanduser('~/TELSetta')

pyrosetta.init("-crystal_refine -cryst::refinable_lattice -score_symm_complex -out:level 300 -out:file:scorefile scores.sc")

pmm = PyMOLMover()
pmm.keep_history(True)

def remake_TELSAM():
	if os.path.exists(os.path.join(base,f'TELSAM_in_9DOC.pdb')):
		os.remove(os.path.join(base,f'TELSAM_in_9DOC.pdb'))
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
	os.remove(os.path.join(base,"STEL.clean.pdb"))

	#Map CA atoms between residues
	S_residues_to_superimpose = range(S_pose.chain_begin(1),S_pose.chain_begin(1)+TELSAM_in_9DOC.total_residue())
	E_residues_to_superimpose = range(TELSAM_in_9DOC.chain_begin(1),TELSAM_in_9DOC.chain_end(1))
	atom_map = AtomID_Map()
	initialize_atomid_map(atom_map, TELSAM_in_9DOC, AtomID())
	for ER, SR in zip(E_residues_to_superimpose,S_residues_to_superimpose):
		E_atom = AtomID(TELSAM_in_9DOC.residue(ER).atom_index("CA"), ER)
		S_atom = AtomID(S_pose.residue(SR).atom_index("CA"), SR)
		atom_map.set(E_atom,S_atom)
	superimpose_pose(TELSAM_in_9DOC,S_pose,atom_map)

	helix_extender = Pose()
	append_subpose_to_pose(helix_extender,TELSAM_in_9DOC,TELSAM_in_9DOC.chain_end(1)-10,TELSAM_in_9DOC.chain_end(1))
	T_residues_to_superimpose = range(TELSAM_in_9DOC.chain_end(1)-3,TELSAM_in_9DOC.chain_end(1)+1)
	H_residues_to_superimpose = range(helix_extender.chain_begin(1),helix_extender.chain_begin(1)+4)
	helix_atom_map = AtomID_Map()
	initialize_atomid_map(helix_atom_map, helix_extender, AtomID())
	for HR, TR in zip(H_residues_to_superimpose,T_residues_to_superimpose):
		H_atom = AtomID(helix_extender.residue(HR).atom_index("CA"), HR)
		T_atom = AtomID(TELSAM_in_9DOC.residue(TR).atom_index("CA"), TR)
		helix_atom_map.set(H_atom,T_atom)
	superimpose_pose(helix_extender,TELSAM_in_9DOC,helix_atom_map)
	#Delete overlap
	delete_region(TELSAM_in_9DOC,TELSAM_in_9DOC.chain_end(1)-3,TELSAM_in_9DOC.chain_end(1))
	#Fuse
	append_pose_to_pose(TELSAM_in_9DOC,helix_extender,new_chain=False)
	TELSAM_in_9DOC.conformation().declare_chemical_bond(TELSAM_in_9DOC.chain_end(1)-helix_extender.total_residue(),"C",TELSAM_in_9DOC.chain_end(1)-helix_extender.total_residue()+1,"N")

	TELSAM_in_9DOC.dump_pdb(os.path.join(base,f'TELSAM_in_9DOC.pdb'))
	last_size = -1
	while True:
		if os.path.exists(os.path.join(base,f'TELSAM_in_9DOC.pdb')):
			size = os.path.getsize(os.path.join(base,f'TELSAM_in_9DOC.pdb'))
			if size == last_size:
				break
			last_size = size
		time.sleep(0.01)

	with open (os.path.join(base,f'TELSAM_in_9DOC.pdb'),'r') as file:
		pdb_sans_cryst = file.read()
	with open(os.path.join(base,'TELSAM_in_9DOC.pdb'),'w') as file:
		with open(os.path.join(base,"STEL.pdb")) as s:
			for line in s:
				if "CRYST1" in line:
					file.write(line)
					break
		file.write(pdb_sans_cryst)
	if os.path.exists(os.path.join(base,"STEL.pdb")):
		os.remove(os.path.join(base,"STEL.pdb"))
	print(f'Remade TELSAM. Stored in: {base}')

if sys.argv[1]=="1TEL":
	if "remake_TELSAM" in sys.argv:
		remake_TELSAM()
	#Get pre-made .pdb:
	if os.path.exists(os.path.join(base,f'TELSAM_in_9DOC.pdb')):
		try:
			TELSAM_in_9DOC = pose_from_file(os.path.join(base,f'TELSAM_in_9DOC.pdb'))
		except Exception:
			remake_TELSAM()
			TELSAM_in_9DOC = pose_from_file(os.path.join(base,f'TELSAM_in_9DOC.pdb'))
	else:
		remake_TELSAM()
		TELSAM_in_9DOC = pose_from_file(os.path.join(base,f'TELSAM_in_9DOC.pdb'))

	#If a target protein was included, fuse it to TELSAM:
	if len(sys.argv)>=3:
		try:
			client = Pose()
			if ".pdb" not in sys.argv[2]:
				client_ID = sys.argv[2].removesuffix(".pdb")
				url = f"https://files.rcsb.org/download/{sys.argv[2]}.pdb"
				pdb_text = requests.get(url).text
				with open(os.path.join(base,f"{sys.argv[2]}.pdb"),"w") as file:
					file.write(pdb_text)
			cleanATOM(os.path.join(base,f"{sys.argv[2]}.pdb"))
			temp_pose = Pose()
			pose_from_pdb(temp_pose,os.path.join(base,f"{sys.argv[2]}.clean.pdb"))
			os.remove(os.path.join(base,f"{sys.argv[2]}.pdb"))
			os.remove(os.path.join(base,f"{sys.argv[2]}.clean.pdb"))
			
			#Extract first 4-aa helical region from target protein to fuse to TELSAM:
			dssp = Dssp(temp_pose)
			dssp.insert_ss_into_pose(temp_pose)
			ss_string = temp_pose.secstruct()
			first_helix = ss_string.find("HHHHH")
			append_subpose_to_pose(client,temp_pose,temp_pose.chain_begin(1)+first_helix,temp_pose.chain_end(1))
			
			#Align the two helices:
			start_residue_to_superimpose = 1
			pdbs = []
			for resi in range(first_helix,first_helix+14):
				TELSAM = TELSAM_in_9DOC.clone()
				start_residue_to_superimpose += 1
				TELSAM_residues_to_superimpose = range(TELSAM.chain_end(1)-start_residue_to_superimpose,TELSAM.chain_end(1)-start_residue_to_superimpose+3)
				client_residues_to_superimpose = range(1,4)
				atom_map = AtomID_Map()
				initialize_atomid_map(atom_map, client, AtomID())

				#Map CA atoms between residues
				for CR, TR in zip(client_residues_to_superimpose,TELSAM_residues_to_superimpose):
					client_atom = AtomID(client.residue(CR).atom_index("CA"), CR)
					TELSAM_atom = AtomID(TELSAM.residue(TR).atom_index("CA"), TR)
					atom_map.set(client_atom,TELSAM_atom)
				superimpose_pose(client,TELSAM,atom_map)

				#Delete overlap
				delete_region(TELSAM,TELSAM.chain_end(1)-start_residue_to_superimpose,TELSAM.chain_end(1))
				#Fuse
				append_pose_to_pose(TELSAM,client,new_chain=False)
				TELSAM.conformation().declare_chemical_bond(TELSAM.chain_end(1)-client.total_residue()-start_residue_to_superimpose,"C",TELSAM.chain_end(1)-client.total_residue()-start_residue_to_superimpose+1,"N")
				#Save
				TELSAM.dump_pdb(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}.pdb'))
				
				if not "skip_symmetry" in sys.argv:
					#Test different unit cell sizes to dock polymers:
					for ucdelta in range(0,40):
						with open(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}_{ucdelta}.pdb'), 'w') as file:
							with open(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}.pdb'), 'r') as s:	
								pdb_sans_cryst = s.read()
								for line in s:
									if "CRYST1" in line:
										p = re.compile(r'\d+\.\d+')
										cryst1_vals = p.findall(line)

										a, b, c = [float(x) for x in cryst1_vals[0:3]]
										alpha, beta, gamma = [float(x) for x in cryst1_vals[3:6]]

										#Shrink unit cell
										a -= 00
										b -= 00
										#c -= 0

										spacegroup = line[55:66].strip()
										z_value = line[66:].strip()
										new_line = (
											f"CRYST1"
											f"{a:9.3f}{b:9.3f}{c:9.3f}"
											f"{alpha:7.2f}{beta:7.2f}{gamma:7.2f} "
											f"{spacegroup:<11}"
											f"{z_value:>4}\n"
										)
										break
							file.write(new_line)
							file.write(pdb_sans_cryst)

				#Make new pose from fusion for clarity
				symm_pose = pose_from_pdb(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}.pdb'))
				
				#Constraints on 1TEL
				movemap = MoveMap()
				movemap.set_bb(False)
				movemap.set_chi(False)

				# Allow movement only in the flexible region
				print(f'FLEXIBLE RANGE: {(TELSAM.chain_end(1)-client.total_residue()-start_residue_to_superimpose,TELSAM.chain_end(1))}')
				for i in range(TELSAM.chain_end(1)-client.total_residue()-start_residue_to_superimpose,TELSAM.chain_end(1)):
					movemap.set_bb(i, True)
					movemap.set_chi(i, True)

				"""
				constraints = pyrosetta.rosetta.core.scoring.constraints.ConstraintSet()
				addVirtualResAsRoot(symm_pose)
				root_atom = AtomID(1,symm_pose.total_residue())
				for i in range(1,TELSAM.chain_end(1)-client.total_residue()-start_residue_to_superimpose):
					atom_id = AtomID(symm_pose.residue(i).atom_index("CA"),i)
					xyz = symm_pose.residue(i).xyz("CA")
					constraint = pyrosetta.rosetta.core.scoring.constraints.CoordinateConstraint(
						atom_id,
						root_atom,
						xyz,
						pyrosetta.rosetta.core.scoring.func.HarmonicFunc(0.0,0.001)
					)
					constraints.add_constraint(constraint)
				symm_pose.constraint_set(constraints)
				"""
				
				#Score
				sf = get_score_function()
				#sf.set_weight(rosetta.core.scoring.coordinate_constraint, 1.0)
				score = sf(symm_pose)
				print(f'SCORE!!!!::::{score}')

				if not "skip_symmetry" in sys.argv:
					#Symmetrize
					interaction_shell_size = 30
					rosetta.basic.options.set_real_option("cryst:interaction_shell", interaction_shell_size)
					makesym = SetupForSymmetryMover("CRYST1")
					makesym.apply(symm_pose)
					symm_pose.dump_pdb(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_symm_{resi}_shellSize_{interaction_shell_size}.pdb'))

				#Relax
				relax = FastRelax()
				relax.set_movemap(movemap)
				relax.set_scorefxn(sf)
				relax.apply(symm_pose)
				if not "skip_symmetry" in sys.argv:
					symm_pose.dump_pdb(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_symmFastRelax_{resi}_shellSize_{interaction_shell_size}.pdb'))
				else:
					symm_pose.dump_pdb(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_FastRelax_{resi}.pdb'))

				#Show in PyMOL
				pmm.apply(symm_pose)
	
		except Exception as e:
			print(e,file=sys.stderr)
			sys.exit()