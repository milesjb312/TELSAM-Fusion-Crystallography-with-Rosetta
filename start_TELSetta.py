

#import and initialize PyRosetta and connect to the PyMOLMover server previously set up.

import sys
import requests
import os
import re
import time
import math

import pyrosetta
from pyrosetta import *
from pyrosetta.toolbox import cleanATOM
from pyrosetta.toolbox.mutants import mutate_residue
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

from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import InitializeFromCommandline, RestrictToRepacking
from pyrosetta.rosetta.protocols.minimization_packing import MinMover

from pyrosetta.rosetta.numeric import xyzMatrix_double_t, xyzVector_double_t

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
	mutate_residue(temp_pose,67,"V",5)
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

				####################SETUP FOR REFINEMENT##########################
				#Get score_function
				sf = get_score_function()
				#Relax mover
				relax = FastRelax()
				relax.set_scorefxn(sf)

				tf = TaskFactory()
				tf.push_back(InitializeFromCommandline())
				tf.push_back(RestrictToRepacking())

				packer = PackRotamersMover(sf)
				packer.task_factory(tf)

				min_mover = MinMover()
				min_mover.score_function(sf)
				min_mover.min_type("lbfgs_armijo_nonmonotone")

				#######################MOVEMAP (Must be re-setup after pose is symmetrized)#################
				#It may be okay to just have it here and then call the min_mover.movemap and relax.set_movemap functions later.
				movemap = MoveMap()
				movemap.set_bb(False)
				movemap.set_chi(False)
				for i in range(TELSAM.chain_end(1)-client.total_residue()-start_residue_to_superimpose,TELSAM.chain_end(1)):
					movemap.set_bb(i, True)
					movemap.set_chi(i, True)
				##############################################################################################

				#Refine gently
				min_mover.movemap(movemap)
				packer.apply(TELSAM)
				min_mover.apply(TELSAM)

				#Save
				TELSAM.dump_pdb(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}.pdb'))

				#Filter:
				score = sf(TELSAM)
				if score<5000:
					if not "skip_symmetry" in sys.argv:
						#Determine largest unit cell:
						with open(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}.pdb'), 'r') as file:
							lines = iter(file)
							furthest_x = 0
							for line in lines:
								if 'CRYST1' in line:
									p = re.compile(r'\d+\.\d+')
									a = float(p.search(line).group())

								if 'ATOM' in line:
									x_coord = float(line[31:39].strip())
									if x_coord>furthest_x:
										furthest_x = x_coord

						######################################## DOCK POLYMERS! ############################
						#Test different unit cell sizes:
						min_score = 100000
						if os.path.exists(os.path.join(base,f'ucdelta_scores.txt')):
							os.remove(os.path.join(base,f'ucdelta_scores.txt'))
						for ucdelta in range(int(a-x_coord*2),int(a-x_coord*2+5)):
							with open(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}_{ucdelta}.pdb'), 'w') as file:
								with open(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}.pdb'), 'r') as s:
									for line in s:
										if "CRYST1" in line:
											p = re.compile(r'\d+\.\d+')
											cryst1_vals = p.findall(line)
											a, b, c = [float(x) for x in cryst1_vals[0:3]]
											alpha, beta, gamma = [float(x) for x in cryst1_vals[3:6]]

											#Shrink unit cell
											a -= ucdelta
											b -= ucdelta
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
											file.write(new_line)
										else:
											file.write(line)

							#Symmetrize
							#Make new pose from fusion for clarity
							symm_pose = pose_from_file(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}_{ucdelta}.pdb'))
							interaction_shell_size = 30
							rosetta.basic.options.set_real_option("cryst:interaction_shell", interaction_shell_size)
							makesym = SetupForSymmetryMover("CRYST1")
							makesym.apply(symm_pose)
							
							#Refine gently
							min_mover.movemap(movemap)
							packer.apply(symm_pose)
							min_mover.apply(symm_pose)

							#relax.set_movemap(movemap)
							#relax.apply(symm_pose)
							#symm_pose.dump_pdb(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_symm_{resi}_shellSize_{interaction_shell_size}.pdb'))

							#Return score
							score = sf(symm_pose)
							if score<min_score:
								min_score = score
							if score<min_score*2:
								passing_pdb = os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}_{ucdelta}_rpkmin.pdb')
								with open(os.path.join(base,f'ucdelta_scores.txt'),'a') as scores:
									scores.write(f'file: {passing_pdb}\n')
									scores.write(f'score: {score}\n')
								symm_pose.dump_pdb(passing_pdb)
							else:
								break

							#Show in PyMOL
							pmm.apply(symm_pose)

							##################################### ROTATE UNIT CELL! ###############################
							if os.path.exists(os.path.join(base,f'deg_scores.txt')):
								os.remove(os.path.join(base,f'deg_scores.txt'))
							for deg in range(1,61):
								#Make new pose from fusion for clarity
								symm_pose = pose_from_file(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}_{ucdelta}.pdb'))
								theta = math.radians(deg)
								R = xyzMatrix_double_t()
								R.xx = math.cos(theta)
								R.xy = -math.sin(theta)
								R.xz = 0.0
								R.yx = math.sin(theta)
								R.yy = math.cos(theta)
								R.yz = 0.0
								R.zx = 0.0
								R.zy = 0.0
								R.zz = 1.0

								# No translation
								v = xyzVector_double_t(0.0, 0.0, 0.0)
								symm_pose.apply_transform_Rx_plus_v(R, v)
								"""
								pose.apply_transform_Rx_plus_v(R, [0,0,0])
								with open(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}_{ucdelta}_{deg}.pdb'), 'w') as new_file:
									with open(passing_pdb,'r') as file:
										lines = iter(file)
										for line in lines:
											if line.startswith("ATOM"):
												prefix = line[:30]
												suffix = line[54:]
												x = float(line[30:38])
												y = float(line[38:46])
												z = float(line[46-54])
												new_x = x*math.cos(math.radians(deg)) - y*math.sin(math.radians(deg))
												new_y = x*math.sin(math.radians(deg)) + y*math.cos(math.radians(deg))
												line = (
												f'{prefix}'
												f'{new_x:8.3f}'
												f'{new_y:8.3f}'
												f'{z:8.3f}'
												f'{suffix:}\n'
											)
											new_file.write(line)
								"""
								#Symmetrize
								interaction_shell_size = 30
								rosetta.basic.options.set_real_option("cryst:interaction_shell", interaction_shell_size)
								makesym = SetupForSymmetryMover("CRYST1")
								makesym.apply(symm_pose)
								
								#Refine gently
								min_mover.movemap(movemap)
								packer.apply(symm_pose)
								min_mover.apply(symm_pose)

								#relax.set_movemap(movemap)
								#relax.apply(symm_pose)
								#symm_pose.dump_pdb(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_symm_{resi}_shellSize_{interaction_shell_size}.pdb'))

								#Return score
								score = sf(symm_pose)
								passing_pdb = os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}_{ucdelta}_{deg}_rpkmin.pdb')
								with open(os.path.join(base,f'deg_scores.txt'),'a') as scores:
									scores.write(f'file: {passing_pdb}\n')
									scores.write(f'score: {score}\n')
								symm_pose.dump_pdb(passing_pdb)

								#Show in PyMOL
								pmm.apply(symm_pose)
							os.remove(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}_{ucdelta}.pdb'))

						"""
						file_scores = {}
						opti_set = []
						if os.path.exists(os.path.join(base,f'ucdelta_scores')):
							with open(os.path.join(base,f'ucdelta_scores'),'r') as scores:						
								lines = iter(scores)
								for line in lines:
									line.strip()
									if line.startswith('file: '):
										file = line.removeprefix('file: ').strip()
										score_line = next(lines).strip()
										if score_line.startswith('score: '):
											score = float(score_line.removeprefix('score: ').strip())
											file_scores[file] = score
							smallest_score = min(file_scores.values())
							for file in file_scores:
								#Optimize the file with the smallest score
								if file_scores[file] == smallest_score:
									opti_set.append(file)
								#Optimize the files with scores slightly higher than the file with the smallest score.
								#This may improve the minimizer's options by increasing the number of contacts it tries to optimize.
								if smallest_score*1.005<file_scores[file]<smallest_score*1.1:
									opti_set.append(file)
						"""

					else:
						asymm_pose = pose_from_file(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}.pdb'))
						relax.set_movemap(movemap)

						#Relax
						relax.apply(asymm_pose)
						asymm_pose.dump_pdb(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_FastRelax_{resi}.pdb'))

						#Return score
						score = sf(symm_pose)
						print(f'SCORE!!!!::::{score}')

						#Show in PyMOL
						pmm.apply(symm_pose)

		except Exception as e:
			print(e,file=sys.stderr)
			sys.exit()