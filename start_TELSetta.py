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
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.scoring import get_score_function
from pyrosetta.rosetta.core.scoring import fa_rep

from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import InitializeFromCommandline, RestrictToRepacking
from pyrosetta.rosetta.protocols.minimization_packing import MinMover

from pyrosetta.rosetta.numeric import xyzMatrix_double_t, xyzVector_double_t

from matplotlib import pyplot as plt

fig = plt.figure()
base = os.path.expanduser('~/TELSetta')

pyrosetta.init("-crystal_refine -cryst::refinable_lattice -score_symm_complex -out:level 300 -out:file:scorefile scores.sc")

pmm = PyMOLMover()
pmm.keep_history(True)
energies_vs_ucab_vs_deg = {'linker':[],'energy':[],'ucab':[],'deg':[]}

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

	################### MOVE 2QAR INTO 9DOC ASYMMETRIC UNIT ##############################
	#Grab a portion of 2QAR
	TELSAM_in_9DOC = Pose()
	temp_pose = pose_from_pdb(os.path.join(base,'ETEL.clean.pdb'))
	mutate_residue(temp_pose,67,"V",5)
	os.remove(os.path.join(base,"ETEL.clean.pdb"))
	append_subpose_to_pose(TELSAM_in_9DOC,temp_pose,temp_pose.chain_begin(1),temp_pose.chain_begin(1)+76)
	S_pose = pose_from_pdb(os.path.join(base,'STEL.clean.pdb'))
	os.remove(os.path.join(base,"STEL.clean.pdb"))

	#Superimpose CA atoms between 2QAR and 9DOC
	S_residues_to_superimpose = range(S_pose.chain_begin(1),S_pose.chain_begin(1)+TELSAM_in_9DOC.total_residue())
	E_residues_to_superimpose = range(TELSAM_in_9DOC.chain_begin(1),TELSAM_in_9DOC.chain_end(1))
	atom_map = AtomID_Map()
	initialize_atomid_map(atom_map, TELSAM_in_9DOC, AtomID())
	for ER, SR in zip(E_residues_to_superimpose,S_residues_to_superimpose):
		E_atom = AtomID(TELSAM_in_9DOC.residue(ER).atom_index("CA"), ER)
		S_atom = AtomID(S_pose.residue(SR).atom_index("CA"), SR)
		atom_map.set(E_atom,S_atom)
	superimpose_pose(TELSAM_in_9DOC,S_pose,atom_map)

	####################################### EXTEND HELIX ###############################################
	helix_extender = Pose()
	#Grab the last 11 residues in TELSAM:
	append_subpose_to_pose(helix_extender,TELSAM_in_9DOC,TELSAM_in_9DOC.chain_end(1)-10,TELSAM_in_9DOC.chain_end(1))
	#Align those residues to the end of the helix over 4 amino acids (effectively copying the helix and shifting it over on top of itself)
	T_residues_to_superimpose = range(TELSAM_in_9DOC.chain_end(1)-3,TELSAM_in_9DOC.chain_end(1)+1)
	H_residues_to_superimpose = range(helix_extender.chain_begin(1),helix_extender.chain_begin(1)+4)
	helix_atom_map = AtomID_Map()
	initialize_atomid_map(helix_atom_map, helix_extender, AtomID())
	for HR, TR in zip(H_residues_to_superimpose,T_residues_to_superimpose):
		H_atom = AtomID(helix_extender.residue(HR).atom_index("CA"), HR)
		T_atom = AtomID(TELSAM_in_9DOC.residue(TR).atom_index("CA"), TR)
		helix_atom_map.set(H_atom,T_atom)
	superimpose_pose(helix_extender,TELSAM_in_9DOC,helix_atom_map)

	#Delete 4-aa overlap
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

def change_cell(read_file,write_file,wa=None,wb=None,wdc=None):
	with open(write_file, 'w') as file:
		with open(read_file, 'r') as s:
			for line in s:
				if "CRYST1" in line:
					p = re.compile(r'\d+\.\d+')
					cryst1_vals = p.findall(line)
					a, b, c = [float(x) for x in cryst1_vals[0:3]]
					if wa!=None:
						a = wa
					if wb!=None:
						b = wb
					if wdc!=None:
						c += wdc
					alpha, beta, gamma = [float(x) for x in cryst1_vals[3:6]]
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

def change_deg(current_ucab_pdb,deg):
	symm_pose = pose_from_file(current_ucab_pdb)
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
	v = xyzVector_double_t(0.0, 0.0, 0.0) #No translation
	symm_pose.apply_transform_Rx_plus_v(R, v)
	makesym.apply(symm_pose)
	return symm_pose

def scilter(symm_pose,min_score,er_cutoff,min_score_pdb,current_pdb,last_pdb,scored_file_name):
	"""scilter scores the symmetric pose, compares the score to the min_score*er_cutoff, and determines whether to update the min_score_pdb or not 
	with the current PDB. It updates the scores_file with the final min_score_pdb of any given setting. Finally, it removes unneeded files.
	It returns a tuple with the min_score,min_score_pdb,and a boolean reporting whether the current PDB exceeded the er_cutoff (True if exceeded).
	"""
	score = sf(symm_pose)
	print(f'score: {score}, min_score: {min_score}, min_score_pdb: {min_score_pdb}, current_pdb: {current_pdb}')
	if 0<score<=min_score*er_cutoff or score<min_score*0.5<0:
		if score<min_score:
			min_score = score
			min_score_pdb = current_pdb
			pmm.apply(symm_pose)
		if last_pdb!=None:
			if os.path.exists(last_pdb) and last_pdb!= min_score_pdb:
				os.remove(last_pdb)
				print(f'removed previous pdb: {last_pdb}')
		return (min_score,min_score_pdb,False)
	else:
		if min_score_pdb != None:
			with open(os.path.join(base,f'scores_file.txt'),'a') as scores:
				scores.write(f'{scored_file_name}: {current_pdb}\n')
				scores.write(f'Score: {score}\n')
				specs = str(current_pdb).split("_")
				for spec in range(len(specs)):
					if str(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}')) in specs[spec]:
						if len(specs)==spec+3:
							energies_vs_ucab_vs_deg['linker'] = specs[spec+1]
							energies_vs_ucab_vs_deg["ucab"] = specs[spec+2]
							energies_vs_ucab_vs_deg["deg"] = specs[spec+3]
							energies_vs_ucab_vs_deg['energy'] = score
			if last_pdb!=None:
				if os.path.exists(last_pdb) and last_pdb!= min_score_pdb:
					os.remove(last_pdb)
		return (min_score,min_score_pdb,True)
	
def chart(linker):
	"""chart is used to make a 3d graph of the relationship between energy, unit cell ab, and degree for each linker length variant.
	"""
	ax = fig.add_subplot(projection='3d')
	aboi = [ucab for ucab, l in zip(energies_vs_ucab_vs_deg['ucab'],energies_vs_ucab_vs_deg['linker']) if l == linker]
	doi = [deg for deg, l in zip(energies_vs_ucab_vs_deg['deg'],energies_vs_ucab_vs_deg['linker']) if l == linker]
	eoi = [energy for energy, l in zip(energies_vs_ucab_vs_deg['energy'],energies_vs_ucab_vs_deg['linker']) if l == linker]
	ax.scatter(aboi,doi,eoi)
	ax.set_title(f'Energies of AB:Degree Combinations for {sys.argv[1]}--{sys.argv[2]}_{linker}')
	ax.set_xlabel('Unit Cell AB Length (Angstroms)')
	ax.set_ylabel('Degree of Polymer Rotation (Degrees)')
	ax.set_zlabel('Energy (REU)')
	plt.show()

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

	if len(sys.argv)>=3:
		########################## CREATE TELSAM FUSION! ######################################
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
			if os.path.exists(os.path.join(base,f'scores_file.txt')):
				os.remove(os.path.join(base,f'scores_file.txt'))
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
				TELSAM.conformation().declare_chemical_bond(TELSAM.chain_end(1)-client.total_residue(),"C",TELSAM.chain_end(1)-client.total_residue()+1,"N")

				#################### SETUP FOR REFINEMENT ##########################
				interaction_shell_size = 30
				rosetta.basic.options.set_real_option("cryst:interaction_shell", interaction_shell_size)
				makesym = SetupForSymmetryMover("CRYST1")

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

				####################### MOVEMAP (Must be re-setup after pose is symmetrized) #################
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

				#relax.set_movemap(movemap)
				#relax.apply(symm_pose)
				#symm_pose.dump_pdb(os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_symm_{resi}_shellSize_{interaction_shell_size}.pdb'))

				#Save
				current_linker_pdb = os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}.pdb')
				TELSAM.dump_pdb(current_linker_pdb)

				#Filter:
				score = sf(TELSAM)
				if score<10000:
					#Record base score of linker variant:
					with open(os.path.join(base,f'scores_file.txt'),'a') as scores:
						scores.write(f'Linker file: {current_linker_pdb}\n')
						scores.write(f'score: {score}\n')
					
					######################################## DOCK POLYMERS! ##################################
					#Determine largest unit cell:
					with open(current_linker_pdb, 'r') as file:
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

					#Test different unit cell sizes:
					min_score = 200000
					min_score_pdb = None
					for ucab in range(int(furthest_x*2),int(furthest_x+30),-1):
						current_ucab_pdb = os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}_{ucab}.pdb')
						change_cell(current_linker_pdb,current_ucab_pdb,wa=ucab,wb=ucab)
						symm_pose = pose_from_file(current_ucab_pdb)
						makesym.apply(symm_pose)
						min_score,min_score_pdb,exceeded = scilter(symm_pose,min_score,1.1,min_score_pdb,current_ucab_pdb,os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}_{ucab+1}.pdb'),'Unit Cell AB File')
						if exceeded:
							if min_score_pdb != None:
								min_ucab_pdb = min_score_pdb
								ucab1 = int(str(min_ucab_pdb).split("_")[-1].removesuffix(".pdb"))
							break
					##################################### ROTATE UNIT CELL AND CONTINUE DOCKING! ###############################
					if min_ucab_pdb!=None:
						ucab = ucab1
						for deg in range(1,21):
							current_deg_pdb = os.path.join(base,f'{min_ucab_pdb.removesuffix(".pdb")}_{deg}.pdb')
							for ucab2 in range(ucab,int(furthest_x+30),-1):
								current_ucab2_deg_pdb = os.path.join(base,f'{min_ucab_pdb.removesuffix(f"{ucab1}.pdb")}{ucab2}_{deg}.pdb')
								change_cell(min_ucab_pdb,current_ucab2_deg_pdb,wa=ucab2,wb=ucab2)
								symm_pose = change_deg(current_ucab2_deg_pdb,deg)
								pmm.apply(symm_pose)
								min_score,min_score_pdb,exceeded = scilter(symm_pose,min_score,1.1,min_score_pdb,current_ucab2_deg_pdb,os.path.join(base,f'{min_ucab_pdb.removesuffix(f"{ucab1}.pdb")}_{ucab2+1}_{deg}.pdb'),'Temp Unit Cell AB2 Deg File')
								if exceeded:
									ucab = ucab2
									break
							min_score,min_score_pdb,exceeded = scilter(symm_pose,min_score,3,min_score_pdb,current_ucab2_deg_pdb,None,'Unit Cell AB2 Deg File')
							#EITHER leave this alone (meaning you will always test all the degrees in the range) OR calculate the angle between the linker vector and
							#the furthest edge of the linker, because that angle has to be the starter angle for whichever degree check is opposite
							#that angle (since proceeding by just 1 degree in the opposing direction will cause the edge of the linker to continually bump and increase
							#in energy)

						ucab = ucab1
						min_score_pdb = min_ucab_pdb
						for deg in range(-1,-21,-1):
							current_deg_pdb = os.path.join(base,f'{min_ucab_pdb.removesuffix(".pdb")}_{deg}.pdb')
							for ucab2 in range(ucab,int(furthest_x+30),-1):
								current_ucab2_deg_pdb = os.path.join(base,f'{min_ucab_pdb.removesuffix(f"{ucab1}.pdb")}{ucab2}_{deg}.pdb')
								change_cell(min_ucab_pdb,current_ucab2_deg_pdb,wa=ucab2,wb=ucab2)
								symm_pose = change_deg(current_ucab2_deg_pdb,deg)
								pmm.apply(symm_pose)
								min_score,min_score_pdb,exceeded = scilter(symm_pose,min_score,1.1,min_score_pdb,current_ucab2_deg_pdb,os.path.join(base,f'{min_ucab_pdb.removesuffix(f"{ucab1}.pdb")}_{ucab2+1}_{deg}.pdb'),'Temp Unit Cell AB2 Deg File')
								if exceeded:
									ucab = ucab2 #This assumes that client proteins are not too concave (pointing in toward the polymer vector)
									break
							min_score,min_score_pdb,exceeded = scilter(symm_pose,min_score,3,min_score_pdb,current_ucab2_deg_pdb,None,'Unit Cell AB2 Deg File')
						chart(resi)
					"""
					##################################### CHANGE HELICAL RISE! #####################################
					#Test different unit cell sizes:
					min_score = 200000
					for ucc in range(1,15):
						change_cell(passing_pdb,os.path.join(base,f'{passing_pdb.removesuffix(".pdb")}_{ucc}.pdb'),wc=ucc)
						symm_pose = pose_from_file(os.path.join(base,f'{passing_pdb.removesuffix(".pdb")}_{ucc}.pdb'))
						makesym.apply(symm_pose)
						min_score,passing_pdb,finished = scilter(symm_pose,min_score,passing_pdb,min_ucab_pdb,os.path.join(base,f'{sys.argv[1]}--{sys.argv[2]}_{resi}_{ucab+1}.pdb'),'Unit Cell C File')
						if finished:
							break
					"""

		except Exception as e:
			print(e,file=sys.stderr)
			sys.exit()