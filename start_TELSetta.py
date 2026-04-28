import sys
import requests
import os
import getopt
import re
import time
import math
import json
import numpy as np

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

from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import InitializeFromCommandline, RestrictToRepacking
from pyrosetta.rosetta.protocols.minimization_packing import MinMover

from pyrosetta.rosetta.numeric import xyzMatrix_double_t, xyzVector_double_t

from matplotlib import pyplot as plt

pyrosetta.init("-crystal_refine -cryst::refinable_lattice -score_symm_complex -out:level 300 -out:file:scorefile scores.sc")

class TELSetta:
	def __init__(self):
		self.TELSAM_version = "1TEL"
		self.client_pdb = None
		self.linker_variant = None
		self.unit_cell_ab = None
		self.degree_rotation = None
		self.remake_TELSAM_bool = False
		self.optimize = False
		self.centroids = True
		self.scores = {}
		try:
			optlist, args = getopt.getopt(sys.argv[1:], "t:c:l:u:d:r:o")
			for o, a in optlist:
				if o == '-t':
					if a != "1TEL":
						print(f'Error: TELSAM_version argument not supported: {a}')
						sys.exit(1)
					else:
						self.TELSAM_version = a
				elif o == '-c':
					self.client_pdb = a
				elif o == '-l':
					if a!="":
						self.linker_variant = int(a)
				elif o =='-u':
					if a!="":
						self.unit_cell_ab = float(a)
				elif o == '-d':
					if a!="":
						self.degree_rotation = int(float(a))
				elif o =='-r':
					self.remake_TELSAM_bool = a.upper() in ['T','TRUE',1]
				elif o == '-o':
					self.optimize = True
				else:
					print(f'Unhandled argument: {o}')
					sys.exit(1)
			if args:
				print(f'Unexpected arguments: {args}')
				sys.exit(1)
			if not self.optimize:
				self.centroids = False
		except getopt.GetoptError as err:
			print(err)
			sys.exit(1)

		self.base = os.path.join(os.path.expanduser('~/TELSetta'),str(self.linker_variant))
		os.makedirs(self.base,exist_ok=True)
		self.pmm = PyMOLMover()
		self.pmm.keep_history(True)
		self.energies_vs_ucab_vs_deg = {'linker':[],'energy':[],'ucab':[],'deg':[]}
		self.furthest_x = 0
		self.to_centroid = SwitchResidueTypeSetMover("centroid")
		self.to_fullatom = SwitchResidueTypeSetMover("fa_standard")
		self.validate_TELSAM()
		
	def remake_TELSAM(self):
		if os.path.exists(os.path.join(self.base,f'TELSAM_in_9DOC.pdb')):
			os.remove(os.path.join(self.base,f'TELSAM_in_9DOC.pdb'))
		#Get 2QAR from the pdb (it's our most convenient engineered TELSAM because of the long helical linker at the end.)
		url = f"https://files.rcsb.org/download/2QAR.pdb"
		pdb_text = requests.get(url).text
		with open(os.path.join(self.base,"ETEL.pdb"),"w") as file:
			file.write(pdb_text)
		cleanATOM(os.path.join(self.base,"ETEL.pdb"))
		os.remove(os.path.join(self.base,"ETEL.pdb"))

		#Get 9DOC from the pdb (it has the proper space group that TELSAM usually fits into)
		url = f"https://files.rcsb.org/download/9DOC.pdb"
		pdb_text = requests.get(url).text
		with open(os.path.join(self.base,"STEL.pdb"),"w") as file:
			file.write(pdb_text)
		cleanATOM(os.path.join(self.base,"STEL.pdb"))
		#Don't remove the STEL.pdb yet because it has crystallographic information we want to get later.

		################### MOVE 2QAR INTO 9DOC ASYMMETRIC UNIT ##############################
		#Grab a portion of 2QAR
		TELSAM_in_9DOC = Pose()
		temp_pose = pose_from_pdb(os.path.join(self.base,'ETEL.clean.pdb'))
		mutate_residue(temp_pose,67,"V",5)
		os.remove(os.path.join(self.base,"ETEL.clean.pdb"))
		append_subpose_to_pose(TELSAM_in_9DOC,temp_pose,temp_pose.chain_begin(2),temp_pose.chain_end(2))
		S_pose = pose_from_pdb(os.path.join(self.base,'STEL.clean.pdb'))
		os.remove(os.path.join(self.base,"STEL.clean.pdb"))

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

		"""
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
		"""

		TELSAM_in_9DOC.dump_pdb(os.path.join(self.base,f'TELSAM_in_9DOC.pdb'))
		last_size = -1
		while True:
			if os.path.exists(os.path.join(self.base,f'TELSAM_in_9DOC.pdb')):
				size = os.path.getsize(os.path.join(self.base,f'TELSAM_in_9DOC.pdb'))
				if size == last_size:
					break
				last_size = size
			time.sleep(0.01)

		with open (os.path.join(self.base,f'TELSAM_in_9DOC.pdb'),'r') as file:
			pdb_sans_cryst = file.read()
		with open(os.path.join(self.base,'TELSAM_in_9DOC.pdb'),'w') as file:
			with open(os.path.join(self.base,"STEL.pdb")) as s:
				for line in s:
					if "CRYST1" in line:
						file.write(line)
						break
			file.write(pdb_sans_cryst)
		if os.path.exists(os.path.join(self.base,"STEL.pdb")):
			os.remove(os.path.join(self.base,"STEL.pdb"))
		print(f'Remade TELSAM. Stored in: {self.base}')

	def change_cell(self,read_file,write_file,wa=None,wb=None,wdc=None):
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

	def change_deg(self,current_ucab_pdb,deg):
		symm_pose = pose_from_file(current_ucab_pdb)
		if self.centroids:
			self.to_centroid.apply(symm_pose)
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
		self.makesym.apply(symm_pose)
		return symm_pose

	def scilter(self,symm_pose,min_score,er_cutoff,min_score_pdb,current_pdb,last_pdb,scored_file_name):
		"""self.scilter scores the symmetric pose, compares the score to the min_score*er_cutoff, and determines whether to update the min_score_pdb or not 
		with the current PDB. It updates the scores_file with the final min_score_pdb of any given setting. Finally, it removes unneeded files.
		It returns a tuple with the min_score,min_score_pdb,and a boolean reporting whether the current PDB exceeded the er_cutoff (True if exceeded).
		"""
		score = self.sf(symm_pose)
		print(f'minimum score and PDB: {min_score}, {min_score_pdb}, current score and PDB: {score}, {current_pdb}')
		if min_score_pdb!=None:
			with open(os.path.join(self.base,f'scores_file.txt'),'a') as scores:
				scores.write(f'{scored_file_name}: {current_pdb}\n')
				scores.write(f'Score: {score}\n')
				specs = str(current_pdb).split("_")
				for spec in range(len(specs)):
					if str(os.path.join(self.base,f'{self.TELSAM_version}--{self.client_pdb}')) in specs[spec]:
						if len(specs)==spec+4:
							self.energies_vs_ucab_vs_deg['linker'].append(int(specs[spec+1]))
							self.energies_vs_ucab_vs_deg["ucab"].append(float(specs[spec+2]))
							self.energies_vs_ucab_vs_deg["deg"].append(float(specs[spec+3].removesuffix(".pdb")))
							self.energies_vs_ucab_vs_deg['energy'].append(float(score))
		if 0<score<=min_score*er_cutoff or score<0:
			if score<min_score:
				#Add these lines back in if you no longer want to see all the previous minimum-scoring PDBs
				#if min_score_pdb!=None:
				#	os.remove(min_score_pdb)
				min_score = score
				min_score_pdb = current_pdb
				if self.centroids:
					self.to_fullatom.apply(symm_pose)
				self.pmm.apply(symm_pose)
			if last_pdb!=None:
				if os.path.exists(last_pdb) and last_pdb!= min_score_pdb:
					os.remove(last_pdb)
					print(f'removed previous pdb: {last_pdb}')
			return (min_score,min_score_pdb,False)
		else:
			if last_pdb!=None:
				if os.path.exists(last_pdb) and last_pdb!= min_score_pdb:
					os.remove(last_pdb)
			#if os.path.exists(current_pdb) and current_pdb!=min_score_pdb:
			#	os.remove(current_pdb)#This may cause issues...
			return (min_score,min_score_pdb,True)
		
	def chart(self,linker):
		"""chart primarily makes a heat map of energies related to unit cell ab and degree of rotation for each linker length variant.
		chart can also be used to make a 3d graph of the relationship between energy, unit cell ab, and degree for each linker length variant.
		"""
		fig = plt.figure()
		data = {
			"aboi":[ucab for ucab, l in zip(self.energies_vs_ucab_vs_deg['ucab'],self.energies_vs_ucab_vs_deg['linker']) if l == linker],
			"doi":[deg for deg, l in zip(self.energies_vs_ucab_vs_deg['deg'],self.energies_vs_ucab_vs_deg['linker']) if l == linker],
			"eoi":[energy for energy, l in zip(self.energies_vs_ucab_vs_deg['energy'],self.energies_vs_ucab_vs_deg['linker']) if l == linker]
		}
		lookup = {
			(a,d): e
			for a,d,e in zip(data['aboi'],data['doi'],data['eoi'])
		}
		array_data = []
		aboi_set_list = list(set(data['aboi']))
		aboi_set_list.sort()
		doi_set_list = list(set(data['doi']))
		doi_set_list.sort()
		for degree_rotation in doi_set_list:
			row = []
			for unit_cell_ab in aboi_set_list:
				row.append(lookup.get((unit_cell_ab,degree_rotation),0))
			array_data.append(row)
		energy_array = np.array(array_data)
		n_rows = len(doi_set_list)
		n_cols = len(aboi_set_list)
		cell_size=3
		fig,ax = plt.subplots(figsize=(n_cols * cell_size,n_rows*cell_size))
		im = ax.imshow(energy_array)
		ax.set_xticks(range(len(aboi_set_list)),labels=aboi_set_list)
		ax.set_yticks(range(len(doi_set_list)),labels=doi_set_list)
		for i in range(len(doi_set_list)):
			for j in range(len(aboi_set_list)):
				text = ax.text(j,i,"{:.3e}".format(energy_array[i,j]),
				   ha='center',va='center',color='w')
		ax.set_title(f'Energies of AB:Degree Combinations for {self.TELSAM_version}--{self.client_pdb}_{linker}')
		fig.savefig(os.path.join(self.base,f"Energies of AB_Degree Combinations for {self.TELSAM_version}--{self.client_pdb}_{linker}"))
		with open(os.path.join(self.base,f"{linker}_chart.json"),"w") as file:
			json.dump(data,file,indent=4)

		"""
		print(data["aboi"][::],data["doi"][::],data["eoi"][::])
		ax = fig.add_subplot(projection='3d')
		ax.scatter(data["aboi"],data["doi"],data["eoi"])
		ax.set_title(f'Energies of AB:Degree Combinations for {self.TELSAM_version}--{self.client_pdb}_{linker}')
		ax.set_xlabel('Unit Cell AB Length (Angstroms)')
		ax.set_ylabel('Degree of Polymer Rotation (Degrees)')
		ax.set_zlabel('Energy (REU)')
		fig.savefig(os.path.join(self.base,f"Energies of AB_Degree Combinations for {self.TELSAM_version}--{self.client_pdb}_{linker}"))
		plt.close(fig)
		"""

	def validate_TELSAM(self):
		if self.remake_TELSAM_bool:
			self.remake_TELSAM()
		#Get pre-made .pdb:
		if os.path.exists(os.path.join(self.base,f'TELSAM_in_9DOC.pdb')):
			try:
				self.TELSAM_in_9DOC = pose_from_file(os.path.join(self.base,f'TELSAM_in_9DOC.pdb'))
			except Exception:
				self.remake_TELSAM()
				self.TELSAM_in_9DOC = pose_from_file(os.path.join(self.base,f'TELSAM_in_9DOC.pdb'))
		else:
			self.remake_TELSAM()
			self.TELSAM_in_9DOC = pose_from_file(os.path.join(self.base,f'TELSAM_in_9DOC.pdb'))

	def fuse(self):
		########################## CREATE TELSAM FUSION! ######################################
		try:
			self.client = Pose()
			if ".pdb" not in self.client_pdb:
				url = f"https://files.rcsb.org/download/{self.client_pdb}.pdb"
				pdb_text = requests.get(url).text
				with open(os.path.join(self.base,f"{self.client_pdb}.pdb"),"w") as file:
					file.write(pdb_text)
			cleanATOM(os.path.join(self.base,f"{self.client_pdb}.pdb"))
			temp_pose = Pose()
			pose_from_pdb(temp_pose,os.path.join(self.base,f"{self.client_pdb}.clean.pdb"))
			os.remove(os.path.join(self.base,f"{self.client_pdb}.pdb"))
			os.remove(os.path.join(self.base,f"{self.client_pdb}.clean.pdb"))
			
			#Extract first 4-aa helical region from target protein to fuse to TELSAM:
			dssp = Dssp(temp_pose)
			dssp.insert_ss_into_pose(temp_pose)
			ss_string = temp_pose.secstruct()
			first_helix = ss_string.find("HHHHH")
			append_subpose_to_pose(self.client,temp_pose,temp_pose.chain_begin(1)+first_helix,temp_pose.chain_end(1))
			
			#Align the two helices:
			self.start_residue_to_superimpose = 2
			if os.path.exists(os.path.join(self.base,f'scores_file.txt')):
				os.remove(os.path.join(self.base,f'scores_file.txt'))
			self.TELSAM = self.TELSAM_in_9DOC.clone()
			if self.linker_variant!=None:
				self.start_residue_to_superimpose+=self.linker_variant
			else:
				self.linker_variant = -1+self.start_residue_to_superimpose
			self.start_residue_to_superimpose = self.start_residue_to_superimpose+self.linker_variant
			TELSAM_residues_to_superimpose = range(self.TELSAM.chain_end(1)-self.start_residue_to_superimpose,self.TELSAM.chain_end(1)-self.start_residue_to_superimpose+3)
			client_residues_to_superimpose = range(1,4)
			atom_map = AtomID_Map()
			initialize_atomid_map(atom_map, self.client, AtomID())

			#Map CA atoms between residues
			for CR, TR in zip(client_residues_to_superimpose,TELSAM_residues_to_superimpose):
				client_atom = AtomID(self.client.residue(CR).atom_index("CA"), CR)
				TELSAM_atom = AtomID(self.TELSAM.residue(TR).atom_index("CA"), TR)
				atom_map.set(client_atom,TELSAM_atom)
			superimpose_pose(self.client,self.TELSAM,atom_map)

			#Delete overlap
			delete_region(self.TELSAM,self.TELSAM.chain_end(1)-self.start_residue_to_superimpose,self.TELSAM.chain_end(1))
			#Fuse
			append_pose_to_pose(self.TELSAM,self.client,new_chain=False)
			self.TELSAM.conformation().declare_chemical_bond(self.TELSAM.chain_end(1)-self.client.total_residue(),"C",self.TELSAM.chain_end(1)-self.client.total_residue()+1,"N")

			#Save so you can extract the furthest_x coordinate and refine correctly
			self.current_linker_pdb = os.path.join(self.base,f'{self.TELSAM_version}--{self.client_pdb}_{self.linker_variant}.pdb')
			self.TELSAM.dump_pdb(self.current_linker_pdb)
			with open(self.current_linker_pdb, 'r') as file:
				lines = iter(file)
				for line in lines:
					if 'CRYST1' in line:
						p = re.compile(r'\d+\.\d+')
						a = float(p.search(line).group())

					if 'ATOM' in line:
						x_coord = float(line[31:39].strip())
						if x_coord>self.furthest_x:
							self.furthest_x = x_coord
			#Refine
			self.refinement(self.TELSAM)
			self.TELSAM.dump_pdb(self.current_linker_pdb)
			
			#Filter:
			if self.centroids:
				self.to_centroid.apply(self.TELSAM)
			score = self.sf(self.TELSAM)
			if score<10000:
				#Record self.base score of linker variant:
				with open(os.path.join(self.base,f'scores_file.txt'),'a') as scores:
					scores.write(f'Linker file: {self.current_linker_pdb}\n')
					scores.write(f'score: {score}\n')
					#Determine largest unit cell:
			
		except Exception as e:
			print(e,file=sys.stderr)
			sys.exit()

	def refinement(self,pose):
		#################### SETUP FOR REFINEMENT ##########################
		self.interaction_shell_size = self.furthest_x*0.65
		print(f'INTERACTION SHELL SIZE: {self.interaction_shell_size}')
		rosetta.basic.options.set_real_option("cryst:interaction_shell", self.interaction_shell_size)
		self.makesym = SetupForSymmetryMover("CRYST1")

		self.sf = None
		#Get score_function
		if self.centroids:
			self.sf = create_score_function("cen_std")
		else:
			self.sf = get_score_function()

		#Relax mover
		self.relax = FastRelax()
		self.relax.set_scorefxn(self.sf)

		self.tf = TaskFactory()
		self.tf.push_back(InitializeFromCommandline())
		self.tf.push_back(RestrictToRepacking())

		if not self.centroids:
			self.packer = PackRotamersMover(self.sf)
			self.packer.task_factory(self.tf)

			self.min_mover = MinMover()
			self.min_mover.score_function(self.sf)
			self.min_mover.min_type("lbfgs_armijo_nonmonotone")

		####################### MOVEMAP REFINEMENT (Must be re-setup after pose is symmetrized) #################
		#It may be okay to just have it here and then call the min_mover.movemap and relax.set_movemap functions later.
		movemap = MoveMap()
		movemap.set_bb(False)
		movemap.set_chi(False)
		for i in range(self.TELSAM.chain_end(1)-self.client.total_residue()-self.start_residue_to_superimpose,self.TELSAM.chain_end(1)):
			movemap.set_bb(i, True)
			movemap.set_chi(i, True)
		#Refine gently
		if not self.centroids:
			self.min_mover.movemap(movemap)
			#self.packer.apply(pose)
			self.min_mover.apply(pose)
			self.relax.set_movemap(movemap)
			self.relax.apply(pose)
		##############################################################################################

	def stepper(self):
		######################################## DOCK POLYMERS! ######################################
		#Test different unit cell sizes:
		min_score = 200000
		min_score_pdb = None
		min_ucab_pdb = None
		ucab_start = None
		ucab_end = None
		ucab1 = None
		if self.optimize:
			if self.unit_cell_ab!=None:
				ucab_start = int(self.unit_cell_ab)
				ucab_end = int(self.unit_cell_ab-30)
			else:
				ucab_start = int(self.furthest_x*2)
				ucab_end = int(self.furthest_x-30)
			if self.degree_rotation!=None:
				deg_start = self.degree_rotation
				deg_end = self.degree_rotation+20
			else:
				deg_start = 1
				deg_end = 21
		else:
			ucab_start = int(self.furthest_x*2)
			ucab_end = int(self.furthest_x-30)
			deg_start = 1
			deg_end = 21

		for ucab in range(ucab_start,ucab_end,-1):
			current_ucab_pdb = os.path.join(self.base,f'{self.TELSAM_version}--{self.client_pdb}_{self.linker_variant}_{ucab}.pdb')
			self.change_cell(self.current_linker_pdb,current_ucab_pdb,wa=ucab,wb=ucab)
			symm_pose = pose_from_file(current_ucab_pdb)
			if self.centroids:
				self.to_centroid.apply(symm_pose)
			self.makesym.apply(symm_pose)
			min_score,min_score_pdb,exceeded = self.scilter(symm_pose,min_score,1.1,min_score_pdb,current_ucab_pdb,os.path.join(self.base,f'{self.TELSAM_version}--{self.client_pdb}_{self.linker_variant}_{ucab+1}.pdb'),'Unit Cell AB File')
			if exceeded:
				if min_score_pdb != None:
					min_ucab_pdb = min_score_pdb
					ucab1 = int(str(min_ucab_pdb).split("_")[-1].removesuffix(".pdb"))
				break
		##################################### ROTATE UNIT CELL AND CONTINUE DOCKING! ###############################
		if min_ucab_pdb!=None:
			ucab = ucab1
			min_score = 200000
			min_score_pdb = None
			for deg in range(deg_start-1,deg_end):
				current_deg_pdb = os.path.join(self.base,f'{min_ucab_pdb.removesuffix(".pdb")}_{deg}.pdb')
				for ucab2 in range(ucab,ucab_end,-1):
					current_ucab2_deg_pdb = os.path.join(self.base,f'{min_ucab_pdb.removesuffix(f"{ucab1}.pdb")}{ucab2}_{deg}.pdb')
					self.change_cell(min_ucab_pdb,current_ucab2_deg_pdb,wa=ucab2,wb=ucab2)
					symm_pose = self.change_deg(current_ucab2_deg_pdb,deg)
					self.makesym.apply(symm_pose)
					min_score,min_score_pdb,exceeded = self.scilter(symm_pose,min_score,1.1,min_score_pdb,current_ucab2_deg_pdb,os.path.join(self.base,f'{min_ucab_pdb.removesuffix(f"{ucab1}.pdb")}_{ucab2+1}_{deg}.pdb'),'Temp Unit Cell AB2 Deg File')
					if self.centroids:
						self.to_fullatom.apply(symm_pose)
					self.pmm.apply(symm_pose)
					if exceeded:
						ucab = ucab2
						break
				if self.centroids:
					self.to_centroid.apply(symm_pose)
				min_score,min_score_pdb,exceeded = self.scilter(symm_pose,min_score,3,min_score_pdb,current_ucab2_deg_pdb,None,'Unit Cell AB2 Deg File')
				#EITHER leave this alone (meaning you will always test all the degrees in the range) OR calculate the angle between the linker vector and
				#the furthest edge of the linker, because that angle has to be the starter angle for whichever degree check is opposite
				#that angle (since proceeding by just 1 degree in the opposing direction will cause the edge of the linker to continually bump and increase
				#in energy)

			ucab = ucab1
			min_score = 200000
			min_score_pdb = None
			for deg in range(-deg_start,-deg_end,-1):
				current_deg_pdb = os.path.join(self.base,f'{min_ucab_pdb.removesuffix(".pdb")}_{deg}.pdb')
				for ucab2 in range(ucab,ucab_end,-1):
					current_ucab2_deg_pdb = os.path.join(self.base,f'{min_ucab_pdb.removesuffix(f"{ucab1}.pdb")}{ucab2}_{deg}.pdb')
					self.change_cell(min_ucab_pdb,current_ucab2_deg_pdb,wa=ucab2,wb=ucab2)
					symm_pose = self.change_deg(current_ucab2_deg_pdb,deg)
					self.makesym.apply(symm_pose)
					min_score,min_score_pdb,exceeded = self.scilter(symm_pose,min_score,1.1,min_score_pdb,current_ucab2_deg_pdb,os.path.join(self.base,f'{min_ucab_pdb.removesuffix(f"{ucab1}.pdb")}_{ucab2+1}_{deg}.pdb'),'Temp Unit Cell AB2 Deg File')
					if self.centroids:
						self.to_fullatom.apply(symm_pose)
					self.pmm.apply(symm_pose)
					if exceeded:
						ucab = ucab2 #This assumes that client proteins are not too concave (pointing in toward the polymer vector)
						break
				if self.centroids:
					self.to_centroid.apply(symm_pose)
				min_score,min_score_pdb,exceeded = self.scilter(symm_pose,min_score,3,min_score_pdb,current_ucab2_deg_pdb,None,'Unit Cell AB2 Deg File')
			
			self.chart(self.linker_variant)

		"""
		##################################### CHANGE HELICAL RISE! #####################################
		#Test different unit cell sizes:
		min_score = 200000
		for ucc in range(1,15):
			self.change_cell(passing_pdb,os.path.join(self.base,f'{passing_pdb.removesuffix(".pdb")}_{ucc}.pdb'),wc=ucc)
			symm_pose = pose_from_file(os.path.join(self.base,f'{passing_pdb.removesuffix(".pdb")}_{ucc}.pdb'))
			self.makesym.apply(symm_pose)
			min_score,passing_pdb,finished = self.scilter(symm_pose,min_score,passing_pdb,min_ucab_pdb,os.path.join(self.base,f'{self.TELSAM_version}--{self.client_pdb}_{self.linker_variant}_{ucab+1}.pdb'),'Unit Cell C File')
			if finished:
				break
		"""

	def picker(self):
		current_ucab_pdb = os.path.join(self.base,f'{self.TELSAM_version}--{self.client_pdb}_{self.linker_variant}_{str(int(float(self.unit_cell_ab)))}.pdb')
		self.change_cell(self.current_linker_pdb,current_ucab_pdb,wa=self.unit_cell_ab,wb=self.unit_cell_ab)
		current_deg_pdb = os.path.join(self.base,f'{current_ucab_pdb.removesuffix(".pdb")}_{self.degree_rotation}.pdb')
		symm_pose = self.change_deg(current_ucab_pdb,self.degree_rotation)
		sequence = symm_pose.sequence()
		if self.centroids:
			self.to_centroid.apply(symm_pose)
		self.makesym.apply(symm_pose)
		score = self.sf(symm_pose)
		print(f'Score of {current_deg_pdb}: {score}')
		with open (f'{current_deg_pdb.removesuffix('.pdb')}.fasta', 'w') as f:
			f.write(">"+current_deg_pdb+", score (REU): "+"{:.3e}".format(score)+"\n"+str(sequence).strip('X'))
		symm_pose.dump_pdb(current_deg_pdb)
		if self.centroids:
			self.to_fullatom.apply(symm_pose)
		self.pmm.apply(symm_pose)

def main():
	TELSetta1 = TELSetta()
	if TELSetta1.TELSAM_version=="1TEL" and TELSetta1.client_pdb!=None:
		TELSetta1.fuse()
		if bool(TELSetta1.unit_cell_ab==None and TELSetta1.degree_rotation==None) or TELSetta1.optimize==True:
			TELSetta1.stepper()
		elif TELSetta1.unit_cell_ab!=None and TELSetta1.degree_rotation!=None and TELSetta1.optimize==False:
			TELSetta1.picker()
		else:
			print(f"Not enough arguments were provided to pick a specific TELSAM--fusion variant to model.")
	else:
		print(f"You didn't pass any client proteins to fuse to TELSAM.")

main()