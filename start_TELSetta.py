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
import pyrosetta.rosetta.core.io
dir(pyrosetta.rosetta.core.io)

from pyrosetta.rosetta.core.pose import append_subpose_to_pose
from pyrosetta.rosetta.core.pose import append_pose_to_pose
from pyrosetta.rosetta.core.pose import initialize_atomid_map
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.id import AtomID_Map_AtomID as AtomID_Map
from pyrosetta.rosetta.core.scoring import superimpose_pose
from pyrosetta.rosetta.protocols.grafting import delete_region

#Movemap Factory and Selectors for interface refinement/scoring
from pyrosetta.rosetta.core.select.movemap import MoveMapFactory, move_map_action
from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.core.select.residue_selector import NeighborhoodResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import OrResidueSelector

from pyrosetta.rosetta.core.scoring.dssp import Dssp
from pyrosetta.rosetta.core.scoring import fa_rep
from pyrosetta.rosetta.core.scoring import fa_atr
#from pyrosetta.rosetta.protocols.simple_filters import ShapeComplementarityFilter
#from pyrosetta.rosetta.core.scoring.sc import ShapeComplementarityCalculator
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

from pyrosetta import PyMOLMover
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.select.movemap import move_map_action
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.scoring import get_score_function
#from pyrosetta.rosetta.core.scoring import get_fa_score_function
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory

from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import InitializeFromCommandline, RestrictToRepacking, OperateOnResidueSubset, RestrictToRepackingRLT
from pyrosetta.rosetta.protocols.minimization_packing import MinMover

from pyrosetta.rosetta.numeric import xyzMatrix_double_t, xyzVector_double_t

from matplotlib import pyplot as plt

#pyrosetta.init("-crystal_refine -cryst::refinable_lattice -score_symm_complex")
pyrosetta.init("-crystal_refine")

class TELSetta:
	def __init__(self):
		self.TELSAM_version = "1TEL"
		self.client_pdb = None
		self.linker_variant = None
		self.unit_cell_ab = None
		self.degree_rotation = None
		self.remake_TELSAM_bool = False
		self.optimize = False
		self.centroids = False
		self.interfaced = False
		self.scores = {}
		self.headers = {
			"User-Agent": (
			"Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
			"AppleWebKit/537.36 (KHTML, like Gecko) "
			"Chrome/125.0 Safari/537.36"
			)
		}
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

		self.base = os.path.join(os.path.expanduser('~/TELSAM-Fusion-Crystallography-with-Rosetta'),str(self.linker_variant))
		os.makedirs(self.base,exist_ok=True)
		self.pmm = PyMOLMover()
		self.pmm.keep_history(True)
		self.energies_vs_ucab_vs_deg = {'linker':[],'energy':[],'ucab':[],'deg':[]}
		self.interfaced_energies_vs_ucab_vs_deg = {'linker':[],'energy':[],'ucab':[],'deg':[]}
		self.furthest_x = 0
		self.to_fullatom = SwitchResidueTypeSetMover("fa_standard")
		self.setup_for_refinement()
		self.validate_TELSAM()

	def setup_for_refinement(self):
		#################### SETUP FOR REFINEMENT ##########################
		self.sf = get_score_function()
		#self.sf = ScoreFunctionFactory.create_score_function("beta_nov16_cart")
		#Relax mover
		self.relax = FastRelax()
		self.relax.set_scorefxn(self.sf)

		tf = TaskFactory()
		tf.push_back(InitializeFromCommandline())
		tf.push_back(RestrictToRepacking())

		#Packer Mover
		self.packer = PackRotamersMover(self.sf)
		self.packer.task_factory(tf)

		self.min_mover = MinMover()
		self.min_mover.score_function(self.sf)
		self.min_mover.min_type("lbfgs_armijo_nonmonotone")

		"""
		#Symmetry stuff
		self.interaction_shell_size = self.furthest_x*0.65
		print(f'INTERACTION SHELL SIZE: {self.interaction_shell_size}')
		rosetta.basic.options.set_real_option("cryst:interaction_shell", self.interaction_shell_size)
		self.makesym = SetupForSymmetryMover("CRYST1")
		"""

		#Shape Complementarity Filter
		self.iam = InterfaceAnalyzerMover("A_B")
		self.iam.set_scorefunction(self.sf)

	def refine(self,pose) -> float:
		"""Creates a movemap, freezing bb and chi angles, then unfreezing the atoms that are outside of the 1TEL module.
		Then, changes self.min_mover to have the movemap that was just created.
		Then, applies self.min_mover to the current pose.
		Then, changes self.relax to have the previous movemap as well.
		Then, initializes a TaskFactory that reinitializes from command line and restricts to repacking???
		Then, changes self.relax by setting the previously made TaskFactory in it.
		Then, applies self.relax to the passed pose.
		Finally, returns the energy of the pose with self.sf(pose)"""
		####################### MOVEMAP REFINEMENT (Must be re-setup after pose is symmetrized) #################
		#It may be okay to just have it here and then call the min_mover.movemap and relax.set_movemap functions later.
		movemap = MoveMap()
		movemap.set_bb(False)
		movemap.set_chi(False)
		for i in range(self.TELSAM.chain_end(1)-self.client.total_residue()-self.start_residue_to_superimpose,self.TELSAM.chain_end(1)):
			movemap.set_bb(i, True)
			movemap.set_chi(i, True)
		#Refine gently
		self.min_mover.movemap(movemap)
		#self.packer.apply(pose)
		self.min_mover.apply(pose)
		self.relax.set_movemap(movemap)
		tf = TaskFactory()
		tf.push_back(InitializeFromCommandline())
		tf.push_back(RestrictToRepacking())
		self.relax.set_task_factory(tf)
		self.relax.apply(pose)
		energy = self.sf(pose)
		return energy
		##############################################################################################

	def interface_refine(self,pose) -> float:
		"""First, creates chain selectors for chains A and B.
		Then, creates NeighborhoodResidueSelectors that accept the chain selectors as arguments...
		Then, creates an OrResidueSelector that accepts both NeighborhoodResidueSelectors as arguments...
		Then, actually generates a selection by applying the previous selector on the pose (this may be broken).
		Creates a movemap_factory that disables movement of the bb and chi angles for all but the interface selector?
		Creates a TaskFactory to restrict to repacking??? Why isn't this the default?
		Changes self.relax by setting the movemap factory and task factory to it.
		Finally, relaxes the passed pose.
		Pushes to PyMOL.
		Tries to run the InterfaceAnalyzerMover("A_B") and return interface dg. On a failure, instead runs normal refinement."""
		#Create selectors to get chains
		chainA_sel = ChainSelector("A")
		chainB_sel = ChainSelector("B")
		#chainC_sel = ChainSelector("C")
		#Create selectors to find neighbors of each chain
		chainA_neighbors_sel = NeighborhoodResidueSelector(chainA_sel,5.0,False)
		chainB_neighbors_sel = NeighborhoodResidueSelector(chainB_sel,5.0,False)
		#chainC_neighbors_sel = NeighborhoodResidueSelector(chainC_sel,5.0,False)
		interfaceAB_sel = OrResidueSelector(chainA_neighbors_sel,chainB_neighbors_sel)
		#interfaceBC_sel = OrResidueSelector(chainB_neighbors_sel,chainC_neighbors_sel)
		#Actually select residues of interest
		interfaceAB = interfaceAB_sel.apply(pose)
		print("Interface-region residues:", [i+1 for i, inter in enumerate(interfaceAB) if inter])
		#interfaceBC = interfaceBC_sel.apply(pose)
		
		#Define residues that can move
		movemap_factory = MoveMapFactory()
		movemap_factory.add_bb_action(move_map_action.mm_disable, TrueResidueSelector())
		movemap_factory.add_chi_action(move_map_action.mm_disable, TrueResidueSelector())
		movemap_factory.add_bb_action(move_map_action.mm_enable, interfaceAB_sel)
		movemap_factory.add_chi_action(move_map_action.mm_enable, interfaceAB_sel)

		#relax
		#restrict residues that can move
		tf = TaskFactory()
		tf.push_back(RestrictToRepacking())
		tf.push_back(OperateOnResidueSubset(RestrictToRepackingRLT(),interfaceAB))
		self.relax.set_movemap_factory(movemap_factory)
		self.relax.set_task_factory(tf)
		self.relax.apply(pose)

		pdbinfo = pose.pdb_info()
		for c in range(1, pose.num_chains()+1):
			print(
				c,
				pose.chain_begin(c),
				pose.chain_end(c),
				pdbinfo.chain(pose.chain_begin(c))
			)
		print(f'fold_tree: {pose.fold_tree()}')
		pose.pdb_info().name("pmm")
		self.pmm.apply(pose)
		self.iam.apply(pose)
		score = self.iam.get_interface_dG()
		"""Use the following lines instead in case of three chains
		relax_interface(interfaceBC_sel)
		self.iam.set_interface("AB_C")
		self.iam.apply(pose)
		score1 = self.iam.get_interface_dG()
		self.iam.set_interface("A_BC")
		self.iam.apply(pose)
		score2 = self.iam.get_interface_dG()
		score = score1+score2
		"""
		#print(f'SASA: {self.iam.get_interface_delta_sasa()}')
		#print(f'unsat H: {self.iam.get_interface_delta_hbond_unsat()}')
		#print(f'Interface energy: {self.iam.get_crossterm_interface_energy()}')
		if not score:
			score = self.refine(pose)
		else:
			print(f'INTERFACE WORKED!')
			if self.interfaced == False:
				self.min_score = score
			self.interfaced = True
		return score

	def get_CRYST1(self,pdb):
		with open(pdb, 'r') as s:
				for line in s:
					if "CRYST1" in line:
						p = re.compile(r'\d+\.\d+')
						cryst1_vals = p.findall(line)
						a, b, c, alpha, beta, gamma = [float(x) for x in cryst1_vals[0:6]]
						print(f'PDB: {pdb}, CRYST1: {a, b, c}')
						return (a, b, c, alpha, beta, gamma)

	def add_CRYST1(self,new_pdb,old_pdb):
		with open (os.path.join(self.base,new_pdb),'r') as file:
			pdb_sans_cryst = file.read()
		with open(os.path.join(self.base,new_pdb),'w') as file:
			with open(os.path.join(self.base,old_pdb)) as s:
				for line in s:
					if "CRYST1" in line:
						file.write(line)
						break
			file.write(pdb_sans_cryst)
		
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

	#Currently, this only works in space group P65, but that's not a problem for this code.
	def rotate_file(self,file,deg):
		deg_pose = pose_from_file(file)
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
		deg_pose.apply_transform_Rx_plus_v(R, v)
		deg_pose.dump_pdb(file)
		return deg_pose
	
	def space_group_symmop_pose_from_pdb(self,pdb,symmop:str):
		"""Performs a symmetry operator on the entire pose, returning a symmop_copy clone of the original pose in the new location."""
		symmop_pose = pose_from_file(pdb)
		#Get Crystal Info:
		cryst1 = self.get_CRYST1(pdb)
		print(f'CRYST1: {cryst1}')
		dir(symmop_pose.conformation())
		dir(symmop_pose.pdb_info())
		#Do translation:
		for res_i in range(1,symmop_pose.total_residue()+1):
			res = symmop_pose.residue(res_i)
			for atom_i in range(1,res.natoms()+1):
				xyz = res.xyz(atom_i)
				abc = ((xyz[0]+xyz[1]*math.tan(math.radians(30)))/cryst1[0],xyz[1]/math.cos(math.radians(30))/cryst1[1],xyz[2]/cryst1[2])
				if symmop == "6":
					#symmop6:
					abc = (abc[1],-abc[0]+abc[1],(1/6)+abc[2])
				elif symmop == "2":
					#symmop2:
					abc = (-abc[1],abc[0]-abc[1],(2/3)+abc[2])
				elif symmop == "4+a":
					#symmop4+a:
					abc = (-abc[0]+1,-abc[1],1/2+abc[2])
				xyz = xyzVector_double_t((abc[0]-abc[1]*math.sin(math.radians(30)))*cryst1[0],abc[1]*cryst1[1]*math.cos(math.radians(30)),abc[2]*cryst1[2])
				symmop_pose.set_xyz(AtomID(atom_i,res_i),xyz)
		return symmop_pose

	def generate_minimal_contact_symmetry_mates(self,pdb,symmops):
		"""Hopefully this function will create rotated and translated copies of a provided monomer such that all lattice contacts are being represented 
		once. This will allow downstream scoring functions to determine the lattice strength."""
		#In the instance of the P65 space group, I'm interested in just two of the symmetry mates: the monomer directly above the root (n+1, or in 
		#symmetry lingo, translate the unit cell across one a axis, then perform the sixth symmetry operation, (y,-x+y,1/6+x)) and the 
		#monomer that's across from the root, in the next polymer over (in symmetry lingo, the second symmetry operation from the root chain).
		symmop_poses = []
		pose = pose_from_file(pdb)
		for symmop in symmops:
			symmop_pose = self.space_group_symmop_pose_from_pdb(pdb,symmop)
			symmop_poses.append(symmop_pose)
		for symmop_pose in symmop_poses:
			append_pose_to_pose(pose,symmop_pose,new_chain=True)
		pdbinfo = pose.pdb_info()
		chain_letters = ["A", "B", "C"]
		for chain_num in range(1, pose.num_chains()+1):
			letter = chain_letters[chain_num-1]
			for res in range(
				pose.chain_begin(chain_num),
				pose.chain_end(chain_num)+1
			):
				pdbinfo.chain(res, letter)
		return pose

	def fill_unit_cell(pose,a,b,c):
		if pose=="P65":
			pass

	def refscilter(self,symm_pose,er_cutoff,current_pdb,last_pdb,scored_file_name) -> bool:
		"""self.refscilter refines and scores the symmetric pose, compares the score to the min_score*er_cutoff, and determines whether to update the min_score_pdb or not 
		with the current PDB. It updates the scores_file with the final min_score_pdb of any given setting. Finally, it removes unneeded files.
		It returns a tuple with the min_score,min_score_pdb,and a boolean reporting whether the current PDB exceeded the er_cutoff (True if exceeded).
		"""
		score = self.interface_refine(symm_pose)
		"""
		sequence = symm_pose.sequence()
		with open (f'{current_pdb.removesuffix('.pdb')}.fasta', 'w') as f:
			f.write(">"+current_pdb+", score (REU): "+"{:.3e}".format(score)+"\n"+"HHHHHHHHHH"+str(sequence).strip('X'))
		self.add_CRYST1(os.path.basename(current_pdb),os.path.basename(current_pdb))
		"""
		print(f'minimum score and PDB: {self.min_score}, {self.min_score_pdb}, current score and PDB: {score}, {current_pdb}')
		if self.interfaced:
			scores_file = os.path.join(self.base,f'interfaced_scores_file.txt')
			if self.min_score_pdb!=None:
				with open(scores_file,'a') as scores:
					scores.write(f'{scored_file_name}: {current_pdb}\n')
					scores.write(f'Score: {score}\n')
					specs = str(current_pdb).split("_")
					for spec in range(len(specs)):
						if str(os.path.join(self.base,f'{self.TELSAM_version}--{self.client_pdb}')) in specs[spec]:
							if len(specs)==spec+4:
								self.interfaced_energies_vs_ucab_vs_deg['linker'].append(int(specs[spec+1]))
								self.interfaced_energies_vs_ucab_vs_deg["ucab"].append(int(specs[spec+2]))
								self.interfaced_energies_vs_ucab_vs_deg["deg"].append(int(specs[spec+3].removesuffix(".pdb")))
								self.interfaced_energies_vs_ucab_vs_deg['energy'].append(float(score))
			symm_pose.pdb_info().name("pmm")
			self.pmm.apply(symm_pose)
			if score<0 or 0<score<=abs(self.min_score)*er_cutoff:
				if score<self.min_score:
					#Add these lines back in if you no longer want to see all the previous minimum-scoring PDBs
					#if self.min_score_pdb!=None:
					#	os.remove(self.min_score_pdb)
					self.min_score = score
					self.min_score_pdb = current_pdb
				if last_pdb!=None:
					if os.path.exists(last_pdb) and last_pdb!= self.min_score_pdb:
						os.remove(last_pdb)
						print(f'removed previous pdb: {last_pdb}')
				print(f'New minimum score and PDB: {self.min_score}, {self.min_score_pdb}')
				return (False)
			else:
				if last_pdb!=None:
					if os.path.exists(last_pdb) and last_pdb!= self.min_score_pdb:
						os.remove(last_pdb)
				#if os.path.exists(current_pdb) and current_pdb!=self.min_score_pdb:
				#	os.remove(current_pdb)#This may cause issues...
				print(f'New minimum score and PDB: {self.min_score}, {self.min_score_pdb}')
				return (True)
		else:
			scores_file = os.path.join(self.base,f'scores_file.txt')
			if self.min_score_pdb!=None:
				with open(scores_file,'a') as scores:
					scores.write(f'{scored_file_name}: {current_pdb}\n')
					scores.write(f'Score: {score}\n')
					specs = str(current_pdb).split("_")
					for spec in range(len(specs)):
						if str(os.path.join(self.base,f'{self.TELSAM_version}--{self.client_pdb}')) in specs[spec]:
							if len(specs)==spec+4:
								self.energies_vs_ucab_vs_deg['linker'].append(int(specs[spec+1]))
								self.energies_vs_ucab_vs_deg["ucab"].append(int(specs[spec+2]))
								self.energies_vs_ucab_vs_deg["deg"].append(int(specs[spec+3].removesuffix(".pdb")))
								self.energies_vs_ucab_vs_deg['energy'].append(float(score))
			symm_pose.pdb_info().name("pmm")
			self.pmm.apply(symm_pose)
			if score<0 or 0<score<=abs(self.min_score)*er_cutoff:
				if score<self.min_score:
					#Add these lines back in if you no longer want to see all the previous minimum-scoring PDBs
					#if self.min_score_pdb!=None:
					#	os.remove(self.min_score_pdb)
					self.min_score = score
					self.min_score_pdb = current_pdb
				if last_pdb!=None:
					if os.path.exists(last_pdb) and last_pdb!= self.min_score_pdb:
						os.remove(last_pdb)
						print(f'removed previous pdb: {last_pdb}')
				print(f'New minimum score and PDB: {self.min_score}, {self.min_score_pdb}')
				return (False)
			else:
				if last_pdb!=None:
					if os.path.exists(last_pdb) and last_pdb!= self.min_score_pdb:
						os.remove(last_pdb)
				#if os.path.exists(current_pdb) and current_pdb!=self.min_score_pdb:
				#	os.remove(current_pdb)#This may cause issues...
				print(f'New minimum score and PDB: {self.min_score}, {self.min_score_pdb}')
				return (True)
		
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

	def remake_TELSAM(self):
		if os.path.exists(os.path.join(self.base,f'TELSAM_in_9DOC.pdb')):
			os.remove(os.path.join(self.base,f'TELSAM_in_9DOC.pdb'))
		#Get 2QAR from the pdb (it's our most convenient engineered TELSAM because of the long helical linker at the end.)
		url = f"https://files.rcsb.org/download/2QAR.pdb"
		pdb_text = requests.get(url,headers=self.headers).text
		with open(os.path.join(self.base,"ETEL.pdb"),"w") as file:
			file.write(pdb_text)
		cleanATOM(os.path.join(self.base,"ETEL.pdb"))
		os.remove(os.path.join(self.base,"ETEL.pdb"))

		#Get 9DOC from the pdb (it has the proper space group that TELSAM usually fits into)
		url = f"https://files.rcsb.org/download/9DOC.pdb"
		pdb_text = requests.get(url,headers=self.headers).text
		with open(os.path.join(self.base,"STEL.pdb"),"w") as file:
			file.write(pdb_text)
		cleanATOM(os.path.join(self.base,"STEL.pdb"))
		#Don't remove the STEL.pdb yet because it has crystallographic information we want to get later.

		################### MOVE 2QAR INTO 9DOC ASYMMETRIC UNIT ##############################
		#Grab a portion of 2QAR
		TELSAM_in_9DOC = Pose()
		temp_pose = pose_from_pdb(os.path.join(self.base,'ETEL.clean.pdb'))
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
		mutate_residue(TELSAM_in_9DOC,34,"R",5)
		mutate_residue(TELSAM_in_9DOC,66,"E",5)
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
		self.add_CRYST1(f'TELSAM_in_9DOC.pdb',f'STEL.pdb')
		if os.path.exists(os.path.join(self.base,"STEL.pdb")):
			os.remove(os.path.join(self.base,"STEL.pdb"))
		print(f'Remade TELSAM. Stored in: {self.base}')

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
				pdb_text = requests.get(url,headers=self.headers).text
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
			self.start_residue_to_superimpose = 3
			if os.path.exists(os.path.join(self.base,f'scores_file.txt')):
				os.remove(os.path.join(self.base,f'scores_file.txt'))
			if os.path.exists(os.path.join(self.base,f'interfaced_scores_file.txt')):
				os.remove(os.path.join(self.base,f'interfaced_scores_file.txt'))
			self.TELSAM = self.TELSAM_in_9DOC.clone()
			if self.linker_variant!=None:
				self.start_residue_to_superimpose+=self.linker_variant
			else:
				self.linker_variant = 1+self.start_residue_to_superimpose
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

			#Refine
			movemap = MoveMap()
			movemap.set_bb(False)
			movemap.set_chi(False)
			for i in range(self.TELSAM.chain_end(1)-self.client.total_residue()-self.start_residue_to_superimpose,self.TELSAM.chain_end(1)):
				movemap.set_bb(i, True)
				movemap.set_chi(i, True)
			#Refine gently
			self.min_mover.movemap(movemap)
			#self.packer.apply(self.TELSAM)
			self.min_mover.apply(self.TELSAM)
			self.relax.set_movemap(movemap)
			self.relax.apply(self.TELSAM)

			#Realign to 9DOC at the polymer extension interface
			TELSAM_residues_to_superimpose = [2,30,31,32,34,53,57,62,65,69]
			atom_map = AtomID_Map()
			initialize_atomid_map(atom_map, self.TELSAM_in_9DOC, AtomID())
			for R in (TELSAM_residues_to_superimpose):
				E_atom = AtomID(self.TELSAM.residue(R).atom_index("CA"), R)
				S_atom = AtomID(self.TELSAM_in_9DOC.residue(R).atom_index("CA"), R)
				atom_map.set(E_atom,S_atom)
			superimpose_pose(self.TELSAM_in_9DOC,self.TELSAM,atom_map)

			#Save so you can extract the furthest_x coordinate and refine correctly
			self.current_linker_pdb = os.path.join(self.base,f'{self.TELSAM_version}--{self.client_pdb}_{self.linker_variant}.pdb')
			self.TELSAM.dump_pdb(self.current_linker_pdb)
			self.add_CRYST1(os.path.basename(self.current_linker_pdb),os.path.basename("TELSAM_in_9DOC.pdb"))
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

			#Filter:
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

	def stepper(self):
		######################################## DOCK POLYMERS! ######################################
		#Test different unit cell sizes:
		self.min_score = 200000
		self.min_score_pdb = None
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
				ucab_end = int(self.furthest_x*2-30)
			if self.degree_rotation!=None:
				deg_start = self.degree_rotation
				deg_end = self.degree_rotation+20
			else:
				deg_start = 1
				deg_end = 21
		else:
			ucab_start = int(self.furthest_x*2)
			ucab_end = int(self.furthest_x*2-30)
			deg_start = 1
			deg_end = 21

		for ucab in range(ucab_start,ucab_end,-1):
			current_ucab_pdb = os.path.join(self.base,f'{self.TELSAM_version}--{self.client_pdb}_{self.linker_variant}_{ucab}.pdb')
			self.change_cell(self.current_linker_pdb,current_ucab_pdb,wa=ucab,wb=ucab)
			symm_pose = self.generate_minimal_contact_symmetry_mates(current_ucab_pdb,("4+a",))
			exceeded = self.refscilter(symm_pose,1.1,current_ucab_pdb,os.path.join(self.base,f'{self.TELSAM_version}--{self.client_pdb}_{self.linker_variant}_{ucab+1}.pdb'),'Unit Cell AB File')
			if exceeded:
				if self.min_score_pdb != None:
					min_ucab_pdb = self.min_score_pdb
					ucab1 = int(str(min_ucab_pdb).split("_")[-1].removesuffix(".pdb"))
				break
		min_ucab_pdb = self.min_score_pdb
		ucab1 = int(str(min_ucab_pdb).split("_")[-1].removesuffix(".pdb"))
		##################################### ROTATE UNIT CELL AND CONTINUE DOCKING! ###############################
		if min_ucab_pdb!=None:
			print(f'min_ucab_pdb: {min_ucab_pdb}')
			ucab = ucab1
			self.interfaced = False
			self.min_score = 200000
			self.min_score_pdb = None
			for deg in range(deg_start-1,deg_end):
				for ucab2 in range(ucab,ucab_end,-1):
					current_ucab2_deg_pdb = os.path.join(self.base,f'{min_ucab_pdb.removesuffix(f"{ucab1}.pdb")}{ucab2}_{deg}.pdb')
					self.change_cell(min_ucab_pdb,current_ucab2_deg_pdb,wa=ucab2,wb=ucab2)
					self.rotate_file(current_ucab2_deg_pdb,deg)
					symm_pose = self.generate_minimal_contact_symmetry_mates(current_ucab2_deg_pdb,("4+a",))
					exceeded = self.refscilter(symm_pose,1.1,current_ucab2_deg_pdb,os.path.join(self.base,f'{min_ucab_pdb.removesuffix(f"{ucab1}.pdb")}_{ucab2+1}_{deg}.pdb'),'Temp Unit Cell AB2 Deg File')
					if exceeded:
						ucab = ucab2
						break
				exceeded = self.refscilter(symm_pose,3,current_ucab2_deg_pdb,None,'Unit Cell AB2 Deg File')
				#EITHER leave this alone (meaning you will always test all the degrees in the range) OR calculate the angle between the linker vector and
				#the furthest edge of the linker, because that angle has to be the starter angle for whichever degree check is opposite
				#that angle (since proceeding by just 1 degree in the opposing direction will cause the edge of the linker to continually bump and increase
				#in energy)

			ucab = ucab1
			self.interfaced = False
			self.min_score = 200000
			self.min_score_pdb = None
			for deg in range(-deg_start,-deg_end,-1):
				for ucab2 in range(ucab,ucab_end,-1):
					current_ucab2_deg_pdb = os.path.join(self.base,f'{min_ucab_pdb.removesuffix(f"{ucab1}.pdb")}{ucab2}_{deg}.pdb')
					self.change_cell(min_ucab_pdb,current_ucab2_deg_pdb,wa=ucab2,wb=ucab2)
					self.rotate_file(current_ucab2_deg_pdb,deg)
					symm_pose = self.generate_minimal_contact_symmetry_mates(current_ucab2_deg_pdb,("4+a",))
					exceeded = self.refscilter(symm_pose,1.1,current_ucab2_deg_pdb,os.path.join(self.base,f'{min_ucab_pdb.removesuffix(f"{ucab1}.pdb")}_{ucab2+1}_{deg}.pdb'),'Temp Unit Cell AB2 Deg File')
					if exceeded:
						ucab = ucab2 #This assumes that client proteins are not too concave (pointing in toward the polymer vector)
						break
				exceeded = self.refscilter(symm_pose,3,current_ucab2_deg_pdb,None,'Unit Cell AB2 Deg File')
			
			self.chart(self.linker_variant)

	def picker(self):
		current_ucab_pdb = os.path.join(self.base,f'{self.TELSAM_version}--{self.client_pdb}_{self.linker_variant}_{str(int(float(self.unit_cell_ab)))}.pdb')
		self.change_cell(self.current_linker_pdb,current_ucab_pdb,wa=self.unit_cell_ab,wb=self.unit_cell_ab)
		self.rotate_file(current_ucab_pdb,self.degree_rotation)
		symm_pose = self.generate_minimal_contact_symmetry_mates(current_ucab_pdb,("4+a",))
		sequence = symm_pose.sequence()
		self.interface_refine(symm_pose)
		score = self.sf(symm_pose)
		print(f'Score of {current_ucab_pdb}: {score}')
		with open (f'{current_ucab_pdb.removesuffix('.pdb')}.fasta', 'w') as f:
			f.write(">"+current_ucab_pdb+", score (REU): "+"{:.3e}".format(score)+"\n"+"HHHHHHHHHH"+str(sequence).strip('X'))
		self.add_CRYST1(os.path.basename(current_ucab_pdb),os.path.basename(current_ucab_pdb))
		symm_pose.pdb_info().name("pmm")
		self.pmm.apply(symm_pose)
		self.chart(self.linker_variant)

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