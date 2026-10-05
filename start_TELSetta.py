import sys
import subprocess
import requests
import os
import getopt
import re
import time
import math
import random
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
from pyrosetta.rosetta.core.select.residue_selector import AndResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import NotResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.core.select.residue_selector import VirtualResidueSelector

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
from pyrosetta.rosetta.core.pose.symmetry import is_symmetric

from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.scoring import get_score_function
#from pyrosetta.rosetta.core.scoring import get_fa_score_function
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory

from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import InitializeFromCommandline, RestrictToRepacking, OperateOnResidueSubset, RestrictToRepackingRLT
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.pack.task.operation import PreventRepackingRLT

from pyrosetta.rosetta.numeric import xyzMatrix_double_t, xyzVector_double_t

from matplotlib import pyplot as plt

pyrosetta.init("-crystal_refine -cryst::refinable_lattice -score_symm_complex -mute all")
rosetta.basic.options.set_real_option("cryst:interaction_shell",12.0)
#pyrosetta.init()

class TELSetta:
	"""Using PyRosetta, creates a 1TEL fusion to a target protein to create a TFC (TELSAM Fusion Construct), docks the TFC against itself in the P65 space group at several\
 rotational offsets about the C-axis of the unit cell and distance offsets with respect to the A and B lengths of the unit cell, scores each docking\
 attempt, and returns a list of scores with their associated settings.\n
	Accepts as arguments (-flag "arg"):\n
	-t "TELSAM_version" (not yet implemented; only accepts "1TEL" or nothing)\n
	-c "client_pdb" (PDB id to fetch or path to a PDB file)\n
	-l "linker_variant" (number of amino acids to place in the rigid helical fusion between 1TEL and the client)\n
	-u "unit_cell_ab" (unit cell size in which to place the TFCs to be docked. Optional; mainly used if you're trying to reproduce a precise docking.)\n
	-d "degree_rotation" (the degree of rotation offset by which to turn the TFCs about the C axis in the crystal's default unit cell)\n
	-r "remake_TELSAM" (bool indicating whether the 1TEL subunit should be remade or reused from a previous session.
	Will be deprecated in favor of recreating the 1TEL subunit every time.)\n
	-o "optimize" (bool indicating whether a precise set of inputs including linker_variant and unit_cell_ab should be further refined.)\n
	-e "exhaustive" (bool indicating whether all poses should be modeled and scored, or whether Monte Carlo sampling should proceed)\n
	-s "client_start_residue" (the residue number in the client protein from which to begin fusion)"""
	def __init__(self):
		self.TELSAM_version = "1TEL"
		self.client_pdb = None
		self.client_start_residue = None
		self.linker_variant = None
		self.unit_cell_ab = None
		self.degree_rotation = None
		self.remake_TELSAM_bool = False
		self.optimize = False
		self.exhaustive = False
		self.centroids = False
		self.symmdef_generating =  True
		self.fa_rep_cutoff = 40000.0
		self.headers = {
			"User-Agent": (
			"Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
			"AppleWebKit/537.36 (KHTML, like Gecko) "
			"Chrome/125.0 Safari/537.36"
			)
		}
		try:
			optlist, args = getopt.getopt(sys.argv[1:], "t:c:l:u:d:r:oe:s:")
			for o, a in optlist:
				if o == '-t':
					if a != "1TEL":
						print(f'Error: TELSAM_version argument not supported: {a}')
						sys.exit(1)
					else:
						self.TELSAM_version = a
				elif o == '-c':
					self.client_pdb = a
				elif o == '-s':
					if a!="":
						self.client_start_residue = int(a)
				elif o == '-l':
					if a!="":
						self.linker_variant = int(a)
						self.start_residue_to_superimpose = 17
						if self.linker_variant>=16:
							self.remake_TELSAM_bool = True
							self.start_residue_to_superimpose = 24
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
				elif o == '-e':
					self.exhaustive = True
				else:
					print(f'Unhandled argument: {o}')
					sys.exit(1)
			if args:
				print(f'Unexpected arguments: {args}')
				sys.exit(1)
			if not self.optimize:
				self.centroids = False
			if self.linker_variant is None:
				self.linker_variant = 0
			if self.start_residue_to_superimpose is None:
				self.start_residue_to_superimpose = 17
		except getopt.GetoptError as err:
			print(err)
			sys.exit(1)

		self.base = os.path.join(os.path.dirname(__file__),str(self.linker_variant))
		os.makedirs(self.base,exist_ok=True)
		self.interfaced = False
		self.scores = {True:200000,False:200000}
		self.scores_files = {True:os.path.join(self.base,f'interfaced_scores_file.txt'),False:os.path.join(self.base,f'scores_file.txt')}
		self.min_score_pdbs = {True:None,False:None}

		self.pmm = PyMOLMover()
		self.pmm.keep_history(True)
		self.energies_vs_ucab_vs_deg = {'linker':[],'energy':[],'ucab':[],'deg':[]}
		self.interfaced_energies_vs_ucab_vs_deg = {'linker':[],'energy':[],'ucab':[],'deg':[]}
		self.furthest_x = 0
		self.to_fullatom = SwitchResidueTypeSetMover("fa_standard")
		self.setup_for_refinement()
		self.validate_TELSAM()

	def setup_for_refinement(self):
		"""Assigns self.sf as get_score_function()\n
		Assigns self.relax as FastRelax()\n
		Changes FastRelax settings by running self.relax.set_scorefxn(self.sf)\n
		Creates a task factory and reinitializes from command line and restricts to repacking???\n
		Assigns self.packer as PackRotamersMover(self.sf)
		Changes PackRotamersMover settings by applying the previously-made task factory???\n
		Assigns self.min_mover as MinMover(), applies the score_function to its settings and applies the "lbfgs_armijo_nonmonotone" setting as its min_type.\n
		"""
		#################### SETUP FOR REFINEMENT ##########################
		self.sf = get_score_function()
		self.sf.set_weight(rosetta.core.scoring.coordinate_constraint, 1.0)
		self.sf.set_weight(rosetta.core.scoring.atom_pair_constraint, 1.0)
		self.interface_sf = self.sf.clone()
		self.interface_sf.set_weight(rosetta.core.scoring.coordinate_constraint, 0.0)
		self.interface_sf.set_weight(rosetta.core.scoring.atom_pair_constraint, 0.0)
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

	def add_coordinate_constraints(self, pose, standard_deviation=8.0):
		"""Restrain TELSAM coordinates and preserve client geometry with weak distance restraints."""
		TELS_func = rosetta.core.scoring.func.HarmonicFunc(0.0, 1.0)
		client_func = rosetta.core.scoring.func.FlatHarmonicFunc(0.0, 1.0, standard_deviation)
		reference_atom = AtomID(1, 1)
		for chain in range(1, pose.num_chains() + 1):
			for residue in range(pose.chain_begin(chain), pose.chain_begin(chain)+self.TELSAM_module_end):
				for atom in range(1, pose.residue(residue).natoms() + 1):
					atom_id = AtomID(atom, residue)
					constraint = rosetta.core.scoring.constraints.CoordinateConstraint(
						atom_id,
						reference_atom,
						pose.xyz(atom_id),
						TELS_func,
					)
					pose.add_constraint(constraint)

		for chain in range(1,pose.num_chains()+1):
			for first_residue in range(pose.chain_begin(chain)+self.TELSAM_module_end, pose.chain_end(chain)+1):
				if pose.residue(first_residue).is_virtual_residue():
					continue
				for second_residue in range(first_residue + 7, pose.chain_end(chain)+1):
					if pose.residue(second_residue).is_virtual_residue():
						continue
					first_atom_id = AtomID(pose.residue(first_residue).atom_index("CA"), first_residue)
					second_atom_id = AtomID(pose.residue(second_residue).atom_index("CA"), second_residue)
					distance = pose.xyz(first_atom_id).distance(pose.xyz(second_atom_id))
					client_func = rosetta.core.scoring.func.FlatHarmonicFunc(
						distance,
						1.0,
						standard_deviation,
					)
					constraint = rosetta.core.scoring.constraints.AtomPairConstraint(
						first_atom_id,
						second_atom_id,
						client_func,
					)
					pose.add_constraint(constraint)

	def passes_fa_rep_filter(self, pose):
		"""Return False when the pose has excessive repulsive energy."""
		self.sf(pose)
		fa_rep_energy = pose.energies().total_energies()[fa_rep]*self.sf.get_weight(fa_rep)
		passes = fa_rep_energy <= self.fa_rep_cutoff
		print(f"fa_rep: {fa_rep_energy:.3f} (cutoff: {self.fa_rep_cutoff:.3f})")
		return passes

	def refine(self,pose) -> float:
		"""Creates a movemap, freezing bb and chi angles, then unfreezing the atoms that are outside of the 1TEL module.\n
		Then, changes self.min_mover to have the movemap that was just created.\n
		Then, applies self.min_mover to the current pose.\n
		Then, changes self.relax to have the previous movemap as well.\n
		Then, initializes a TaskFactory that reinitializes from command line and restricts to repacking???\n
		Then, changes self.relax by setting the previously made TaskFactory in it.\n
		Then, applies self.relax to the passed pose.\n
		Finally, returns the energy of the pose from self.sf(pose)"""
		####################### MOVEMAP REFINEMENT (Must be re-setup after pose is symmetrized) #################
		#It may be okay to just have it here and then call the min_mover.movemap and relax.set_movemap functions later.
		movemap = MoveMap()
		movemap.set_bb(False)
		movemap.set_chi(False)
		#The following line specifies the entirety of the client protein as well as the residues that were connected to the client's first alpha helix
		#for every chain. Later, this information is passed into the min_mover as the list of residues that are allowed to move and in what ways they can move.
		#Theoretically, the 1TEL subunit itself should not change much.
		for chain in range(1,pose.num_chains()+1):
			for i in range(pose.chain_end(chain)-self.client.total_residue()-self.start_residue_to_superimpose,pose.chain_end(chain)+1):
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

	def bound_separated_interface_energy(self, pose):
		"""Score chain A against the remaining explicit crystal chains without repacking."""
		explicit_pose = pose.clone()
		if is_symmetric(explicit_pose):
			rosetta.core.pose.symmetry.make_asymmetric_pose(explicit_pose)

		virtual_residues = [
			residue
			for residue in range(1, explicit_pose.total_residue() + 1)
			if explicit_pose.residue(residue).is_virtual_residue()
		]
		for residue in reversed(virtual_residues):
			explicit_pose.delete_residue_slow(residue)

		if explicit_pose.num_chains() < 2:
			raise ValueError("Bound-versus-separated scoring requires at least two chains")

		chains = explicit_pose.split_by_chain()
		chain_a = chains[1]
		other_chains = Pose()
		for chain in range(2, len(chains) + 1):
			append_pose_to_pose(other_chains, chains[chain], True)

		bound_score = self.interface_sf(explicit_pose)
		separated_score = self.interface_sf(chain_a) + self.interface_sf(other_chains)
		interface_energy = bound_score - separated_score
		print(f"Bound score: {bound_score:.3f}")
		print(f"Separated score: {separated_score:.3f}")
		print(f"Bound - separated interface energy: {interface_energy:.3f}")
		return interface_energy

	def interface_refine(self,pose) -> float:
		"""First, creates chain selectors for all chains.\n
		Then, creates NeighborhoodResidueSelectors that accept the chain selectors as arguments...\n
		Then, creates an OrResidueSelector that accepts both NeighborhoodResidueSelectors as arguments...\n
		Then, actually generates a selection by applying the previous selector on the pose.\n
		Creates a movemap_factory that disables movement of the bb and chi angles for all but the interface selector\n
		Creates a TaskFactory to restrict to repacking\n
		Changes self.relax by setting the movemap factory and task factory to it.\n
		Relaxes the passed pose.\n
		Pushes to PyMOL.\n
		Tries to run the InterfaceAnalyzerMover("A_BCDEFG...") and return interface dg. On a failure, instead runs normal 
		refinement and returns the program to non-interfaced mode.\n
		In the future, the migration between non-interfaced and interfaced modes may be a good indicator that we
		should switch into a higher-granularity docking simulation."""
		alphanumeric_dict = {1:"A",2:"B",3:"C",4:"D",5:"E",6:"F",7:"G",8:"H",9:"I",10:"J",11:"K",12:"L",13:"M",14:"N",15:"O",16:"P",
					   17:"Q",18:"R",19:"S",20:"T",21:"U",22:"V",23:"W",24:"X",25:"Y",26:"Z",27:"0",28:"1",29:"2",30:"3",31:"4",32:"5",33:"6",34:"7",35:"8",36:"9",37:"10",38:"11",39:"12",40:"13",41:"14",42:"15",43:"16",44:"17",45:"18",46:"19",47:"20",48:"21",49:"22",50:"23",51:"24",52:"25",53:"26",54:"27",55:"28",56:"29",57:"30"}
		chain_sels = []
		for chain in range(1,pose.num_chains()+1):#walk through the chains in the pose
			print("chain", chain, "begin/end", pose.chain_begin(chain), pose.chain_end(chain), "label", pose.pdb_info().chain(pose.chain_begin(chain)))
			chain_sel = ChainSelector(alphanumeric_dict[chain])#create a chain selector for each chain
			chain_sels.append(chain_sel)#add the chain selector tool to the chain_sels list
		interface_sels = []
		chain_A_neighbor_sel = NeighborhoodResidueSelector(chain_sels[0],4,False)#Create a neighborhood residue selector for chain A
		for chain_sel in chain_sels[1:]:#Walk through all the chain selectors in the chain_sels list, but skip the chain A selector
			neighbor_sel = NeighborhoodResidueSelector(chain_sel, 4, False)#Create a neighborhood residue selector for every chain
			interface_sel = OrResidueSelector(chain_A_neighbor_sel,neighbor_sel)#Create an interface selector; basically, select all atoms that are in both chain A and the current neighborhood selector
			interface_sels.append(interface_sel)#Add each interface selector to the interface_sels list

		# Combine all interfaces
		all_interface_sel = interface_sels[0]#Initialize an interface_selection that includes only the selector for A:B, then...
		for sel in interface_sels[1:]:#Walk through all the other interface_selections
			all_interface_sel = OrResidueSelector(all_interface_sel,sel)#Change the all_interface_sel so that it counts both those previously mentioned and any in the current interface_selection
		module_ranges = []
		for chain in range(1, pose.num_chains() + 1):
			first_module_residue = pose.chain_begin(chain)
			last_module_residue = min(
				first_module_residue + self.TELSAM_module_end - 1,
				pose.chain_end(chain),
			)
			if first_module_residue <= last_module_residue:
				module_ranges.append(f"{first_module_residue}-{last_module_residue}")
		module_selector = ResidueIndexSelector(",".join(module_ranges))
		excluded_selector = OrResidueSelector(module_selector, VirtualResidueSelector())
		all_interface_sel = AndResidueSelector(
			all_interface_sel,
			NotResidueSelector(excluded_selector),
		)
		non_interface_sel = NotResidueSelector(all_interface_sel)#Finally, create a selector that includes all the atoms not in all_interface_sel

		# MoveMap DOESN'T SEEM TO BE WORKING
		movemap_factory = MoveMapFactory()
		movemap_factory.add_bb_action(move_map_action.mm_disable,TrueResidueSelector())
		movemap_factory.add_chi_action(move_map_action.mm_disable,TrueResidueSelector())
		movemap_factory.add_bb_action(move_map_action.mm_enable,all_interface_sel)
		movemap_factory.add_chi_action(move_map_action.mm_enable,all_interface_sel)

		# TaskFactory DOESN'T SEEM TO BE WORKING
		tf = TaskFactory()
		tf.push_back(RestrictToRepacking())
		tf.push_back(OperateOnResidueSubset(PreventRepackingRLT(),non_interface_sel))
		self.relax.set_movemap_factory(movemap_factory)
		self.relax.set_task_factory(tf)
		self.relax.apply(pose)
		try:
			score = self.bound_separated_interface_energy(pose)
			self.interfaced = True
		except Exception as e:
			print(f'Bound-versus-separated scoring failed: {e}')
			score = self.refine(pose)
			self.interfaced = False
			return score

		iam_string = "A_"
		for chain in range(2,pose.num_chains()+1):
			iam_string = iam_string+alphanumeric_dict[chain]
		print(f'IAM_STRING: {iam_string}')
		self.iam = InterfaceAnalyzerMover(iam_string)
		self.iam.set_scorefunction(self.interface_sf)
		self.iam.set_compute_separated_sasa(True)
		self.iam.set_calc_dSASA(True)
		self.iam.set_compute_interface_energy(True)
		try:
			self.iam.apply(pose)
			fixed_chains = self.iam.get_fixed_chains()
			print(f'fixed_chains: {fixed_chains}')
			interface_dG = self.iam.get_interface_dG()
			print(f'InterfaceAnalyzerMover interface dG: {interface_dG}')
			cenergy = self.iam.get_complex_energy()
			print(f'Complex Energy: {cenergy}')
			csasa = self.iam.get_complexed_sasa()
			print(f'csasa: {csasa}')
			dsasa = self.iam.get_interface_delta_sasa()
			print(f'dsasa: {dsasa}')
			interface_set = self.iam.get_interface_set()
			print(f'interface_set: {interface_set}')
		except Exception as e:
			print(f'InterfaceAnalyzerMover diagnostics failed: {e}')
		return score

	def get_CRYST1(self,pdb):
		"""Accepts a PDB file and returns the CRYST1 line from that file as a tuple with (a,b,c,alpha,beta,gamma)."""
		with open(pdb, 'r') as s:
				for line in s:
					if "CRYST1" in line:
						p = re.compile(r'\d+\.\d+')
						cryst1_vals = p.findall(line)
						a, b, c, alpha, beta, gamma = [float(x) for x in cryst1_vals[0:6]]
						return (a, b, c, alpha, beta, gamma)

	def add_CRYST1(self,new_pdb,old_pdb):
		"""Copies CRYST1 data from the old_pdb into the new_pdb."""
		with open (os.path.join(self.base,new_pdb),'r') as file:
			pdb_sans_cryst = file.read()
		with open(os.path.join(self.base,new_pdb),'w') as file:
			with open(os.path.join(self.base,old_pdb)) as s:
				for line in s:
					if "CRYST1" in line:
						file.write(line)
						break
			file.write(pdb_sans_cryst)
		
	def change_cell(self,read_file,write_file,wa=None,wb=None,wc=None):
		"""Copies the read_file PDB into a new write_file PDB at the given path, supplying it with an altered CRYST1 line as the user
		specifies. wa becomes the new a, wb becomes the new b, wc becomes the new c. The rest of the PDB remains unchanged."""
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
						if wc!=None:
							c = wc
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
	def rotate_pose(self,pose,deg):
		"""Rotates a pose about the global z-axis."""
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
		pose.apply_transform_Rx_plus_v(R, v)
		return pose
		
	def chart(self,linker:int):
		"""Primarily makes a heat map of energies related to unit cell ab and degree of rotation for each linker length variant.\n
		Can also be used to make a 3d graph of the relationship between energy, unit cell ab, and degree for each linker length variant.\n
		Relies heavily on matplotlib.
		"""
		fig = plt.figure()
		data = {
			"aboi":[ucab for ucab, l in zip(self.interfaced_energies_vs_ucab_vs_deg['ucab'],self.interfaced_energies_vs_ucab_vs_deg['linker']) if l == linker],
			"doi":[deg for deg, l in zip(self.interfaced_energies_vs_ucab_vs_deg['deg'],self.interfaced_energies_vs_ucab_vs_deg['linker']) if l == linker],
			"eoi":[energy for energy, l in zip(self.interfaced_energies_vs_ucab_vs_deg['energy'],self.interfaced_energies_vs_ucab_vs_deg['linker']) if l == linker]
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
		ax.set_title(f'Energies of UCAB:Degree Combinations for {self.TELSAM_version}--{self.client_pdb}_{linker}')
		fig.savefig(os.path.join(self.base,f"Energies of UCAB_Degree Combinations for {self.TELSAM_version}--{self.client_pdb}_{linker}"))
		with open(os.path.join(self.base,f"{linker}_chart.json"),"w") as file:
			json.dump(data,file,indent=4)
		
		#print(data["aboi"][::],data["doi"][::],data["eoi"][::])
		ax = fig.add_subplot(projection='3d')
		ax.scatter(data["aboi"],data["doi"],data["eoi"])
		ax.set_title(f'Energies of UCAB:Degree Combinations for {self.TELSAM_version}--{self.client_pdb}_{linker}')
		ax.set_xlabel('Unit Cell AB Length (Angstroms)')
		ax.set_ylabel('Degree of Polymer Rotation (Degrees)')
		ax.set_zlabel('Energy (REU)')
		fig.savefig(os.path.join(self.base,f"Energies of UCAB_Degree Combinations for {self.TELSAM_version}--{self.client_pdb}_{linker}"))
		plt.close(fig)

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

		if self.linker_variant>=16:
			####################################### EXTEND HELIX ###############################################
			#Extend TELSAM's helix by 7 amino acids.
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
			mutate_residue(TELSAM_in_9DOC,90,"A",5)
			mutate_residue(TELSAM_in_9DOC,92,"K",5)
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
			else:
				client_path = self.client_pdb
				self.client_pdb = os.path.basename(client_path)
				subprocess.run(["cp",client_path,os.path.join(self.base,f"{self.client_pdb}.pdb")],check=True)
			cleanATOM(os.path.join(self.base,f"{self.client_pdb}.pdb"))
			temp_pose = pose_from_pdb(os.path.join(self.base,f"{self.client_pdb}.clean.pdb"))
			os.remove(os.path.join(self.base,f"{self.client_pdb}.pdb"))
			os.remove(os.path.join(self.base,f"{self.client_pdb}.clean.pdb"))
			
			if self.client_start_residue is None:
				#Extract first 4-aa helical region from target protein to fuse to TELSAM:
				dssp = Dssp(temp_pose)
				dssp.insert_ss_into_pose(temp_pose)
				ss_string = temp_pose.secstruct()
				first_helix = ss_string.find("HHHHH")
				append_subpose_to_pose(self.client,temp_pose,temp_pose.chain_begin(1)+first_helix,temp_pose.chain_end(1))
			else:
				append_subpose_to_pose(self.client,temp_pose,temp_pose.chain_begin(1)+self.client_start_residue-1,temp_pose.chain_end(1))
			client_start = self.client.chain_begin(1)
			if self.client.residue(client_start).has_variant_type(rosetta.core.chemical.LOWERTERM_TRUNC_VARIANT):
				rosetta.core.pose.remove_variant_type_from_pose_residue(
					self.client,
					rosetta.core.chemical.LOWERTERM_TRUNC_VARIANT,
					client_start,
				)
			#Align the two helices:
			if os.path.exists(os.path.join(self.base,f'scores_file.txt')):
				os.remove(os.path.join(self.base,f'scores_file.txt'))
			if os.path.exists(os.path.join(self.base,f'interfaced_scores_file.txt')):
				os.remove(os.path.join(self.base,f'interfaced_scores_file.txt'))
			self.TELSAM = self.TELSAM_in_9DOC.clone()
			self.start_residue_to_superimpose-=self.linker_variant
			self.TELSAM_module_end = self.TELSAM.chain_end(1)-self.start_residue_to_superimpose-1
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
			self.add_coordinate_constraints(self.TELSAM)

			#Refine
			movemap = MoveMap()
			movemap.set_bb(False)
			movemap.set_chi(False)
			for i in range(self.TELSAM.chain_end(1)-self.client.total_residue()-self.start_residue_to_superimpose,self.TELSAM.chain_end(1)+1):
				movemap.set_bb(i, True)
				movemap.set_chi(i, True)
			#Refine gently
			self.min_mover.movemap(movemap)
			#self.packer.apply(self.TELSAM)
			self.min_mover.apply(self.TELSAM)
			self.relax.set_movemap(movemap)
			self.relax.apply(self.TELSAM)
			if not self.passes_fa_rep_filter(self.TELSAM):
				raise RuntimeError("Initial fusion rejected by the fa_rep filter")

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
			self.linker_pdb = os.path.join(self.base,f'{self.TELSAM_version}--{self.client_pdb}_{self.linker_variant}.pdb')
			self.TELSAM.dump_pdb(self.linker_pdb)
			self.add_CRYST1(os.path.basename(self.linker_pdb),os.path.basename("TELSAM_in_9DOC.pdb"))
			with open(self.linker_pdb, 'r') as file:
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
					scores.write(f'Linker file: {self.linker_pdb}\n')
					scores.write(f'score: {score}\n')
					pdb_path = os.path.join(self.base,f'{self.TELSAM_version}--{self.client_pdb}_{self.linker_variant}.pdb')
					sequence = self.TELSAM.sequence()	
					with open (f'{pdb_path.removesuffix('.pdb')}.fasta', 'w') as f:
						f.write(">"+pdb_path+", score (REU): "+"{:.3e}".format(score)+"\n"+"HHHHHHHHHH"+str(sequence).strip('X'))
			self.TELSAM.pdb_info().name("Fusion")
			self.pmm.apply(self.TELSAM)
			self.pmm.send_energy(self.TELSAM)

		except Exception as e:
			print(e,file=sys.stderr)
			sys.exit()

	###This is where I'm replacing the decision tree with Monte Carlo methodology.
	def monte_carlo_interface_minimization(self):
		"""monte_carlo_interface_minimization runs a Monte Carlo algorithm to change the degree of rotation of the pose\
about the unit cell C-axis and to change the distance between subunits by altering the unit cell AB distance.\
It does this semi-randomly, with a gaussian step in degree of rotation (centered around 0, st.dev of 10)\
(many rotations can be tested, and they are tested widely, with little preference for one or another),\
and a gaussian step in unit cell size (centered around -2.5, st. dev 10) (most attempts are to shrink the unit cell)\n
Currently, the kT of this Monte Carlo method is set at 3.3. Notably, the method is not foolproof. If the lowest achievable energy is relatively\
close to higher interfaced energies, you will not be able to isolate it. It is not unusual for the score to dip and rise continuously, unless\
a clear energy well can be accessed and the gaussian steps and kT are small enough to limit the pose's movement to the inside of that energy well."""
		ucab = self.get_CRYST1(self.linker_pdb)[0]
		for u_sample in range(15):
			test_ucab = round(ucab + random.gauss(-5,5),3)
			ucab_pdb = os.path.join(self.base,f'{self.TELSAM_version}--{self.client_pdb}_{self.linker_variant}_{test_ucab}.pdb')
			self.change_cell(self.linker_pdb,ucab_pdb,wa=test_ucab,wb=test_ucab)
			symm_pose = pose_from_pdb(ucab_pdb)
			makesym = SetupForSymmetryMover("CRYST1")
			makesym.apply(symm_pose)
			self.add_coordinate_constraints(symm_pose)
			deg = 0
			starting_energy = self.interface_refine(symm_pose)#this has to come before self.interfaced, because it changes self.interfaced based on whether the interface scoring works.
			if not self.passes_fa_rep_filter(symm_pose):
				ucab=ucab+3
				u_sample = u_sample-1
				symm_pose.pdb_info().name("failed_ucab_pose")
				self.pmm.apply(symm_pose)
				self.pmm.send_energy(symm_pose)
				continue
			if self.interfaced:
				for d_sample in range(4):
					test_deg = deg + random.gauss(0,10)
					symm_pose = pose_from_pdb(ucab_pdb)
					symm_pose.pdb_info().name("pmm")
					self.rotate_pose(symm_pose,test_deg)
					makesym = SetupForSymmetryMover("CRYST1")
					makesym.apply(symm_pose)
					self.add_coordinate_constraints(symm_pose)
					energy = self.interface_refine(symm_pose)#this has to come before self.interfaced, because it changes self.interfaced based on whether the interface scoring works.
					if not self.passes_fa_rep_filter(symm_pose):
						symm_pose.pdb_info().name("failed_degree_pose")
						self.pmm.apply(symm_pose)
						self.pmm.send_energy(symm_pose)
						d_sample = d_sample-1
						continue
					symm_pose.dump_pdb(os.path.join(self.base,f'{self.TELSAM_version}--{self.client_pdb}_{self.linker_variant}_{test_ucab}_{test_deg}_symmetric.pdb'))
					if self.interfaced:
						self.pmm.apply(symm_pose)
						self.pmm.send_energy(symm_pose)
						delta_e = energy - self.scores[self.interfaced]
						if delta_e<=0 or random.random()<math.exp(-delta_e/3.3):
							self.min_score_pdbs[self.interfaced] = os.path.join(self.base,f'{self.TELSAM_version}--{self.client_pdb}_{self.linker_variant}_{test_ucab}_{test_deg}.pdb')
							self.scores[self.interfaced] = energy
							ucab = test_ucab
							deg = test_deg
							#RECORD FOR LATER
							self.interfaced_energies_vs_ucab_vs_deg['linker'].append(int(self.linker_variant))
							self.interfaced_energies_vs_ucab_vs_deg["ucab"].append(int(float(test_ucab)))
							self.interfaced_energies_vs_ucab_vs_deg["deg"].append(int(float(test_deg)))
							self.interfaced_energies_vs_ucab_vs_deg['energy'].append(float(energy))
							scores_file = self.scores_files[self.interfaced]
							with open(scores_file,'a') as scores:
								scores.write(f'{self.linker_variant}_{test_ucab}_{test_deg}\n')
								scores.write(f'Score: {energy}\n')
					else:
						symm_pose.pdb_info().name("failed_interface_pose")
						self.pmm.apply(symm_pose)
						self.pmm.send_energy(symm_pose)
						d_sample = d_sample-1
			else:
				u_sample = u_sample-1
		#self.chart(self.linker_variant)

	def picker(self):
		ucab_pdb = os.path.join(self.base,f'{self.TELSAM_version}--{self.client_pdb}_{self.linker_variant}_{str(int(float(self.unit_cell_ab)))}.pdb')
		self.change_cell(self.linker_pdb,ucab_pdb,wa=self.unit_cell_ab,wb=self.unit_cell_ab)
		symm_pose = pose_from_file(ucab_pdb)
		self.rotate_pose(symm_pose,self.degree_rotation)
		sequence = symm_pose.sequence()
		self.interface_refine(symm_pose)
		score = self.sf(symm_pose)
		with open (f'{ucab_pdb.removesuffix('.pdb')}.fasta', 'w') as f:
			f.write(">"+ucab_pdb+", score (REU): "+"{:.3e}".format(score)+"\n"+"HHHHHHHHHH"+str(sequence).strip('X'))
		self.add_CRYST1(os.path.basename(ucab_pdb),os.path.basename(ucab_pdb))
		symm_pose.pdb_info().name("pmm")
		self.pmm.apply(symm_pose)
		self.pmm.send_energy(symm_pose)
		self.chart(self.linker_variant)

def main():
	TELSetta1 = TELSetta()
	if TELSetta1.TELSAM_version=="1TEL" and TELSetta1.client_pdb!=None:
		TELSetta1.fuse()
		if bool(TELSetta1.unit_cell_ab==None and TELSetta1.degree_rotation==None) or TELSetta1.optimize==True:
			TELSetta1.monte_carlo_interface_minimization()
		elif TELSetta1.unit_cell_ab!=None and TELSetta1.degree_rotation!=None and TELSetta1.optimize==False:
			TELSetta1.picker()
		else:
			print(f"Not enough arguments were provided to pick a specific TELSAM--fusion variant to model.")
	else:
		print(f"You didn't pass any client proteins to fuse to TELSAM.")

main()