#Docker_Tutorial
from pyrosetta import *
from pyrosetta.rosetta.protocols.docking import setup_foldtree
import pyrosetta.rosetta.protocols.rigid as rigid_moves

pyrosetta.init()
pose = pose_from_file("/home/milesjb/TELSetta/1/1TEL--7TCY_1.pdb")
starting_pose = pose.clone()
cen_pose = pose.clone()
cen_switch = SwitchResidueTypeSetMover("centroid")
cen_switch.apply(cen_pose)
starting_cen_pose = cen_pose.clone()
print(pose.fold_tree())
setup_foldtree(pose, "A_B", Vector1([1]))
print(pose.fold_tree)
print(pose.jump(1).get_rotation())
print(pose.jump(1).get_translation())
pert_mover = rigid_moves.RigidBodyPerturbMover(1,8,3)
