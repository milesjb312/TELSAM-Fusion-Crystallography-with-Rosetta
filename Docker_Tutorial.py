#Docker_Tutorial
from pyrosetta import *
from pyrosetta.rosetta.protocols.docking import setup_foldtree

pyrosetta.init()
pose = pose_from_file("/home/milesjb/TELSetta/1/1TEL--7TCY_1.pdb")
starting_pose = pose.clone()
cen_pose = pose.clone()
cen_switch = SwitchResidueTypeSetMover("centroid")
cen_switch.apply(cen_pose)
starting_cen_pose = cen_pose.clone()
print(pose.fold_tree())
setup_foldtree(pose, "A_B", Vector1([1]))