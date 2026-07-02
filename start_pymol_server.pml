python
import pyrosetta, os
server_path = os.path.join(os.path.dirname(pyrosetta.__file__), "PyMOLRosettaServer.py")
cmd.do(f'run {server_path}')
python end