TELSAM_Type=$1
Client_Protein=$2

pymol ~/TELSetta/start_pymol_server.pml &
sleep 2

python -i ~/TELSetta/start_pyrosetta.py $TELSAM_Type $Client_Protein

