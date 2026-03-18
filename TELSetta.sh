TELSAM_Type=$1
Client_Protein=$2

pymol ~/TELSetta/start_pymol_server.pml &
sleep 2

python ~/TELSetta/start_TELSetta.py "$TELSAM_Type" "$Client_Protein"

mapfile -t pdbs < "$HOME/TELSetta/pdb_list.txt"
i=0
for pdb in "${pdbs[@]}"; do
    echo "The current PDB is: $pdb"
    symm="$HOME/TELSetta/${Client_Protein}_${i}.symm"
    make_symmdef_file.pl -p "$pdb" -m cryst>"$symm"
    echo symm: $symm
    python ~/TELSetta/TELSetta_pmm.py "$pdb" "$symm"
    rm "$pdb"
    ((i++))
done