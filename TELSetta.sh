#The default command to run this program on a given client is:
#bash TELSetta.sh -c "7TCY" -r true -l "1"
#This remakes the 1TEL subunit, fuses it to your client once, and tests the different alignments of the polymers to see which is best.

TELSAM_version="1TEL"
pymol_setting="true"
linker_variant=""
degree_rotation=""
remake_TELSAM="false"
optimize_TELSAM="false"
exhaustive="false"
max_parallel_jobs="${TELSETTA_MAX_JOBS:-2}"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
gene_designer="$script_dir/GeneDesigner2.exe"

if [[ ! "$max_parallel_jobs" =~ ^[1-9][0-9]*$ ]]; then
    echo "TELSETTA_MAX_JOBS must be a positive integer: $max_parallel_jobs" >&2
    exit 1
fi

if [ ! -x "$gene_designer" ]; then
    echo "GeneDesigner executable not found or not executable: $gene_designer" >&2
    exit 1
fi

while getopts "pt:c:s:l:u:d:r:oe" flag
do
    case "${flag}" in
        p) pymol_setting="false";;
        t) TELSAM_version="${OPTARG}";;
        c) client="${OPTARG}";;
        s) client_start_residue="${OPTARG}";;
        l) linker_variant="${OPTARG}";;
        u) unit_cell_ab="${OPTARG}";;
        d) degree_rotation="${OPTARG}";;
        r) remake_TELSAM="${OPTARG}";;
        o) optimize_TELSAM="true";;
        e) exhaustive="true";;
        \?) echo "Invalid option: -$OPTARG" >&2; exit 1;;
    esac
done

client_name="${client##*/}"

cmd_base=(
    python ~/TELSAM-Fusion-Crystallography-with-Rosetta/start_TELSetta.py \
    -t "$TELSAM_version" \
    -c "$client" \
    -s "$client_start_residue" \
    -u "$unit_cell_ab" \
    -d "$degree_rotation" \
    -r "$remake_TELSAM"
)

#If no linker variant of interest is provided, test all of them in parallel without posting them to PyMOL. Then, with the lowest-energy combinations from
#each of the fourteen linker variants, run a symmetric refinement and save the pdb.
if [ "$linker_variant" = "" ]; then
    echo "Running start_TELSetta for each linker variant (0/16)"

    for linker_variant in {0..16}; do
        (
            echo "Launching $linker_variant/16"
            cmd=("${cmd_base[@]}" -l "$linker_variant" -o)
            if [ "$exhaustive" = "true" ]; then
                cmd=("${cmd_base[@]}" -l "$linker_variant" -o -e)              
            fi
            "${cmd[@]}"
            file="$HOME/TELSAM-Fusion-Crystallography-with-Rosetta/${linker_variant}/${linker_variant}_chart.json"
            result=$(jq '
                [
                    range(0; (.eoi | length)) as $i
                    | {
                        eoi: .eoi[$i],
                        aboi: .aboi[$i],
                        doi: .doi[$i]
                        }
                ]
                | min_by(.eoi)
                ' "$file")
            min_e=$(echo "$result" | jq -r '.eoi')
            min_ab=$(echo "$result" | jq -r '.aboi')
            min_d=$(echo "$result" | jq -r '.doi')
            mcmd=(
                python ~/TELSAM-Fusion-Crystallography-with-Rosetta/start_TELSetta.py \
                -t "$TELSAM_version" \
                -c "$client" \
                -s "$client_start_residue" \
                -u "$min_ab" \
                -d "$min_d" \
                -r "$remake_TELSAM" \
                -l "$linker_variant" 
            )
            "${mcmd[@]}"
            fasta="$HOME/TELSAM-Fusion-Crystallography-with-Rosetta/${linker_variant}/${TELSAM_version}--${client_name}_${linker_variant}.fasta"
            fastout="$HOME/TELSAM-Fusion-Crystallography-with-Rosetta/${TELSAM_version}--${client_name}_${linker_variant}_gene.fasta"
            echo "fasta:$fasta fastout:$fastout"
            "$gene_designer" "$fasta" "$fastout" "None"
        ) &

        while (( $(jobs -rp | wc -l) >= max_parallel_jobs )); do
            wait -n
        done

    done

    wait

#If a linker variant is provided, connect to PyMOL. Run the stepper program if optimize is enabled. Otherwise, run the picker program.
else
    if [ "$pymol_setting" = "true" ]; then
        pymol ~/TELSAM-Fusion-Crystallography-with-Rosetta/start_pymol_server.pml &
        sleep 2
    fi

    cmd=("${cmd_base[@]}" -l "$linker_variant")
    if [ "$optimize_TELSAM" = "true" ]; then
        cmd=("${cmd_base[@]}" -l "$linker_variant" -o)
        if [ "$exhaustive" = "true" ]; then
            cmd=("${cmd_base[@]}" -l "$linker_variant" -o -e)
        fi
    else
        if [ "$exhaustive" = "true" ]; then
            cmd=("${cmd_base[@]}" -l "$linker_variant" -e)
        fi
    fi
    "${cmd[@]}"
    fasta="$HOME/TELSAM-Fusion-Crystallography-with-Rosetta/${linker_variant}/${TELSAM_version}--${client_name}_${linker_variant}.fasta"
    fastout="$HOME/TELSAM-Fusion-Crystallography-with-Rosetta/${linker_variant}/${TELSAM_version}--${client_name}_${linker_variant}_gene.fasta"
    echo "fasta:$fasta fastout:$fastout"
    "$gene_designer" "$fasta" "$fastout" "None"
fi