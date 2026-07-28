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

while getopts "pt:c:l:u:d:r:oe" flag
do
    case "${flag}" in
        p) pymol_setting="false";;
        t) TELSAM_version="${OPTARG}";;
        c) client="${OPTARG}";;
        l) linker_variant="${OPTARG}";;
        u) unit_cell_ab="${OPTARG}";;
        d) degree_rotation="${OPTARG}";;
        r) remake_TELSAM="${OPTARG}";;
        o) optimize_TELSAM="true";;
        e) exhaustive="true";;
        \?) echo "Invalid option: -$OPTARG" >&2; exit 1;;
    esac
done

cmd_base=(
    python ~/TELSAM-Fusion-Crystallography-with-Rosetta/start_TELSetta.py \
    -t "$TELSAM_version" \
    -c "$client" \
    -u "$unit_cell_ab" \
    -d "$degree_rotation" \
    -r "$remake_TELSAM"
)

#If no linker variant of interest is provided, test all of them in parallel without posting them to PyMOL. Then, with the lowest-energy combinations from
#each of the fourteen linker variants, run a symmetric refinement and save the pdb.
if [ "$linker_variant" = "" ]; then
    echo "Running start_TELSetta for each linker variant (0/14)"

    for linker_variant in {1..14}; do
        (
            echo "Launching $linker_variant/14"
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
                -u "$min_ab" \
                -d "$min_d" \
                -r "$remake_TELSAM" \
                -l "$linker_variant" 
            )
            "${mcmd[@]}"
            #fasta="$HOME/TELSAM-Fusion-Crystallography-with-Rosetta/${linker_variant}/${TELSAM_version}--${client}_${linker_variant}_${min_ab}_${min_d}.fasta"
            #fastout="$HOME/TELSAM-Fusion-Crystallography-with-Rosetta/${TELSAM_version}--${client}_${linker_variant}_${min_ab}_${min_d}_gene.fasta"
            #echo "fasta:$fasta fastout:$fastout"
            #GeneDesigner2.exe "$fasta" "$fastout" "None"
        ) &

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
    #fasta="$HOME/TELSAM-Fusion-Crystallography-with-Rosetta/${linker_variant}/${TELSAM_version}--${client}_${linker_variant}_${min_ab}_${min_d}.fasta"
    #fastout="$HOME/TELSAM-Fusion-Crystallography-with-Rosetta/${TELSAM_version}--${client}_${linker_variant}_${min_ab}_${min_d}_gene.fasta"
    #echo "fasta:$fasta fastout:$fastout"
    #GeneDesigner2.exe "$fasta" "$fastout" "None"
fi