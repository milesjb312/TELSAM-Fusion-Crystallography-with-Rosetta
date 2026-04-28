TELSAM_version="1TEL"
pymol_setting="true"
linker_variant=""
degree_rotation=""
remake_TELSAM="false"
optimize_TELSAM="false"

while getopts "p:t:c:l:u:d:r:o" flag
do
    case "${flag}" in
        p) pymol_setting="${OPTARG}";;
        t) TELSAM_version="${OPTARG}";;
        c) client="${OPTARG}";;
        l) linker_variant="${OPTARG}";;
        u) unit_cell_ab="${OPTARG}";;
        d) degree_rotation="${OPTARG}";;
        r) remake_TELSAM="${OPTARG}";;
        o) optimize_TELSAM="true";;
        \?) echo "Invalid option: -$OPTARG" >&2; exit 1;;
    esac
done

cmd_base=(
    python ~/TELSetta/start_TELSetta.py \
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
            "${cmd[@]}"
            file="$HOME/TELSetta/${linker_variant}/${linker_variant}_chart.json"
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
            min_e=$(echo "$result" | jq '.eoi')
            min_ab=$(echo "$result" | jq '.aboi')
            min_d=$(echo "$result" | jq '.doi')
            mcmd=(
                python ~/TELSetta/start_TELSetta.py \
                -t "$TELSAM_version" \
                -c "$client" \
                -u "$min_ab" \
                -d "$min_d" \
                -r "$remake_TELSAM" \
                -l "$linker_variant" 
            )
            "${mcmd[@]}"
            fasta="$HOME/TELSetta/${linker_variant}/${TELSAM_version}--${client}_${linker_variant}_${min_ab}_${min_d}.fasta"
            fastout="$HOME/TELSetta/${TELSAM_version}--${client}_${linker_variant}_${min_ab}_${min_d}_gene.txt"
            echo "fasta:$fasta fastout:$fastout"
            GeneDesigner2.exe $fasta $fastout "None"
        ) &

    done

    wait

#If a linker variant is provided, connect to PyMOL. Run the stepper program if optimize is enabled. Otherwise, run the picker program.
else
    if [ "$pymol_setting" = "true" ]; then
        pymol ~/TELSetta/start_pymol_server.pml &
        sleep 2
    fi

    if [ $optimize_TELSAM = "true" ]; then
        cmd=("${cmd_base[@]}" -l "$linker_variant" -o)
    else
        cmd=("${cmd_base[@]}" -l "$linker_variant")
    fi
    "${cmd[@]}"
fi