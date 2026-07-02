#The default command to run this program on a given client is:
#bash TELSetta.sh -c "7TCY" -r true -l "1"
#This remakes the 1TEL subunit, fuses it to your client once, and tests the different alignments of the polymers to see which is best.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${TELSETTA_DATA_DIR:-$SCRIPT_DIR/data}"
mkdir -p "$DATA_DIR"
TELSAM_version="1TEL"
pymol_setting="true"
linker_variant=""
degree_rotation=""
remake_TELSAM="false"
optimize_TELSAM="false"
TELSAM_url=""

cleanup_port() {
    pkill -f "PyMOLRosettaServer.py" 2>/dev/null
    pkill -f "pymol.*start_pymol_server" 2>/dev/null
    sleep 1
}

while getopts "p:t:c:l:u:d:r:z:o" flag
do
    case "${flag}" in
        p) pymol_setting="${OPTARG}";;
        t) TELSAM_version="${OPTARG}";;
        c) client="${OPTARG}";;
        l) linker_variant="${OPTARG}";;
        u) unit_cell_ab="${OPTARG}";;
        d) degree_rotation="${OPTARG}";;
        r) remake_TELSAM="${OPTARG}";;
        z) TELSAM_url="${OPTARG}";;
        o) optimize_TELSAM="false";;
        \?) echo "Invalid option: -$OPTARG" >&2; exit 1;;
    esac
done

cmd_base=(
    python "$SCRIPT_DIR/start_TELSetta.py"\
    -t "$TELSAM_version" \
    -c "$client" \
    -u "$unit_cell_ab" \
    -d "$degree_rotation" \
    -r "$remake_TELSAM" \
    -z "$TELSAM_url"
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
            file="$DATA_DIR/${linker_variant}/${linker_variant}_chart.json"
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
                python "$SCRIPT_DIR/start_TELSetta.py" \
                -t "$TELSAM_version" \
                -c "$client" \
                -u "$min_ab" \
                -d "$min_d" \
                -r "$remake_TELSAM" \
                -l "$linker_variant" 
            )
            "${mcmd[@]}"
            fasta="${DATA_DIR}/${linker_variant}/${TELSAM_version}--${client}_${linker_variant}_${min_ab}_${min_d}.fasta"
            fastout="${DATA_DIR}/${TELSAM_version}--${client}_${linker_variant}_${min_ab}_${min_d}_gene.fasta"
            echo "fasta:$fasta fastout:$fastout"
            GeneDesigner2.exe "$fasta" "$fastout" "None"
        ) &

    done

    wait

#If a linker variant is provided, connect to PyMOL. Run the stepper program if optimize is enabled. Otherwise, run the picker program.
else
    if [ "$pymol_setting" = "true" ]; then
        cleanup_port
        cd "$SCRIPT_DIR" || exit
        pymol "start_pymol_server.pml" > /tmp/pymol_log.txt 2>&1 &
        cd - > /dev/null
        sleep 2
    fi

    if [ $optimize_TELSAM = "true" ]; then
        cmd=("${cmd_base[@]}" -l "$linker_variant" -o)
    else
        cmd=("${cmd_base[@]}" -l "$linker_variant")
    fi
    "${cmd[@]}"
    fasta="${DATA_DIR}/${linker_variant}/${TELSAM_version}--${client}_${linker_variant}_${min_ab}_${min_d}.fasta"
    fastout="${DATA_DIR}/${TELSAM_version}--${client}_${linker_variant}_${min_ab}_${min_d}_gene.fasta"
    echo "fasta:$fasta fastout:$fastout"
    wine "${SCRIPT_DIR}/bin/GeneDesigner2.exe" "$fasta" "$fastout" "None"
fi