TELSAM_version="1TEL"
pymol_setting="true"
linker_size=""
degree_rotation=""
remake_TELSAM="false"

while getopts p:t:c:l:u:d:r: flag
do
    case "${flag}" in
        p) pymol_setting="${OPTARG}";;
        t) TELSAM_version="${OPTARG}";;
        c) client="${OPTARG}";;
        l) linker_size="${OPTARG}";;
        u) unit_cell_ab="${OPTARG}";;
        d) degree_rotation="${OPTARG}";;
        r) remake_TELSAM="${OPTARG}";;
        \?) echo "Invalid option: -$OPTARG" >&2; exit 1;;
    esac
done

if [ "$pymol_setting" = "true" ]; then
    pymol ~/TELSetta/start_pymol_server.pml &
    sleep 2
fi

if [ "$TELSAM_version" = "1TEL" ]; then
    python ~/TELSetta/start_TELSetta.py \
    -t "$TELSAM_version" \
    -c "$client" \
    -l "$linker_size" \
    -u "$unit_cell_ab" \
    -d "$degree_rotation" \
    -r "$remake_TELSAM"
fi