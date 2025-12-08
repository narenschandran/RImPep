CURR_SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")
conf_file="${CURR_SCRIPT_DIR}/../.conf.sh"

if [ ! -f "$conf_file" ]; then
    echo "Unable to find configuration file at: [$conf_file]" >&2
    echo "Exiting..." >&2
    kill -9 $$
else
    source "$conf_file"
fi

RUNF="${PROJROOT}/run.sh"
for d in "${PROJROOT}/results/03-multiallele-benchmark/"*; do
    if [ -d "$d" ]; then
        echo $(basename "$d")
        bash "$RUNF" "$d"
    fi
done
