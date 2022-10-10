#! /bin/sh

set -e
set -u

function run() {
    echo "[run mode]"
    f="${2:-0000}"
    echo -e "\n> test case:$f"
    python main.py < given_tools/in/"$f".txt
    python main.py < given_tools/in/"$f".txt | clip
}

function debug() {
    echo "[debug mode]"
    f="${2:-0000}"
    echo -e "\n> test case:$f"
    python main.py debug < given_tools/in/"$f".txt
    python main.py < given_tools/in/"$f".txt | clip
}

function save() {
    echo "[save mode]"
    echo "save at $(date '+%Y-%m-%d-%H-%M-%S')"
    cat main.py > "archive/$(date '+%Y-%m-%d-%H-%M-%S').py"
}

if [ $# = 0 ]; then
    echo "Usage: sh $0 [run|debug|save] <option>"
    exit 0
fi

mode="$1"
if [ "$mode" = "run" ]; then
    run $@
    break
elif [ "$mode" = "debug" ]; then
    debug $@
    break
elif [ "$mode" = "save" ]; then
    save
    break
else
    echo "Usage: sh $0 [run|debug|save] <option>"
fi
