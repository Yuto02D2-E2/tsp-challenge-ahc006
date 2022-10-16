#! /bin/sh

set -e
set -u

function compile() {
    # FIXME: clangが使えない
    g++ solution/main.cpp -o solution/main.bin -std=c++17 -Wall -Wextra -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC $1
    # clang++ solution/main.cpp -o solution/main.bin -std=c++17 -Wall -Wextra -fsanitize=undefined,address -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC $1
    return;
}

function run() {
    echo "[run mode]"
    f="${3:-0000}"
    echo -e "\n> test case:$f"
    if [ "$2" = "python" ]; then
        python solution/main.py < given_tools/in/"$f".txt | tee out/python/"$f".txt
        cat out/python/"$f".txt | clip
    else
        compile ""
        solution/main.bin < given_tools/in/"$f".txt | tee out/cpp/"$f".txt
        cat out/cpp/"$f".txt | clip
    fi
    return;
}

function debug() {
    echo "[debug mode]"
    f="${3:-0000}"
    echo -e "\n> test case:$f"
    if [ "$2" = "python" ]; then
        python solution/main.py debug < given_tools/in/"$f".txt | tee out/python/"$f".txt
        cat out/python/"$f".txt | clip
    else
        compile "-D_LOCAL -g"
        solution/main.bin < given_tools/in/"$f".txt | tee out/cpp/"$f".txt
        cat out/cpp/"$f".txt | clip
    fi
    return;
}

function test() {
    echo "[test mode]"
    if [ "$2" = "cpp" ]; then
        compile "-D_LOCAL"
    fi
    read -p "target script:" src_
    python solution/analysis.py --src $src_
    return;
}

function save() {
    echo "[save mode]"
    now="$(date '+%Y-%m-%d-%H-%M-%S')"
    echo "save at $now"
    mkdir archives/$now
    cp -r solution/ archives/$now/
    cp -r tests/ archives/$now/
    return;
}

if [ $# = 0 ]; then
    echo "Usage: sh $0 [run|debug|test|save] [python|cpp] <test case e.g. 0000>"
    exit 0
fi

mode="$1"
if [ "$mode" = "run" ]; then
    run $@
    break
elif [ "$mode" = "debug" ]; then
    debug $@
    break
elif [ "$mode" = "test" ]; then
    test $@
    break
elif [ "$mode" = "save" ]; then
    save
    break
else
    echo "Usage: sh $0 [run|debug|test|save] [python|cpp] <test case e.g. 0000>"
fi
