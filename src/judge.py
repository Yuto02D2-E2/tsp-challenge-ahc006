import argparse
import os
import subprocess


def compile(src: str) -> None:
    # c++fileのコンパイル
    base, ext = src.split(".")
    assert ext in ["cpp", "cc"], (src, ext)
    base = os.path.join("solution", base)
    src = os.path.join("solution", src)
    cmd_ = [
        "clang++",
        f"{src}",
        "-o",
        f"{base}.bin",
        "-std=c++17",
        "-Wall",
        "-Wextra",
        "-fsanitize=undefined,address",
        "-D_GLIBCXX_DEBUG",
        "-D_GLIBCXX_DEBUG_PEDANTIC",
    ]
    print("compile command:", " ".join(cmd_))
    print("compile... ", end="")
    subprocess.run(cmd_)
    print("done")
    return


def judge(input_file_path: str, outputs: list) -> int:
    first_line = outputs[0].strip().split(" ")
    second_line = outputs[1].strip().split(" ")
    # m = int(first_line[0])
    r = list(map(int, first_line[1:]))
    n = int(second_line[0])
    buffer = second_line[1:]
    xy = list()
    for i in range(0, 2 * n - 1, 2):
        xy.append((int(buffer[i]), int(buffer[i + 1])))
    if xy[0] != (400, 400) or xy[-1] != (400, 400):
        print("[error]")
        print(f"start/goal should be (400,400). but start:{xy[0]} / goal:{xy[-1]}")
        print("outputs:", "\n".join(outputs))
        return -1
    coord = set()
    from_ = dict()
    to_ = dict()
    with open(input_file_path, encoding="UTF-8") as fp:
        for id, line in enumerate(fp.readlines(), start=1):
            a, b, c, d = map(int, line.split(" "))
            coord.add((a, b))
            coord.add((c, d))
            from_[id] = (a, b)
            to_[id] = (c, d)
    obj = 0
    # (x,y)に何番目に訪れたか
    first_visit = dict()
    last_visit = dict()
    for i in range(len(xy)):
        px, py = xy[i]
        qx, qy = xy[(i + 1) % len(xy)]
        obj += (abs(px - qx) + abs(py - qy))
        if xy[i] not in first_visit:
            first_visit[xy[i]] = i
        last_visit[xy[i]] = i
        if px < 0 or 800 < px or py < 0 or 800 < py:
            print("[error]")
            print(f"coord must be 0<=x,y<=800. but xy[{i}]:{xy[i]}")
            print("outputs:", "\n".join(outputs))
            return -1
    for i in r:
        if last_visit[to_[i]] < first_visit[from_[i]]:
            print("[error]")
            print(f"{i}-th delivery has not been completed")
            print(f"reached 'to'@{last_visit[to_[i]]} before 'from'@{first_visit[from_[i]]}")
            print("outputs:", "\n".join(outputs))
            return -1
    score = 10**8 // (1000 + obj)
    print(f"-> obj:{obj} / score:{score}")
    return score


def main(args: argparse.Namespace) -> None:
    src = os.path.join("solution", args.src)
    print("src:", src)
    base, ext = args.src.split(".")
    base = os.path.join("solution", base)
    if ext in ["cpp", "cc"]:
        if not args.no_compile:
            compile(args.src)
        cmd_ = [f"{base}.bin"]
    else:
        cmd_ = ["python", src]
    scores = list()
    invalid_cases = list()
    IN = os.path.join("given_tools", "in")
    for i in range(*eval(args.range)):
        test_case = str(i).zfill(4)
        print(f"\n#case:{test_case}")
        input_file_path = os.path.join(IN, f"{test_case}.txt")
        with open(input_file_path, encoding="UTF-8") as fp:
            res = subprocess.run(
                cmd_,
                encoding="UTF-8",
                stdin=fp,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        outputs = res.stdout.strip().split("\n")[-2:]  # 最後の出力だけを抜き出し
        score = judge(input_file_path, outputs)
        if score == -1:
            invalid_cases.append(test_case)
            scores.append(0)
            break  # DEBUG:
        else:
            scores.append(score)

    print("\n\n[result]")
    print(f"total score: {sum(scores):,} / {10**5 * len(scores):,}")
    print(f"max score: {max(scores):,} at {scores.index(max(scores))}")
    print(f"min score: {min(scores):,} at {scores.index(min(scores))}")
    print(f"invalid cases:{invalid_cases}")
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser("LOCAL Judge system")
    parser.add_argument(
        "--src",
        type=str,
        required=True,
        help="target script name (include file ext. e.g. 'main.py' )"
    )
    parser.add_argument(
        "--no_compile",
        action="store_true",
        help="this is a flag whether or not to compile."
    )
    parser.add_argument(
        "--range",
        default="0,100",
        type=str,
        help="set range of the dataset that you want to check (e.g. '0,100' means 0,1,...,99 )"
    )
    args = parser.parse_args()
    main(args)
