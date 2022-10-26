import argparse
import os
import subprocess


def compile(src: str) -> None:
    # c++fileのコンパイル
    print("compile... ", end="")
    base, ext = src.split(".")
    assert ext in ["cpp", "cc"], (src, ext)
    base = os.path.join("solution", base)
    src = os.path.join("solution", src)
    subprocess.run(
        [
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
    )
    print("done")
    return


def judge(input_file_path: str, outputs: list) -> int:
    first_line = outputs[0].strip().split(" ")
    second_line = outputs[1].strip().split(" ")
    m = int(first_line[0])
    # r = list(map(int, first_line[1:]))
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
    from_ = list()
    to_ = list()
    with open(input_file_path, encoding="UTF-8") as fp:
        for line in fp.readlines():
            a, b, c, d = map(int, line.split(" "))
            from_.append((a, b))
            to_.append((c, d))
    # FIXME: objの計算バグってる
    obj = 0
    pending = set()  # fromで受け取ったid集合
    assert len(xy) == n, (xy, len(xy), n)
    for i in range(n):
        px, py = xy[i]
        qx, qy = xy[(i + 1) % m]
        obj += (abs(px - qx) + abs(py - qy))
        # check
        if 0 < i < n - 1:
            # not start/goal
            if xy[i] in from_:
                pending.add(from_.index(xy[i]))
                continue
            else:
                # FIXME: ここの判定もバグってるっぽい
                if to_.index(xy[i]) not in pending:
                    print("[error]")
                    print("reached 'to' before 'from'")
                    print(f"id:{to_.index(xy[i])} / from:{from_[to_.index(xy[i])]} to:{xy[i]} / i:{i}")
                    print("outputs:", "\n".join(outputs))
                    return -1
    score = 10**8 // (1000 + obj)
    print(f"-> obj:{obj} / score:{score}")
    return score


def main(args) -> None:
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
    for i in range(20):
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
        else:
            scores.append(score)

    print("\n\n[result]")
    print(f"total score: {sum(scores):,} / 10,000,000")
    print(f"max score: {max(scores):,} at {scores.index(max(scores))}")
    print(f"min score: {min(scores):,} at {scores.index(min(scores))}")
    print(f"invalid cases:{invalid_cases}")
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser("LOCAL Judge system")
    parser.add_argument(
        "--src",
        default="main.cpp",
        type=str,
        help="target script name (include file ext. e.g. main.py)"
    )
    parser.add_argument(
        "--no_compile",
        action="store_true",
        help="this is a flag whether or not to compile."
    )
    args = parser.parse_args()
    main(args)
