import argparse
import os
import subprocess
import matplotlib.pyplot as plt


def main(script: str) -> None:
    SCRIPT = os.path.join("solution", script)
    print("script:", SCRIPT)
    if SCRIPT[-2:] == "py":
        cmd__ = ["python", SCRIPT]
    else:
        cmd__ = [SCRIPT]
    scores = list()
    IN = os.path.join("given_tools", "in")
    for i in range(100):
        input = str(i).zfill(4)
        INPUT = os.path.join(IN, f"{input}.txt")
        # print("input:", INPUT)
        with open(INPUT, encoding="UTF-8") as fp:
            cmd_ = subprocess.run(
                cmd__,
                encoding="UTF-8",
                stdin=fp,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        # print("stdout:", cmd_.stdout)
        print(f"case:{input}.txt -> score:{cmd_.stderr}")
        scores.append(float(cmd_.stderr.strip()))
    print("\n[test result]")
    print("total score:", sum(scores))
    print(f"max score:{max(scores)} at {scores.index(max(scores))}")
    print(f"min score:{min(scores)} at {scores.index(min(scores))}")
    plt.scatter(range(len(scores)), scores)
    plt.ylim(0, 20000)
    plt.xticks(range(0, len(scores) + 1, 10))
    plt.yticks(range(0, 20000 + 1, 1000))
    plt.xlabel("test case")
    plt.ylabel("score")
    plt.title(f"score distribution\nscript:{script}")
    plt.text(
        # (x,y,string)
        0,
        100,
        f"\
        total score:{sum(scores):.3f}\n\
        max score:{max(scores):.3f} at {scores.index(max(scores))}\n\
        low score:{min(scores):.3f} at {scores.index(min(scores))}"
    )
    plt.show()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--src",
        default="main.bin",
        type=str,
        help="target script name (include file ext. e.g. main.py)"
    )
    args = parser.parse_args()
    main(args.src)
