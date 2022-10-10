"""
1年前にリアルタイム参加した時に作ったもの
無知なりに頑張っている
"""


from heapq import heappop, heappush


# input
ab = [list() for _ in range(1000)]
cd = [list() for _ in range(1000)]
ManhattanDist = [tuple for _ in range(1000)]  # (dist, idx)
for idx in range(1000):
    buffer = list(map(int, input().split()))
    ab[idx] = buffer[:2]
    cd[idx] = buffer[2:]
    ManhattanDist[idx] = (abs(ab[idx][0] - 400) + abs(ab[idx][1] - 400) +
                          abs(cd[idx][0] - 400) + abs(cd[idx][1] - 400), idx)


# process
ManhattanDist.sort()
pq = list()  # priority queue := (priority, idx)
R = list()  # ans
xy = list()  # ans


def createPath(mode: str):
    while pq:
        idx = heappop(pq)[1]
        if mode == "get":
            xy.append(ab[idx])
        elif mode == "post":
            xy.append(cd[idx])
    return


for _, idx in ManhattanDist[:50]:
    R.append(idx + 1)
    x, y = ab[idx]
    if x >= y:
        heappush(pq, (x + y, idx))
    else:
        heappush(pq, (10000 + (800 - x + 800 - y), idx))
xy.append([400, 400])  # start
createPath("get")

for _, idx in ManhattanDist[:50]:
    x, y = cd[idx]
    if x >= y:
        heappush(pq, (-(x + y), idx))
    else:
        heappush(pq, (-(10000 + (800 - x + 800 - y)), idx))
createPath("post")
xy.append([400, 400])  # goal


# output
print(len(R), *R)
print(len(xy), end=" ")
for out in xy:
    print(*out, end=" ")
print("")
