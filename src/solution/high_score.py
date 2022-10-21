"""
ハイスコアを取れたプログラム
結果の比較用
"""

import sys
import time
import random


class Data:
    def __init__(self) -> None:
        # office coordinate: 本社の座標/start/goal -> tuple(x, y)
        self.OFFICE_COORD: tuple = (400, 400)
        # threshold: 閾値
        self.THRESHOLD = 200
        self.ALL_ORDER_NUM: int = 1000
        self.ORDER_NUM: int = 50
        # tour length = order_num * 2(from/to) + 2(start/end)
        self.TOUR_LEN: int = self.ORDER_NUM * 2 + 2
        # _from[id] := id番目の注文の受取先 -> tuple(x,y)
        self._from: list = [tuple() for _ in range(self.ALL_ORDER_NUM)]
        # _to[id] := id番目の注文の配達先 -> tuple(x,y)
        self._to: list = [tuple() for _ in range(self.ALL_ORDER_NUM)]
        # or optのsub pathのmax length
        self.MIN_OR_OPT_LEN: int = 1
        self.MAX_OR_OPT_LEN: int = 3
        return

    def read_data(self) -> None:
        for id in range(self.ALL_ORDER_NUM):
            a, b, c, d = map(int, input().split())
            self._from[id] = (a, b)
            self._to[id] = (c, d)
        return


class Job:
    """ 解の単位 """

    def __init__(self, data: Data) -> None:
        """
        - idとi,j,k...の使い分け
            - idは注文のid(1000件の内何番目か)
            - i,j,kは経路内のindex
        """
        # selected orders id: 受ける注文のid列(1-indexed)
        self.orders_id: list = [None for _ in range(data.ORDER_NUM)]
        # tour: 移動経路 -> tuple(order_id, order_destination:from/to). これが解
        self.tour: list = [(-1, None) for _ in range(data.TOUR_LEN)]
        # objective value: total tour distance 配達ルートの総移動距離 -> minimize
        self.obj: int = float("inf")
        return


class Solver:
    """ インターフェース的なクラス """

    def __init__(self) -> None:
        self.TIME_LIMIT: float = 1.95  # [sec]
        self.START_TIME: float = time.time()
        self.data: Data = Data()  # 問題で与えられた情報
        self.data.read_data()  # 入力の受け取り
        self.job: Job = Job(self.data)  # 解
        return

    def check_time_limit(self, msg: str = None) -> bool:
        return (time.time() - self.START_TIME) < self.TIME_LIMIT

    def eval_dist(self, p: tuple, q: tuple) -> int:
        """ p(x,y)とq(x,y)のマンハッタン距離を返す
        ### 注意：ユークリッド距離では無い
        """
        return (abs(p[0] - q[0]) + abs(p[1] - q[1]))

    def calc_obj(self, job: Job) -> None:
        """jobの目的値(objective value):総移動距離を計算．
        """
        job.obj = 0
        for i in range(self.data.TOUR_LEN):
            job.obj += self.eval_dist(
                self.get_coord_by_id(*self.job.tour[i]),
                self.get_coord_by_id(*self.job.tour[(i + 1) % self.data.TOUR_LEN])
            )
        return

    def get_current_score(self) -> float:
        self.calc_obj(self.job)
        return (10**8 / (1000 + self.job.obj))

    def get_coord_by_id(self, id: int, dest: str) -> tuple:
        """ order_idとorder_destからcoord(座標)を取得
        """
        if id == -1:
            return self.data.OFFICE_COORD
        if dest == "from":
            return self.data._from[id]
        elif dest == "to":
            return self.data._to[id]

    def select_order(self) -> None:
        """受ける注文を選択(1000 -> 50)
        ### 方針
        [中心(office)からfrom,toへのマンハッタン距離の和+ペナルティ]が小さい上位50件を採用する．
        遠すぎるものにペナルティを課すことで，(fromが近くてtoが遠い)といったものを排除している．
        """
        # manhattan_dist[i]=tuple(dist,id)
        manhattan_dist = [tuple() for _ in range(self.data.ALL_ORDER_NUM)]
        for id in range(self.data.ALL_ORDER_NUM):
            from_dist = self.eval_dist(self.data.OFFICE_COORD, self.data._from[id])
            to_dist = self.eval_dist(self.data.OFFICE_COORD, self.data._to[id])
            # 閾値よりも遠いものにはペナルティを課す
            if self.data.THRESHOLD < from_dist:
                from_dist *= 2
            if self.data.THRESHOLD < to_dist:
                to_dist *= 2
            manhattan_dist[id] = (from_dist + to_dist, id)
        manhattan_dist.sort()
        self.job.orders_id = [id for _, id in manhattan_dist[:self.data.ORDER_NUM]]
        return

    def nearest_neighbor(self) -> None:
        """最近近傍法(貪欲法) - 初期解の生成．fromを全部回ってからtoを全部回る
        - TODO: GRASP(近いものの上位n件から確率で選ぶ)にしても良いかも
        """
        # from
        done = set()  # 確定済み
        for i in range(self.data.ORDER_NUM):  # fromの全てのiについて
            min_dist, min_dist_j = float("inf"), None
            for j in range(self.data.ORDER_NUM):  # i+1の候補をjで全探索
                if j in done:
                    continue
                cur_dist = self.eval_dist(
                    self.get_coord_by_id(*self.job.tour[i]),
                    self.get_coord_by_id(self.job.orders_id[j], "from")
                )
                if cur_dist < min_dist:
                    min_dist = cur_dist
                    min_dist_j = j
            if min_dist_j is not None:
                self.job.tour[i + 1] = (self.job.orders_id[min_dist_j], "from")
                done.add(min_dist_j)
            # self.print()
        assert len(done) == self.data.ORDER_NUM,\
            (len(done), self.data.ORDER_NUM, set(range(self.data.ORDER_NUM)) - done)
        # to
        done = set()  # 確定済み
        for i in range(self.data.ORDER_NUM):  # toの全てのiについて
            min_dist, min_dist_j = float("inf"), None
            for j in range(self.data.ORDER_NUM):  # i+1の候補をjで全探索
                if j in done:
                    continue
                cur_dist = self.eval_dist(
                    self.get_coord_by_id(*self.job.tour[self.data.ORDER_NUM + i]),
                    self.get_coord_by_id(self.job.orders_id[j], "to")
                )
                if cur_dist < min_dist:
                    min_dist = cur_dist
                    min_dist_j = j
            if min_dist_j is not None:
                self.job.tour[self.data.ORDER_NUM + i + 1] = (self.job.orders_id[min_dist_j], "to")
                done.add(min_dist_j)
            # self.print()
        assert len(done) == self.data.ORDER_NUM,\
            (len(done), self.data.ORDER_NUM, set(range(self.data.ORDER_NUM)) - done)
        return

    def local_search(self) -> None:
        """ 局所探索法
        """
        self.print()
        while self.check_time_limit():
            self.two_opt_search()
            self.or_opt_search()
        return

    def two_opt_search(self) -> None:
        """ 2-opt
        - とりあえず，from同士，to同士で考える
        """

        def eval_diff(i: int, j: int) -> int:
            """
            i -> i+1 -> ... -> j -> j+1 を
            i -> j -> ... -> i+1 -> j+1 に変えたときのobjの差分を返す
            """
            cur_dist = self.eval_dist(
                self.get_coord_by_id(*self.job.tour[i]),
                self.get_coord_by_id(*self.job.tour[i + 1])
            ) + self.eval_dist(
                self.get_coord_by_id(*self.job.tour[j]),
                self.get_coord_by_id(*self.job.tour[j + 1])
            )
            new_dist = self.eval_dist(
                self.get_coord_by_id(*self.job.tour[i]),
                self.get_coord_by_id(*self.job.tour[j])
            ) + self.eval_dist(
                self.get_coord_by_id(*self.job.tour[i + 1]),
                self.get_coord_by_id(*self.job.tour[j + 1])
            )
            return new_dist - cur_dist

        def update_tour(i: int, j: int) -> None:
            """
            i -> i+1 -> ... -> j -> j+1 を
            i -> j -> ... -> i+1 -> j+1 に変える
            """
            # 半開区間だから分かりずらいけど，区間[i+1,j]をreverseしている
            self.job.tour[i + 1:j + 1] = list(reversed(self.job.tour[i + 1:j + 1]))
            return

        from_nbhd = (  # fromだけ
            (i, j)
            for i in range(1, self.data.ORDER_NUM - 1)
            for j in range(i + 2, self.data.ORDER_NUM - 1)
        )
        to_nbhd = (  # toだけ
            (i, j)
            for i in range(self.data.ORDER_NUM, self.data.TOUR_LEN - 2)
            for j in range(i + 2, self.data.TOUR_LEN - 2)
        )
        nbhd = list(from_nbhd) + list(to_nbhd)  # len(nbhd)=2304
        for i, j in nbhd:
            diff = eval_diff(i, j)
            if diff < 0:
                # 改善するなら即時変更する
                self.job.obj += diff  # decrease obj -> increase score
                update_tour(i, j)
                assert self.job.tour[0] == (-1, None), (i, j)
                assert self.job.tour[-1] == (-1, None), (i, j)
                # break
            if not self.check_time_limit("2-opt"):
                # タイムリミット
                return
        return

    def or_opt_search(self) -> None:
        """ or-opt
        - とりあえず，from同士，to同士で考える
        """

        def eval_diff(l: int, i: int, j: int) -> tuple:
            """
            i-1 -> i -> ... -> i+l-1 -> i+l -> ... -> j -> j+1 を
            i-1 -> i+l -> ... -> j -> i -> ... -> i+l-1 -> j+1 に変えたときのobjの差分を返す
            """
            prev_p = self.job.tour[i - 1]
            head_p = self.job.tour[i]
            tail_p = self.job.tour[i + l - 1]
            next_p = self.job.tour[i + l]
            q = self.job.tour[j]
            next_q = self.job.tour[j + 1]
            cur_dist = self.eval_dist(
                self.get_coord_by_id(*prev_p),
                self.get_coord_by_id(*head_p)
            ) + self.eval_dist(
                self.get_coord_by_id(*tail_p),
                self.get_coord_by_id(*next_p)
            ) + self.eval_dist(
                self.get_coord_by_id(*q),
                self.get_coord_by_id(*next_q)
            )
            new_dist_fwd = self.eval_dist(
                self.get_coord_by_id(*prev_p),
                self.get_coord_by_id(*next_p)
            ) + self.eval_dist(
                self.get_coord_by_id(*q),
                self.get_coord_by_id(*head_p)
            ) + self.eval_dist(
                self.get_coord_by_id(*tail_p),
                self.get_coord_by_id(*next_q)
            )
            new_dist_bwd = self.eval_dist(
                self.get_coord_by_id(*prev_p),
                self.get_coord_by_id(*next_p)
            ) + self.eval_dist(
                self.get_coord_by_id(*q),
                self.get_coord_by_id(*tail_p)
            ) + self.eval_dist(
                self.get_coord_by_id(*head_p),
                self.get_coord_by_id(*next_q)
            )
            if new_dist_fwd <= new_dist_bwd:
                return (new_dist_fwd - cur_dist, "fwd")
            else:
                return (new_dist_bwd - cur_dist, "bwd")

        def update_tour(fb: str, l: int, i: int, j: int) -> None:
            """
            i-1 -> i -> ... -> i+l-1 -> i+l -> ... -> j -> j+1 を
            i-1 -> i+l -> ... -> j -> i -> ... -> i+l-1 -> j+1 に変えたときのobjの差分を返す
            """
            assert 1 <= i <= self.data.TOUR_LEN - l - 1, (1, i, self.data.TOUR_LEN - l - 1)
            assert i + l <= j <= self.data.TOUR_LEN - 1, (i + l, j, self.data.TOUR_LEN - 1)
            # temporary save
            if fb == "fwd":
                subpath = self.job.tour[i:i + l]
            elif fb == "bwd":
                subpath = list(reversed(self.job.tour[i:i + l]))
            # shift
            for k in range(i + l, j + 1):
                self.job.tour[k - l] = self.job.tour[k]
            # shift
            self.job.tour[j - l + 1:j + 1] = subpath[:]
            return

        from_nbhd = (  # fromだけ
            (l, i, j)
            for l in range(self.data.MIN_OR_OPT_LEN, self.data.MAX_OR_OPT_LEN + 1)
            for i in range(1 + l, self.data.ORDER_NUM - l - 1)
            for j in range(i + l, self.data.ORDER_NUM - 1)
        )
        to_nbhd = (  # toだけ
            (l, i, j)
            for l in range(self.data.MIN_OR_OPT_LEN, self.data.MAX_OR_OPT_LEN + 1)
            for i in range(1 + self.data.ORDER_NUM, self.data.TOUR_LEN - l - 2)
            for j in range(i + l, self.data.TOUR_LEN - 2)
        )
        nbhd = list(from_nbhd) + list(to_nbhd)  # OR_LEN=3のとき，len(nbhd)=6912
        for l, i, j in nbhd:
            diff, fb = eval_diff(l, i, j)
            if diff < 0:
                # 改善するなら即時変更する
                self.job.obj += diff  # decrease obj -> increase score
                update_tour(fb, l, i, j)
                assert self.job.tour[0] == (-1, None), (fb, l, i, j)
                assert self.job.tour[-1] == (-1, None), (fb, l, i, j)
                return
            if not self.check_time_limit("or-opt"):
                # タイムリミット
                return
        return

    def print(self) -> None:
        """ 答えの出力
        """
        r = [i + 1 for i in self.job.orders_id]  # MEMO: 1-indexedに変換
        print(len(r), *r)
        tour_coord = [self.get_coord_by_id(id, dest) for id, dest in self.job.tour]
        tour_coord = [f"{x} {y}" for x, y in tour_coord]
        print(len(tour_coord), " ".join(tour_coord))
        return


def main() -> None:
    random.seed(0)  # debugのために固定する
    solver = Solver()
    solver.select_order()
    solver.nearest_neighbor()
    solver.local_search()
    solver.print()
    print(solver.get_current_score(), file=sys.stderr)


if __name__ == '__main__':
    args = sys.argv
    DEBUG = bool("debug" in args)
    main()
