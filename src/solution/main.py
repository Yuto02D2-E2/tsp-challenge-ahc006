import sys
import time
import random
import dataclasses

# k-means test
import numpy as np
from sklearn.cluster import KMeans
# import seaborn as sns
# import matplotlib.pyplot as plt


class Data:
    def __init__(self) -> None:
        # map size (width=height)
        self.MAP_SIZE: int = 800
        # office coordinate: 本社の座標/start/goal -> tuple(x, y)
        self.OFFICE_COORD: tuple = (400, 400)
        # threshold: 閾値
        self.ALL_ORDER_NUM: int = 1000
        self.ORDER_NUM: int = 50
        # tour length = order_num * 2(from/to) + 2(start/end)
        self.TOUR_LEN: int = self.ORDER_NUM * 2 + 2
        # _from[id] := id番目の注文の受取先 -> tuple(x,y)
        self._from: list = [tuple() for _ in range(self.ALL_ORDER_NUM)]
        # _to[id] := id番目の注文の配達先 -> tuple(x,y)
        self._to: list = [tuple() for _ in range(self.ALL_ORDER_NUM)]
        # or optのsub pathのlength
        self.MIN_OR_OPT_LEN: int = 1
        self.MAX_OR_OPT_LEN: int = 3
        return

    def read_data(self) -> None:
        for id in range(self.ALL_ORDER_NUM):
            a, b, c, d = map(int, input().split())
            self._from[id] = (a, b)
            self._to[id] = (c, d)
        return


@dataclasses.dataclass
class Order:
    id: int = -1
    dest: str = "office"


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
        # tour: 移動経路 -> Order(order_id, order_destination:from/to). これが解
        self.tour: list = [Order for _ in range(data.TOUR_LEN)]
        # objective value: total tour distance 配達ルートの総移動距離 -> minimize
        self.obj: int = float("inf")
        return


class Solver:
    """ インターフェース的なクラス """

    def __init__(self) -> None:
        # self.TIME_LIMIT: float = 30.0  # [sec]
        self.TIME_LIMIT: float = 1.95  # [sec]
        self.START_TIME: float = time.time()
        self.data: Data = Data()  # 問題で与えられた情報
        self.data.read_data()  # 入力の受け取り
        self.job: Job = Job(self.data)  # 最適解
        return

    def check_time_limit(self, msg: str = None) -> bool:
        if DEBUG and msg is not None:
            print(f"current time:{time.time() - self.START_TIME:.10f}[sec], msg:{msg}")
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
                self.get_coord_by_order(self.job.tour[i]),
                self.get_coord_by_order(self.job.tour[(i + 1) % self.data.TOUR_LEN])
            )
        return

    def get_current_score(self) -> float:
        self.calc_obj(self.job)
        return (10**8 / (1000 + self.job.obj))

    def get_coord_by_order(self, order: Order) -> tuple:
        """ orderからcoord(座標)を取得
        """
        if order.dest == "from":
            return self.data._from[order.id]
        elif order.dest == "to":
            return self.data._to[order.id]
        else:
            return self.data.OFFICE_COORD

    def select_order(self) -> None:
        """受ける注文を選択(1000 -> 50)
        ### 方針
        [中心(office)からfrom,toへのマンハッタン距離の和+ペナルティ]が小さい上位50件を採用する．
        遠すぎるものにペナルティを課すことで，(fromが近くてtoが遠い)といったものを排除している．
        """
        threshold = 200
        # manhattan_dist[i]=tuple(dist,id)
        manhattan_dist = [tuple() for _ in range(self.data.ALL_ORDER_NUM)]
        for id in range(self.data.ALL_ORDER_NUM):
            from_dist = self.eval_dist(self.data.OFFICE_COORD, self.data._from[id])
            to_dist = self.eval_dist(self.data.OFFICE_COORD, self.data._to[id])
            # 閾値よりも遠いものにはペナルティを課す
            if threshold < from_dist:
                from_dist *= 5
            if threshold < to_dist:
                to_dist *= 5
            manhattan_dist[id] = (from_dist + to_dist, id)
        manhattan_dist.sort()
        self.job.orders_id = [id for _, id in manhattan_dist[:self.data.ORDER_NUM]]
        return

    def k_means(self) -> tuple:
        # from/toそれぞれをk=8くらいのk-meansで分類したい
        # ただ，atcoderではscikit-learnは使えなさそうなので，自力で実装しないとダメそう．C++使う？
        # 速度的にはC++じゃないとダメそうだけど，実装かなりしんどいので，
        # 有用性を確かめるために一旦scikit-learn使ってPythonで解いてビジュアライズしてみるかな
        # → 良さそう！！！C++で書いてみるか．
        K = 8
        init_centers = np.array([
            [200, 200],
            [200, 400],
            [200, 600],
            [400, 200],
            [400, 600],
            [600, 200],
            [600, 400],
            [600, 600]
        ])
        kmeans = KMeans(
            n_clusters=K,
            max_iter=1000,
            random_state=1,
            init=init_centers,
            n_init=1,
            tol=0.001,
        )
        # from
        from_coords = np.array([self.get_coord_by_order(Order(id, "from")) for id in self.job.orders_id])
        kmeans.fit(from_coords)
        from_labels = kmeans.labels_
        from_cluster = [set() for _ in range(K)]
        for i, id in enumerate(self.job.orders_id):
            from_cluster[from_labels[i]].add(id)
        # to
        to_coords = np.array([self.get_coord_by_order(Order(id, "to")) for id in self.job.orders_id])
        kmeans.fit(to_coords)
        to_labels = kmeans.labels_
        to_cluster = [set() for _ in range(K)]
        for i, id in enumerate(self.job.orders_id):
            to_cluster[to_labels[i]].add(id)
        return (from_cluster, to_cluster)

    def nearest_neighbor(self) -> None:
        """最近近傍法(貪欲法) - 初期解の生成．fromを全部回ってからtoを全部回る
        - TODO: GRASP(近いものの上位n件から確率で選ぶ)にしても良いかも
        - TODO: k-meansで分類して，各集合内でnearest_neightborをする．最後にそれぞれを時計回りで繋いで完成
        """
        from_cluster, to_cluster = self.k_means()
        tour_i = 1
        # from
        done = set()
        # from_quadrants = self.divide_to_quadrants(self.data._from)
        for one_quadrant in from_cluster:
            for i_id in one_quadrant:  # fromの全てのiについて
                min_dist, min_dist_j_id = float("inf"), None
                for j_id in one_quadrant:  # i+1の候補をjで全探索
                    if j_id in done:
                        # if i_id == j_id:
                        continue
                    cur_dist = self.eval_dist(
                        self.data._from[i_id],
                        self.data._from[j_id]
                    )
                    if cur_dist < min_dist:
                        min_dist = cur_dist
                        min_dist_j_id = j_id
                if min_dist_j_id is not None:
                    self.job.tour[tour_i] = Order(min_dist_j_id, "from")
                    tour_i += 1
                    done.add(min_dist_j_id)
                self.print()
        assert tour_i == self.data.ORDER_NUM + 1, (tour_i, self.data.ORDER_NUM)
        # to
        done = set()
        # to_quadrants = self.divide_to_quadrants(self.data._to)
        for one_quadrant in reversed(to_cluster):
            for i_id in one_quadrant:  # toの全てのiについて
                min_dist, min_dist_j_id = float("inf"), None
                for j_id in one_quadrant:  # i+1の候補をjで全探索
                    if j_id in done:
                        # if i_id == j_id:
                        continue
                    cur_dist = self.eval_dist(
                        self.data._to[i_id],
                        self.data._to[j_id],
                    )
                    if cur_dist < min_dist:
                        min_dist = cur_dist
                        min_dist_j_id = j_id
                if min_dist_j_id is not None:
                    self.job.tour[tour_i] = Order(min_dist_j_id, "to")
                    tour_i += 1
                    done.add(min_dist_j_id)
                self.print()
        assert tour_i == self.data.TOUR_LEN - 1, (tour_i, self.data.TOUR_LEN - 1)
        return

    def local_search(self) -> None:
        """ 局所探索法
        """
        self.print()
        return
        while self.check_time_limit():
            self.two_opt_search()
            self.or_opt_search()
            # self.update_order()
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
                self.get_coord_by_order(self.job.tour[i]),
                self.get_coord_by_order(self.job.tour[i + 1])
            ) + self.eval_dist(
                self.get_coord_by_order(self.job.tour[j]),
                self.get_coord_by_order(self.job.tour[j + 1])
            )
            new_dist = self.eval_dist(
                self.get_coord_by_order(self.job.tour[i]),
                self.get_coord_by_order(self.job.tour[j])
            ) + self.eval_dist(
                self.get_coord_by_order(self.job.tour[i + 1]),
                self.get_coord_by_order(self.job.tour[j + 1])
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
                self.print()
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
                self.get_coord_by_order(prev_p),
                self.get_coord_by_order(head_p)
            ) + self.eval_dist(
                self.get_coord_by_order(tail_p),
                self.get_coord_by_order(next_p)
            ) + self.eval_dist(
                self.get_coord_by_order(q),
                self.get_coord_by_order(next_q)
            )
            new_dist_fwd = self.eval_dist(
                self.get_coord_by_order(prev_p),
                self.get_coord_by_order(next_p)
            ) + self.eval_dist(
                self.get_coord_by_order(q),
                self.get_coord_by_order(head_p)
            ) + self.eval_dist(
                self.get_coord_by_order(tail_p),
                self.get_coord_by_order(next_q)
            )
            new_dist_bwd = self.eval_dist(
                self.get_coord_by_order(prev_p),
                self.get_coord_by_order(next_p)
            ) + self.eval_dist(
                self.get_coord_by_order(q),
                self.get_coord_by_order(tail_p)
            ) + self.eval_dist(
                self.get_coord_by_order(head_p),
                self.get_coord_by_order(next_q)
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
                self.print()
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
        if DEBUG:
            print("[answer]")
            print("format: m r_1 ... r_m \t(1-indexed)")
        r = [i + 1 for i in self.job.orders_id]  # MEMO: 1-indexedに変換
        print(len(r), *r)
        tour_coord = [self.get_coord_by_order(order) for order in self.job.tour]
        if DEBUG:
            print("format: n x_1 y_1 ... x_y y_n")
            print(len(tour_coord), *tour_coord)
            return
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
    if DEBUG:
        print("final score:", solver.get_current_score())
        print("final objective value:", solver.job.obj)
        solver.check_time_limit("end")
    print(solver.get_current_score(), file=sys.stderr)


if __name__ == '__main__':
    args = sys.argv
    DEBUG = bool("debug" in args)
    main()
