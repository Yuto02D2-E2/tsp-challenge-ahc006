import sys
import time
import random


class Data:
    """ 問題文で与えられたデータ """

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
        # pos[id][dest]=i -> id,destからtour内のpos(=i)を復元. tourの逆
        self.pos: dict = dict()
        # objective value: total tour distance 配達ルートの総移動距離 -> minimize
        self.obj: int = float("inf")
        return


class Solver:
    """ インターフェース的なクラス """

    def __init__(self) -> None:
        self.TIME_LIMIT: float = 1.8  # [sec]
        self.START_TIME: float = time.time()
        self.data: Data = Data()  # 問題で与えられた情報
        self.data.read_data()  # 入力の受け取り
        self.best_job: Job = Job(self.data)  # 最適解
        return

    def check_time_limit(self) -> bool:
        if DEBUG:
            print(f"current time:{time.time() - self.START_TIME}[sec]")
        return (time.time() - self.START_TIME) < self.TIME_LIMIT

    def calc_dist(self, p: tuple, q: tuple) -> int:
        """ p(x,y)とq(x,y)のマンハッタン距離を返す
        """
        return (abs(p[0] - q[0]) + abs(p[1] - q[1]))

    def calc_obj(self, job: Job) -> None:
        """jobの目的値(objective value):総移動距離を計算．
        """
        job.obj = 0
        for i in range(self.data.TOUR_LEN):
            job.obj += self.calc_dist(
                self.get_coord_by_id(*self.best_job.tour[i]),
                self.get_coord_by_id(*self.best_job.tour[(i + 1) % self.data.TOUR_LEN])
            )
        return job.obj

    def get_current_score(self) -> float:
        self.calc_obj(self.best_job)
        return (10**8 / (1000 + self.best_job.obj))

    def get_coord_by_id(self, id: int, dest: str) -> tuple:
        """ order_idとorder_destからcoordを取得
        """
        if id == -1:
            return self.data.OFFICE_COORD
        if dest == "from":
            return self.data._from[id]
        elif dest == "to":
            return self.data._to[id]

    def select_order(self) -> None:
        """受ける注文を選択(1000 -> 50)
        - 方針：中心(office)からfrom,toへのマンハッタン距離の和が小さい上位50件を採用する
        """
        # manhattan_dist[i]=tuple(dist,id)
        manhattan_dist = [tuple() for _ in range(self.data.ALL_ORDER_NUM)]
        for id in range(self.data.ALL_ORDER_NUM):
            from_dist = self.calc_dist(self.data.OFFICE_COORD, self.data._from[id])
            to_dist = self.calc_dist(self.data.OFFICE_COORD, self.data._to[id])
            # 余りにも遠いものにはペナルティを課す
            if self.data.THRESHOLD < from_dist:
                from_dist *= 2
            if self.data.THRESHOLD < to_dist:
                to_dist *= 2
            manhattan_dist[id] = (from_dist + to_dist, id)
        manhattan_dist.sort()
        self.best_job.orders_id = [id for _, id in manhattan_dist[:self.data.ORDER_NUM]]
        return

    def nearest_neighbor(self) -> None:
        """最近近傍法(貪欲法) - 初期解の生成
        """
        # from
        done = set()  # 確定済み
        for i in range(self.data.ORDER_NUM):  # fromの全てのiについて
            min_dist, min_dist_j = float("inf"), None
            for j in range(self.data.ORDER_NUM):  # i+1の候補をjで全探索
                if j in done:
                    continue
                cur_dist = self.calc_dist(
                    self.get_coord_by_id(*self.best_job.tour[i]),
                    self.get_coord_by_id(self.best_job.orders_id[j], "from")
                )
                if cur_dist < min_dist:
                    min_dist = cur_dist
                    min_dist_j = j
            if min_dist_j is not None:
                self.best_job.tour[i + 1] = (self.best_job.orders_id[min_dist_j], "from")
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
                cur_dist = self.calc_dist(
                    self.get_coord_by_id(*self.best_job.tour[self.data.ORDER_NUM + i]),
                    self.get_coord_by_id(self.best_job.orders_id[j], "to")
                )
                if cur_dist < min_dist:
                    min_dist = cur_dist
                    min_dist_j = j
            if min_dist_j is not None:
                self.best_job.tour[self.data.ORDER_NUM + i + 1] = (self.best_job.orders_id[min_dist_j], "to")
                done.add(min_dist_j)
            # self.print()
        assert len(done) == self.data.ORDER_NUM,\
            (len(done), self.data.ORDER_NUM, set(range(self.data.ORDER_NUM)) - done)
        return

    def local_search(self) -> None:
        """ 局所探索法
        - 方針：初期解を元に2-opt()とかupdate_order()を使って改善していく
        """
        # TODO:
        if DEBUG:
            print("solver.ls() is not supported yet :(")
        return

    def two_opt_search(self, cur_tour: Job) -> None:
        """ 2-opt探索
        """
        # TODO:
        if DEBUG:
            print("solver.2-opt() is not supported yet :(")
        return

    def update_order(self) -> None:
        # TODO:
        # randomでmanhattan distが小さくてまだ使われていない( and 確認したことが無い)ものとswapする？
        # tabuリスト的なかんじ？
        # random要素を入れてもいいかも
        return

    def print(self) -> None:
        """ 答えの出力
        """
        # TODO:
        if DEBUG:
            print("[answer]")
            print("format: m r_1 ... r_m \t(1-indexed)")
        r = [i + 1 for i in self.best_job.orders_id]  # MEMO: 1-indexedに変換
        print(len(r), *r)
        tour_coord = [self.get_coord_by_id(id, dest) for id, dest in self.best_job.tour]
        if DEBUG:
            print("format: n x_1 y_1 ... x_y y_n")
            print(len(tour_coord), *tour_coord)
            return
        print(len(tour_coord), end=" ")
        for tc in tour_coord:
            print(*tc, end=" ")
        print("")
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
        print("final objective value:", solver.best_job.obj)
        solver.check_time_limit()


if __name__ == '__main__':
    args = sys.argv
    DEBUG = bool("debug" in args)
    main()
