import sys
import math
import time
import enum
import dataclasses
import heapq
import random


@dataclasses.dataclass
class Point:
    """ 0 <= x,y < 800 """
    x: int = -1
    y: int = -1


class Data:
    def __init__(self) -> None:
        # map size (width=height)
        self.MAP_SIZE: int = 800
        # office coordinate: 本社の座標/start/goal -> Point(x, y)
        self.OFFICE_POINT: Point = Point(400, 400)
        # threshold: 閾値
        self.ALL_ORDER_NUM: int = 1000
        self.ORDER_NUM: int = 50
        # tour length = order_num * 2(from/to) + 2(start/end)
        self.TOUR_LEN: int = self.ORDER_NUM * 2 + 2
        # from_[id] := id番目の注文の受取先 -> tuple(x,y)
        self.from_: list = [Point() for _ in range(self.ALL_ORDER_NUM)]
        # to_[id] := id番目の注文の配達先 -> tuple(x,y)
        self.to_: list = [Point() for _ in range(self.ALL_ORDER_NUM)]
        return

    def read_data(self) -> None:
        for id in range(self.ALL_ORDER_NUM):
            a, b, c, d = map(int, input().split())
            self.from_[id] = Point(a, b)
            self.to_[id] = Point(c, d)
        return


class Dest(enum.Enum):
    """ destination """
    OFFICE = -1
    FROM = 0
    TO = 1


@dataclasses.dataclass
class Order:
    id: int = -1
    dest: Dest = Dest.OFFICE


class Job:
    """ 解の単位 """

    def __init__(self, data: Data) -> None:
        """
        - idとi,j,k...の使い分け
            - idは注文のid(1000件の内何番目か)
            - i,j,kは経路内のindex
        """
        # selected orders id: 受ける注文のid列(1-indexed)
        self.orders_id: list = [-1 for _ in range(data.ORDER_NUM)]
        # tour: 移動経路 -> Order(order_id, order_destination:from/to). これが解
        self.tour: list = [Order() for _ in range(data.TOUR_LEN)]
        # objective value: total tour distance 配達ルートの総移動距離 -> minimize
        self.obj: int = float("inf")
        return


class Solver:
    """ インターフェース的なクラス """

    """
    public function below
    """

    def __init__(self) -> None:
        self.TIME_LIMIT: float = 1.95  # [sec]
        self.__start_time: float = time.time()
        self.data: Data = Data()  # 問題で与えられた情報
        self.data.read_data()  # 入力の受け取り
        self.job: Job = Job(self.data)  # 最適解
        return

    def get_elapsed_time(self) -> float:
        return (time.time() - self.__start_time)

    def get_current_score(self) -> float:
        return (10**8 / (1000 + self.job.obj))

    def select_order(self) -> None:
        """受ける注文を選択(1000 -> 50)
        ### 方針
        [中心(office)からfrom,toへのマンハッタン距離の和+ペナルティ]が小さい上位50件を採用する．
        遠すぎるものにペナルティを課すことで，(fromが近くてtoが遠い)といったものを排除している．
        """
        # TODO: 半円ごとにfrom/toを選んでみる
        threshold = self.data.MAP_SIZE / 4 * 1.1
        # manhattan_dist[i]=tuple(dist,id)
        manhattan_dist = [tuple() for _ in range(self.data.ALL_ORDER_NUM)]
        for id in range(self.data.ALL_ORDER_NUM):
            from_dist = self.__get_dist_by_point(self.data.OFFICE_POINT, self.data.from_[id])
            to_dist = self.__get_dist_by_point(self.data.OFFICE_POINT, self.data.to_[id])
            # 閾値よりも遠いものにはペナルティを課す
            if threshold < from_dist:
                from_dist *= 5
            if threshold < to_dist:
                to_dist *= 5
            manhattan_dist[id] = (from_dist + to_dist, id)
        manhattan_dist.sort()
        self.job.orders_id = [id for _, id in manhattan_dist[:self.data.ORDER_NUM]]
        return

    def build_init_solution(self) -> None:
        """ 最近近傍法(貪欲法) - 初期解の生成．fromを全部回ってからtoを全部回る
        - TODO: GRASP(近いものの上位n件から確率で選ぶ)にしても良いかも
        """
        from_degs = list()
        to_degs = list()
        for id in self.job.orders_id:
            heapq.heappush(from_degs, (
                self.__get_rad_by_order(Order(id, Dest.FROM)),
                id
            ))
            heapq.heappush(to_degs, (
                self.__get_rad_by_order(Order(id, Dest.TO)),
                id
            ))
        tour_index = 1  # index=0,-1 -> office
        while from_degs:
            _, id = heapq.heappop(from_degs)
            self.job.tour[tour_index] = Order(id, Dest.FROM)
            tour_index += 1
        assert tour_index == self.data.ORDER_NUM + 1, tour_index
        while to_degs:
            _, id = heapq.heappop(to_degs)
            self.job.tour[tour_index] = Order(id, Dest.TO)
            tour_index += 1
        assert tour_index == self.data.TOUR_LEN - 1, tour_index
        self.__calc_obj(self.job)
        return

    def local_search(self) -> None:
        """ 局所探索法
        """
        # while self.__check_time_limit():
        #     self.__two_opt_search()
        return

    def print(self) -> None:
        """ 解の出力
        """
        r = [i + 1 for i in self.job.orders_id]  # MEMO: 1-indexedに変換
        tour_points = [self.__get_point_by_order(order) for order in self.job.tour]
        self.__calc_obj(self.job)
        if DEBUG:
            print("[answer]")
            print("format: m r_1 ... r_m \t(1-indexed)")
            print(len(r), *r)
            print("format: n x_1 y_1 ... x_y y_n")
            print(len(tour_points), *tour_points)
            print("current objective value:", self.job.obj)
            print("current score:", self.get_current_score())
        else:
            print(len(r), *r)
            tour_points = [f"{p.x} {p.y}" for p in tour_points]
            print(len(tour_points), " ".join(tour_points))
        return

    """
    private function below
    """

    def __check_time_limit(self) -> bool:
        return self.get_elapsed_time() < self.TIME_LIMIT

    def __get_dist_by_point(self, p: Point, q: Point) -> int:
        """ p(x,y)とq(x,y)のマンハッタン距離を返す
        ### 注意：ユークリッド距離では無い
        """
        assert (type(p) is Point and type(q) is Point), (p, q)
        return (abs(p.x - q.x) + abs(p.y - q.y))

    def __get_point_by_order(self, order: Order) -> Point:
        assert type(order) is Order, order
        if order.dest == Dest.FROM:
            return self.data.from_[order.id]
        elif order.dest == Dest.TO:
            return self.data.to_[order.id]
        else:
            return self.data.OFFICE_POINT

    def __get_rad_by_order(self, order: Order) -> float:
        """ orderの座標ベクトル(point)が中心(400, 400)を原点としてx軸となす角rad=[-pi,pi]を返す
        """
        assert type(order) is Order, order
        p = self.__get_point_by_order(order)
        # atan2(y, x) = arc tan(y/x) = theta[rad]. theta=[-pi,pi]
        # 中心(400,400)を原点とみなすために全ての点を-400,-400する
        return math.atan2(p.y - 400, p.x - 400)

    def __calc_obj(self, job: Job) -> None:
        """ jobの目的値(objective value):総移動距離を計算．
        """
        assert type(job) is Job, job
        job.obj = 0
        for i in range(self.data.TOUR_LEN):
            job.obj += self.__get_dist_by_point(
                self.__get_point_by_order(job.tour[i]),
                self.__get_point_by_order(job.tour[(i + 1) % self.data.TOUR_LEN]),
            )
        return

    def __two_opt_search(self) -> None:
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
            if not self.__check_time_limit("2-opt"):
                # タイムリミット
                return
        return


def main() -> None:
    random.seed(0)  # debugのために固定する
    solver = Solver()
    solver.select_order()
    solver.build_init_solution()
    # solver.print()
    solver.local_search()
    solver.print()
    if DEBUG:
        print("process time:", solver.get_elapsed_time())
    print(solver.get_current_score(), file=sys.stderr)


if __name__ == '__main__':
    args = sys.argv
    DEBUG = bool("debug" in args)
    main()
