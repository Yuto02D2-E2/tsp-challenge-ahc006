#include <bits/stdc++.h>
using std::cin;
using std::cout;
using std::endl;
namespace chrono = std::chrono;

#ifdef _LOCAL
const bool LOCAL = true;
#else
const bool LOCAL = false;
#endif

// print out iterator object (sep=",",  end="\n")
#define dump(iter_) print__(#iter_, iter_)
namespace {
template <typename T>
void print__(const std::string& name_, const T& iter_) {
    std::cout << name_ << ":{";
    std::copy(iter_.begin(), iter_.end(), std::ostream_iterator<decltype(*iter_.begin())>(std::cout, ","));
    std::cout << "}" << std::endl;
}
}  // namespace

// 座標(0-indexed)
struct Point {
    int x, y;
    Point() { x = -1, y = -1; }
    Point(const int x_, const int y_) : x(x_), y(y_) {}
};

// destination: 行き先
enum class Dest {
    from,    // 受け取り先のレストラン
    to,      // 届け先
    office,  // start/goal
};

// 注文
struct Order {
    int id;
    Dest dest;
    Order(const int id_, const Dest dest_) : id(id_), dest(dest_) {}
};

// 問題で与えられた/問題から得られるデータ
struct Data {
    const int MAP_SIZE = 800;
    const Point OFFICE_POINT = Point(400, 400);
    const int ALL_ORDER_NUM = 1000;
    const int ORDER_NUM = 50;
    const int TOUR_LEN = ORDER_NUM * 2 + 2;
    std::vector<Point> from = std::vector<Point>(ALL_ORDER_NUM);
    std::vector<Point> to = std::vector<Point>(ALL_ORDER_NUM);
    Data() {}

    void read_data() {
        int a, b, c, d;
        for (int id = 0; id < ALL_ORDER_NUM; id++) {
            cin >> a >> b >> c >> d;
            assert(0 <= id && id < int(from.size()));
            assert(0 <= id && id < int(to.size()));
            from[id] = Point(a, b);
            to[id] = Point(c, d);
        }
        return;
    }
};

// 解の単位
struct Job {
    std::vector<int> orders_id;
    std::vector<Order> tour;
    int obj;  // objective value -> minimize
    Job() {}

    void init(const Data& data) {
        orders_id = std::vector<int>(data.ORDER_NUM, -1);
        tour = std::vector<Order>(data.TOUR_LEN, Order(-1, Dest::office));
        obj = (1 << 30);
        return;
    }
};

// インターフェース
class Solver {
private:
    const double TIME_LIMIT = 1950.0;  // [msec]
    chrono::system_clock::time_point start_time, cur_time;
    Data data;
    Job job;

public:
    Solver() {
        start_time = chrono::system_clock::now();
        data.read_data();
        job.init(data);
    }

    void select_order() {
        // officeからのマンハッタン距離が近そうな上位50件を受け付ける．
        // fromだけ,toだけが近いものを排除するために，遠すぎる座標に対してペナルティを課している
        std::vector<std::pair<int, int>> manhattan_dist(data.ALL_ORDER_NUM);
        const int threshold = data.MAP_SIZE / 4;
        for (int id = 0; id < data.ALL_ORDER_NUM; id++) {
            int from_dist = eval_dist(data.OFFICE_POINT, data.from[id]);
            int to_dist = eval_dist(data.OFFICE_POINT, data.to[id]);
            // 遠すぎる(閾値よりも遠い)やつにペナルティ
            if (threshold < from_dist) from_dist *= 5;
            if (threshold < to_dist) to_dist *= 5;
            manhattan_dist[id] = std::make_pair(from_dist + to_dist, id);
        }
        std::sort(manhattan_dist.begin(), manhattan_dist.end());
        for (int i = 0; i < data.ORDER_NUM; i++) {
            job.orders_id[i] = manhattan_dist[i].second;
        }
        return;
    }

    void nearest_neighbor() {
        // 方針：
        // fromを全て回ってから，toを全て回る
        // from/toそれぞれをk-meansでクラスタに分けて，それぞれの中でnnをする
        int tour_i = 1;
        // from
        const std::vector<std::vector<int>> from_clusters = k_means(data.from);
        std::multiset<int> from_done;
        // clusterを回る順番
        std::vector<int> rotate_indexes = {4, 1, 2, 5, 8, 7, 6, 3, 0};
        for (const int ri : rotate_indexes) {
            std::vector<int> one_cluster = from_clusters[ri];
            for (const int i_id : one_cluster) {
                int min_dist = std::numeric_limits<int>::max();
                int min_dist_j_id = -1;
                for (const int j_id : one_cluster) {
                    if (from_done.find(j_id) != from_done.end()) continue;
                    int cur_dist = eval_dist(data.from[i_id], data.from[j_id]);
                    if (cur_dist < min_dist) {
                        min_dist = cur_dist;
                        min_dist_j_id = j_id;
                    }
                }
                if (min_dist_j_id != -1) {
                    job.tour[tour_i] = Order(min_dist_j_id, Dest::from);
                    tour_i++;
                    from_done.insert(min_dist_j_id);
                }
                // print();
            }
        }
        // to
        const std::vector<std::vector<int>> to_clusters = k_means(data.to);
        std::multiset<int> to_done;
        // fromとは逆順に回る
        std::reverse(rotate_indexes.begin(), rotate_indexes.end());
        for (const int ri : rotate_indexes) {
            std::vector<int> one_cluster = to_clusters[ri];
            for (const int i_id : one_cluster) {
                int min_dist = std::numeric_limits<int>::max();
                int min_dist_j_id = -1;
                for (const int j_id : one_cluster) {
                    if (to_done.find(j_id) != to_done.end()) continue;
                    int cur_dist = eval_dist(data.to[i_id], data.to[j_id]);
                    if (cur_dist < min_dist) {
                        min_dist = cur_dist;
                        min_dist_j_id = j_id;
                    }
                }
                if (min_dist_j_id != -1) {
                    job.tour[tour_i] = Order(min_dist_j_id, Dest::to);
                    tour_i++;
                    to_done.insert(min_dist_j_id);
                }
                // print();
            }
        }
        return;
    }

    void local_search() {
        // TODO: 実装する
        // メモ：回り方を変えるか，注文を変えるかどっちにしよう？
        // while (check_time_limit()) {
        //   two_opt_search();
        //   or_opt_search();
        // }
        return;
    }

    void print() {
        std::vector<int> orders_id_;
        for (int id : job.orders_id) {
            orders_id_.emplace_back(id + 1);  // 1-indexed
        }
        std::vector<int> tour_;
        for (auto t : job.tour) {
            Point p = get_point_by_order(t);
            tour_.emplace_back(p.x);
            tour_.emplace_back(p.y);
        }
        if (LOCAL) {
            cout << "[result]" << endl;
            cout << "orders len:" << orders_id_.size() << endl;
            dump(orders_id_);
            cout << "tour len:" << job.tour.size() << endl;
            dump(tour_);
            cout << "score:" << get_current_score() << endl;
            cout << "objective value:" << job.obj << endl;
            std::cerr << get_current_score() << endl;
        } else {
            cout << orders_id_.size() << " ";
            std::copy(orders_id_.begin(), orders_id_.end(),
                      std::ostream_iterator<decltype(*orders_id_.begin())>(std::cout, " "));
            cout << endl;
            cout << job.tour.size() << " ";
            std::copy(tour_.begin(), tour_.end(), std::ostream_iterator<decltype(*tour_.begin())>(std::cout, " "));
            cout << endl;
        }
        return;
    }

private:
    bool check_time_limit() {
        cur_time = chrono::system_clock::now();
        double elapsed_time = chrono::duration_cast<chrono::milliseconds>(cur_time - start_time).count();
        return (elapsed_time < TIME_LIMIT);
    }

    int eval_dist(const Point& p, const Point& q) const {
        // p(x,y)とq(x,y)のマンハッタン距離を返す
        return (abs(p.x - q.x) + abs(p.y - q.y));
    }

    void calc_obj() {
        job.obj = 0;
        for (int i = 0; i < data.TOUR_LEN; i++) {
            job.obj +=
                eval_dist(get_point_by_order(job.tour[i]), get_point_by_order(job.tour[(i + 1) % data.TOUR_LEN]));
        }
        return;
    }

    double get_current_score() {
        calc_obj();
        return (1e8 / (1000.0 + job.obj));
    }

    Point get_point_by_order(const Order& order) const {
        // assert(-1 <= order.id && order.id < data.ORDER_NUM);
        if (order.dest == Dest::from) {
            return data.from[order.id];
        } else if (order.dest == Dest::to) {
            return data.to[order.id];
        } else {
            return data.OFFICE_POINT;
        }
    }

    std::vector<std::vector<int>> k_means(const std::vector<Point>& points) {
        // この関数は厳密にはk-meansでは無い．k-meansっぽいことをしてクラスターに分けている．
        constexpr int K = 9;  // クラスターの個数
        // return: clusters[i] := i番目のクラスター内に含まれる注文idの集合
        std::vector<std::vector<int>> clusters(K, std::vector<int>());
        // 左上を0として
        // 0,1,2
        // 3,4,5
        // 6,7,8
        // とラベリングする
        // 各ラベルの中心を基準点としてクラスタリング．
        // NNでは4(start)->1->2->5->8->7->6->3->0の順に(toはreverseする)回ることで良い感じになる予定(予定..)
        // 基準座標(x, y) (reference position) -> i番目のクラスターの基準座標
        std::vector<Point> ref_pos = {
            Point(300, 300), Point(400, 300), Point(500, 300), Point(300, 400), Point(400, 400),
            Point(500, 400), Point(300, 500), Point(400, 500), Point(500, 500),
        };
        for (int id : job.orders_id) {
            int nearest_k = -1, nearest_k_dist = std::numeric_limits<int>::max();
            for (int k = 0; k < K; k++) {
                int cur_dist = eval_dist(ref_pos[k], points[id]);
                if (cur_dist < nearest_k_dist) {
                    nearest_k = k;
                    nearest_k_dist = cur_dist;
                }
            }
            assert(0 <= nearest_k && nearest_k < K);
            clusters[nearest_k].emplace_back(id);
        }
        if (LOCAL) {
            for (auto one_cluster : clusters) {
                dump(one_cluster);
            }
        }
        return clusters;
    }

    void two_opt_search() {}
    void or_opt_search() {}
};

int main() {
    // std I/Oの高速化
    std::cin.tie(nullptr);
    std::ios::sync_with_stdio(false);

    Solver solver;
    solver.select_order();
    solver.nearest_neighbor();
    solver.local_search();
    solver.print();
    return 0;
}
