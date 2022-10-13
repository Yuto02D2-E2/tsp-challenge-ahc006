/*
分かってはいたけど，やっぱりPythonだと速度的に限界がある．
60secくらいあればいいけど，2secだと殆ど何もできない
一旦C++でPythonの実装をコピーしてみる
*/


#include <bits/stdc++.h>
using std::cin;
using std::cout;
using std::endl;
namespace chrono = std::chrono;


// 座標(0-indexed)
struct Point {
    int x, y;
    Point() {
        x = -1, y = -1;
    }
    Point(const int x_, const int y_) : x(x_), y(y_) {}
};

// destination: 行き先
enum class Dest {
    from, // 受け取り先のレストラン
    to, // 届け先
    office, // start/goal
};

// 注文
struct Order {
    int id;
    Dest dest;
    Order(const int id_, const Dest dest_) : id(id_), dest(dest_) {}
};

// 問題で与えられた/問題から得られるデータ
struct Data {
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
    int obj; // objective value -> minimize
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
    const double TIME_LIMIT = 1950.0; // [msec]
    chrono::system_clock::time_point start_time_, cur_time_;
    Data data;
    Job job;

public:
    Solver() {
        start_time_ = chrono::system_clock::now();
        data.read_data();
        job.init(data);
    }

    void select_order() {
        std::vector<std::pair<int, int>> manhattan_dist(data.ALL_ORDER_NUM);
        const int threshold = 200;
        for (int id = 0; id < data.ALL_ORDER_NUM; id++) {
            int from_dist = eval_dist(data.OFFICE_POINT, data.from[id]);
            int to_dist = eval_dist(data.OFFICE_POINT, data.to[id]);
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

    void nearest_neighbor() {}
    void local_search() {}
    void print() {
        cout << "[answer]" << endl;
        cout << "orders len:" << job.orders_id.size() << endl;
        cout << "orders id:" << endl;
        for (int id : job.orders_id) {
            cout << "\tid:" << id << endl;
        }
        cout << "score:" << get_current_score() << endl;
        cout << "objective value:" << job.obj << endl;
        return;
    }

private:
    bool check_time_limit() {
        cur_time_ = chrono::system_clock::now();
        double elapsed_time = chrono::duration_cast<chrono::milliseconds>(cur_time_ - start_time_).count();
        return (elapsed_time < TIME_LIMIT);
    }

    int eval_dist(const Point& p, const Point& q) const {
        // p(x,y)とq(x,y)のマンハッタン距離を返す
        return (abs(p.x - q.x) + abs(p.y - q.y));
    }

    void calc_obj() {
        job.obj = 0;
        for (int i = 0; i < data.TOUR_LEN; i++) {
            job.obj += eval_dist(
                get_point_by_order(job.tour[i]),
                get_point_by_order(job.tour[(i + 1) % data.TOUR_LEN])
            );
        }
        return;
    }

    double get_current_score() {
        calc_obj();
        return (1e8 / 1000.0 + job.obj);
    }

    Point get_point_by_order(const Order& order) const {
        if (order.dest == Dest::office) {
            return data.OFFICE_POINT;
        } else if (order.dest == Dest::from) {
            return data.from[order.id];
        } else if (order.dest == Dest::to) {
            return data.to[order.id];
        } else {
            cout << "[get_point_by_order] invalid args" << endl;
            exit(0);
        }
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
