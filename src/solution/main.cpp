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


// 座標を表す
struct Point {
    int x, y;
    Point() {
        x = -1, y = -1;
    }
    Point(const int x_, const int y_) : x(x_), y(y_) {}
};

// 注文を表す
struct Order {
    int id;
    std::string dest; // destination
    Order(const int id_, const std::string dest_) : id(id_), dest(dest_) {}
};

//
struct Data {
    const Point OFFICE_POINT = Point(400, 400);
    const int ALL_ORDER_NUM = 1000;
    const int ORDER_NUM = 50;
    const int TOUR_LEN = ORDER_NUM * 2 + 2;
    std::vector<Point> from = std::vector<Point>(ALL_ORDER_NUM);
    std::vector<Point> to = std::vector<Point>(ALL_ORDER_NUM);

    // Data() {}

    void read_data() {
        int a, b, c, d;
        for (int id = 0; id < ALL_ORDER_NUM; id++) {
            cin >> a >> b >> c >> d;
            assert(0 <= id && id < int(from.size()));
            assert(0 <= id && id < int(to.size()));
            from[id] = Point(a, b);
            to[id] = Point(c, d);
        }
    }
};


// 解の単位
struct Job {
    std::vector<int> orders_id;
    std::vector<Order> tour;
    int obj; // objective value -> minimize

    // Job() {}
    Job(const Data data) {
        orders_id = std::vector<int>(data.ORDER_NUM, -1);
        tour = std::vector<Order>(data.ORDER_NUM, Order(-1, "none"));
        obj = (1 << 30);
    }
};


class Solver {
private:
    const double TIME_LIMIT = 1950.0; // [msec]
    chrono::system_clock::time_point start_time_, cur_time_;
    Data data;
    // MEMO: structを書かないとエラーになる．無くてもいけた気がするけど無理なの？
    // FIXME: solver内からjobを参照できない．job.orders_idとするとエラーになる．なぜ...
    Job job(struct data);

public:
    Solver() {
        start_time_ = chrono::system_clock::now();
        data.read_data();
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
        // for (int i = 0; i < data.ORDER_NUM; i++) {
        //     job.orders_id[i] = manhattan_dist[i].second;
        // }
    }

    void nearest_neighbor() {}
    void local_search() {}
    void print() {
        cout << "[answer]" << endl;
        // MEMO: エラー出る．わけ分からん．C++難しい....
        // cout << "orders len:" << job.orders_id.size() << endl;
        // cout << "orders id:" << endl;
        // for (int id : job.orders_id) {
        //     cout << "id:" << id << endl;
        // }
    }

private:
    bool check_time_limit() {
        cur_time_ = chrono::system_clock::now();
        double elapsed_time = chrono::duration_cast<chrono::milliseconds>(cur_time_ - start_time_).count();
        return (elapsed_time < TIME_LIMIT);
    }

    int eval_dist(const Point p, const Point q) {
        // p(x,y)とq(x,y)のマンハッタン距離を返す
        return (abs(p.x - q.x) + abs(p.y - q.y));
    }

    void calc_obj() {
        // job.obj = 0;
        // for (int i = 0; i < data.TOUR_LEN; i++) {
        //     job.obj += eval_dist(
        //         get_point_by_id(job.tour[i]),
        //         get_point_by_id(job.tour[(i + 1) % data.TOUR_LEN])
        //     );
        // }
        return;
    }

    // double get_current_score() {
    //     calc_obj();
    //     return (1e8 / 1000.0 + job.obj);
    // }

    Point get_point_by_id(const int id, const std::string dest) {
        if (id == -1) {
            return data.OFFICE_POINT;
        }
        if (dest == "from") {
            return data.from[id];
        } else if (dest == "to") {
            return data.to[id];
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
