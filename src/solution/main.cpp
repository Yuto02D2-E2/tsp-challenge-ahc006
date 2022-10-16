// #include <bits/stdc++.h>
// 上記はg++だと使えるが，clang++だと使えない．
// 代わりに，bits/stdc++.hに書かれている中身を*.cppファイルの先頭に書けば良い
// 面倒くさいしg++で良いかと思うが，-fsanitize=undefined,addressはclang++しか使えない？っぽいのでclangを使いたい

// C
#ifndef _GLIBCXX_NO_ASSERT
#include <cassert>
#endif
#include <cctype>
#include <cerrno>
#include <cfloat>
#include <ciso646>
#include <climits>
#include <clocale>
#include <cmath>
#include <csetjmp>
#include <csignal>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cwchar>
#include <cwctype>

#if __cplusplus >= 201103L
#include <ccomplex>
#include <cfenv>
#include <cinttypes>
#include <cstdbool>
#include <cstdint>
#include <ctgmath>
#include <cuchar>
#endif

// C++
#include <algorithm>
#include <bitset>
#include <complex>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <iterator>
#include <limits>
#include <list>
#include <locale>
#include <map>
#include <memory>
#include <new>
#include <numeric>
#include <ostream>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <typeinfo>
#include <utility>
#include <valarray>
#include <vector>

#if __cplusplus >= 201103L
#include <array>
#include <atomic>
#include <chrono>
#include <codecvt>
#include <condition_variable>
#include <forward_list>
#include <future>
#include <initializer_list>
#include <mutex>
#include <random>
#include <ratio>
#include <regex>
#include <scoped_allocator>
#include <system_error>
#include <thread>
#include <tuple>
#include <typeindex>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#endif

#if __cplusplus >= 201402L
#include <shared_mutex>
#endif

#if __cplusplus >= 201703L
#include <any>
#include <charconv>
#include <execution>
#include <filesystem>
#include <optional>
#include <string_view>
#include <variant>
#endif

#if __cplusplus >= 202002L
#include <barrier>
#include <bit>
#include <compare>
#include <concepts>
#if __cpp_impl_coroutine
# include <coroutine>
#endif
#include <latch>
#include <numbers>
#include <ranges>
#include <span>
#include <stop_token>
#include <semaphore>
#include <source_location>
#include <syncstream>
#include <version>
#endif

#if __cplusplus > 202002L
#include <expected>
#include <spanstream>
#if __has_include(<stacktrace>)
# include <stacktrace>
#endif
#include <stdatomic.h>
#endif

// end of bits/stdc++

using std::cin;
using std::cout;
using std::endl;
namespace chrono = std::chrono;

#ifdef _LOCAL
const bool LOCAL = true;
#else
const bool LOCAL = false;
#endif

namespace writer {
/* print (iter, sep=",", end="\n") */
template <typename T>
void print(const T& iter_, const std::string& sep_ = ",", const std::string& end_ = "\n") {
    std::copy(iter_.begin(), iter_.end(), std::ostream_iterator<decltype(*iter_.begin())>(std::cout, sep_.c_str()));
    std::cout << end_;
    std::flush(std::cout);
}
/* print full (iter, msg="", sep=",") */
template <typename T>
void printf(const T& iter_, const std::string& info_ = "", const std::string& sep_ = ",") {
    std::cout << info_ << ":{";
    writer::print(iter_, sep_, "");
    std::cout << "}" << std::endl;
}
/* print initialize list (init_list, msg="", sep=",") */
template <typename T>
void printil(const std::initializer_list<T>& iter_, const std::string& info_ = "", const std::string& sep_ = ",") {
    std::cout << info_ << ":{";
    writer::print(iter_, sep_, "");
    std::cout << "}" << std::endl;
}
}  // namespace writer

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
    const int MAP_SIZE = 800;  // width=height=MAP_SIZE
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
    std::vector<Order> tour;        // index -> Order
    std::map<Order, int> order2id;  // Order -> index
    int obj;                        // objective value -> minimize
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
    // (dist, id)
    std::vector<std::pair<int, int>> manhattan_dist;

public:
    Solver() {
        start_time = chrono::system_clock::now();
        data.read_data();
        job.init(data);
    }

    void select_order() {
        // officeからのマンハッタン距離が近そうな上位50件を受け付ける．
        // fromだけ,toだけが近いものを排除するために，遠すぎる座標に対してペナルティを課している
        const int threshold = data.MAP_SIZE / 4;
        for (int id = 0; id < data.ALL_ORDER_NUM; id++) {
            int from_dist = get_dist_by_order(data.OFFICE_POINT, data.from[id]);
            int to_dist = get_dist_by_order(data.OFFICE_POINT, data.to[id]);
            // 遠すぎる(閾値よりも遠い)やつにペナルティ
            if (threshold < from_dist) from_dist *= 5;
            if (threshold < to_dist) to_dist *= 5;
            manhattan_dist.emplace_back(std::make_pair(from_dist + to_dist, id));
        }
        std::sort(manhattan_dist.begin(), manhattan_dist.end());
        for (int i = 0; i < data.ORDER_NUM; i++) {
            job.orders_id[i] = manhattan_dist[i].second;
        }
        return;
    }

    void build_init_solution() {
        // 初期解の作成
        // pop(degs) -> min({deg, id})
        std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>> from_degs, to_degs;
        for (int i = 0; i < data.ORDER_NUM; i++) {
            int id = job.orders_id[i];
            from_degs.push(std::make_pair(get_deg_by_order(Order(id, Dest::from)), id));
            to_degs.push(std::make_pair(get_deg_by_order(Order(id, Dest::to)), id));
        }
        int i = 1;  // tour index
        while (!from_degs.empty()) {
            auto [deg, id] = from_degs.top();  // c++17以降で使えるunpack方法．c++14以前ならstd::tie()を使う
            from_degs.pop();
            job.tour[i++] = Order(id, Dest::from);
        }
        assert(from_degs.empty());
        while (!to_degs.empty()) {
            auto [deg, id] = to_degs.top();
            to_degs.pop();
            job.tour[i++] = Order(id, Dest::to);
        }
        assert(to_degs.empty());
        return;
    }

    void local_search() {
        // TODO: 実装する
        // メモ：回り方を変えるか，注文を変えるかどっちにしよう？
        // while (check_time_limit()) {
        //     two_opt_search();
        //     or_opt_search();
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
            writer::printf(orders_id_, "orders id");
            cout << "tour len:" << tour_.size() / 2 << endl;
            writer::printf(tour_, "tour");
            cout << "score:" << get_current_score() << endl;
            cout << "objective value:" << job.obj << endl;
        } else {
            cout << orders_id_.size() << " ";
            writer::print(orders_id_, " ");
            cout << tour_.size() / 2 << " ";
            writer::print(tour_, " ");
        }
        return;
    }

    double get_current_score() {
        calc_obj();
        return (1e8 / (1000.0 + job.obj));
    }

private:
    bool check_time_limit() {
        cur_time = chrono::system_clock::now();
        double elapsed_time = chrono::duration_cast<chrono::milliseconds>(cur_time - start_time).count();
        return (elapsed_time < TIME_LIMIT);
    }

    int get_dist_by_order(const Point& p, const Point& q) const {
        // p(x,y)とq(x,y)のマンハッタン距離を返す
        return (abs(p.x - q.x) + abs(p.y - q.y));
    }

    void calc_obj() {
        job.obj = 0;
        for (int i = 0; i < data.TOUR_LEN; i++) {
            job.obj += get_dist_by_order(get_point_by_order(job.tour[i]),
                                         get_point_by_order(job.tour[(i + 1) % data.TOUR_LEN]));
        }
        return;
    }

    Point get_point_by_order(const Order& order) const {
        if (order.dest == Dest::from) {
            return data.from[order.id];
        } else if (order.dest == Dest::to) {
            return data.to[order.id];
        } else {
            return data.OFFICE_POINT;
        }
    }

    double get_deg_by_order(const Order& order) {
        // orderの座標ベクトル(point)が中心(400, 400)を原点としてx軸となす角[rad]を返す
        Point p = get_point_by_order(order);
        // 中心(400,400)を原点とみなすために全ての点を-400,-400する
        p.x -= 400, p.y -= 400;
        // atan2(y, x) = arc tan(y/x) = theta[rad]. theta=[-pi,pi]
        return atan2(double(p.y), double(p.x));
    }

    void two_opt_search() {
        // TODO: 実装する
        // pythonで回り方を変えるのをやったけど芳しくなかったので，注文を変える方針を試してみる
        // manhattan_distの51番目から順番に，from/toそれぞれを最近近傍法で挿入してみて，良さげだったら採用する
        // 一番Edgeが長くて，officeから離れているnodeを消す
        // SA(アニーリング;所謂焼きなまし)をしても良いかも
        // for (int i = data.ORDER_NUM + 1; i < data.ALL_ORDER_NUM; i++) {
        //     int i_id = manhattan_dist[i].second;
        //     int del_id = -1;
        //     int del_dist = std::numeric_limits<int>::min();  // officeからのdist(from + to) + edge_len
        //     for (int j_id : job.orders_id) {
        //         int cur_dist = get_dist_by_order(data.OFFICE_POINT, get_point_by_order(Order(j_id, Dest::from))) +
        //                        get_dist_by_order(data.OFFICE_POINT, get_point_by_order(Order(j_id, Dest::to)));
        //         cur_dist += get_dist_by_order(get_point_by_order());
        //     }
        // }
        return;
    }
    void or_opt_search() {}
};

int main() {
    // std I/Oの高速化
    std::cin.tie(nullptr);
    std::ios::sync_with_stdio(false);
    // init random seed
    std::srand(0);

    Solver solver;
    solver.select_order();
    solver.build_init_solution();
    solver.local_search();
    solver.print();
    std::cerr << solver.get_current_score() << endl;  // pythonスクリプトによるスコア計算用
    return 0;
}
