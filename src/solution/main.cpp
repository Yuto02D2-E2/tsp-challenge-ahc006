// #include <bits/stdc++.h>
// 上記はg++だと使えるが，clang++だと使えない．
// 代わりに，bits/stdc++.hに書かれている中身を*.cppファイルの先頭に書けば良い

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
#include <type_traits>
#include <typeindex>
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
#include <coroutine>
#endif
#include <latch>
#include <numbers>
#include <ranges>
#include <semaphore>
#include <source_location>
#include <span>
#include <stop_token>
#include <syncstream>
#include <version>
#endif

#if __cplusplus > 202002L
#include <expected>
#include <spanstream>
#if __has_include(<stacktrace>)
#include <stacktrace>
#endif
#include <stdatomic.h>
#endif

// end of bits/stdc++

#ifdef _LOCAL
const bool LOCAL = true;
#else
const bool LOCAL = false;
#endif

using std::cin;
using std::cout;
using std::endl;
namespace chrono = std::chrono;

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

// 実行時間計測
struct Timer {
private:
    const std::int64_t TIME_LIMIT = 1950;  // [msec]
    chrono::system_clock::time_point start_time, cur_time;

public:
    Timer() {
        start_time = chrono::system_clock::now();
    }

    std::int64_t get_elapsed_time() {
        // 現在までの経過時間elapsed -> (int64_t)time[msec]を返す
        cur_time = chrono::system_clock::now();
        return chrono::duration_cast<chrono::milliseconds>(cur_time - start_time).count();
    }

    bool check_time_limit() {
        return (get_elapsed_time() < TIME_LIMIT);
    }
};

// 座標(0-indexed)
struct Point {
public:
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
public:
    int id;
    Dest dest;
    Order(const int id_, const Dest dest_) : id(id_), dest(dest_) {}
};

// 問題で与えられた/問題から得られるデータ
struct Data {
public:
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
private:
    // あるidのorderのtour内での順番
    std::vector<int> from_id_pos;
    std::vector<int> to_id_pos;

public:
    // 選択したorderのid集合
    std::vector<int> orders_id;
    // 訪れる座標(point)順序付き集合
    std::vector<Order> tour;
    int obj;  // objective value -> minimize
    Job() {}

    void init(const Data& data) {
        obj = (1 << 30);
        orders_id = std::vector<int>(data.ORDER_NUM, -1);
        tour = std::vector<Order>(data.TOUR_LEN, Order(-1, Dest::office));
        from_id_pos = std::vector<int>(data.ALL_ORDER_NUM, -1);
        to_id_pos = std::vector<int>(data.ALL_ORDER_NUM, -1);
        return;
    }

    void set_copy(const Job& one_job) {
        // one_jobの内容をコピーする
        obj = one_job.obj;
        for (int i = 0; i < int(orders_id.size()); i++) {
            orders_id[i] = one_job.orders_id[i];
        }
        for (int i = 0; i < int(tour.size()); i++) {
            tour[i] = one_job.tour[i];
        }
        for (int i = 0; i < int(from_id_pos.size()); i++) {
            from_id_pos[i] = one_job.from_id_pos[i];
            to_id_pos[i] = one_job.to_id_pos[i];
        }
        return;
    }

    void set_index_by_order(const Order& order, const int& pos) {
        // あるorderのtour内の順番(index)をセットする
        assert(0 <= pos && pos < int(tour.size()));
        if (order.dest == Dest::from) {
            from_id_pos[order.id] = pos;
        } else if (order.dest == Dest::to) {
            to_id_pos[order.id] = pos;
        } else {
            return;
        }
    }

    int get_index_by_order(const Order& order) const {
        // あるorderのtour内の順番(index)を返す
        if (order.dest == Dest::from) {
            assert(from_id_pos[order.id] != -1);
            return from_id_pos[order.id];
        } else if (order.dest == Dest::to) {
            assert(to_id_pos[order.id] != -1);
            return to_id_pos[order.id];
        } else {
            return 0;  // office
        }
    }

    void insert(const int& target_index, const Order& order) {
        // FIXME: バグってる？
        shift_tour(target_index, get_index_by_order(order));
        tour[target_index] = order;
        return;
    }

private:
    void shift_tour(const int& top_index, const int& tail_index) {
        // tourの[top,tail]を1つずつ後ろにずらす．tour[tail]を捨てる
        for (int i = tail_index; i >= top_index; i--) {
            // 後ろから前にずらす
            tour[i] = tour[i - 1];
        }
        return;
    }
};

// インターフェース
class Solver {
private:
    Timer timer = Timer();  // 実行時間計測
    Data data = Data();     // data set
    Job best_job = Job();   // 最適解
    // (dist, id)
    std::vector<std::pair<int, int>> manhattan_dist;
    const double PI = std::acos(-1.0);

public:
    Solver() {
        data.read_data();
        best_job.init(data);
    }

    void select_order() {
        // officeからのマンハッタン距離が近そうな上位50件を受け付ける．
        // fromだけ,toだけが近いものを排除するために，遠すぎる座標に対してペナルティを課している
        const int threshold = data.MAP_SIZE / 4 * 1.1;
        for (int id = 0; id < data.ALL_ORDER_NUM; id++) {
            int from_dist = get_dist_by_point(data.OFFICE_POINT, data.from[id]);
            int to_dist = get_dist_by_point(data.OFFICE_POINT, data.to[id]);
            // 遠すぎる(閾値よりも遠い)やつにペナルティ
            if (threshold < from_dist) from_dist *= 5;
            if (threshold < to_dist) to_dist *= 5;
            manhattan_dist.emplace_back(std::make_pair(from_dist + to_dist, id));
        }
        std::sort(manhattan_dist.begin(), manhattan_dist.end());
        for (int i = 0; i < data.ORDER_NUM; i++) {
            best_job.orders_id[i] = manhattan_dist[i].second;
        }
        return;
    }

    void build_init_solution() {
        // 初期解の作成
        // pop(degs) -> min({deg, id})
        std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>> from_degs, to_degs;
        for (int i = 0; i < data.ORDER_NUM; i++) {
            int id = best_job.orders_id[i];
            from_degs.push(std::make_pair(get_rad_by_order(Order(id, Dest::from)), id));
            to_degs.push(std::make_pair(get_rad_by_order(Order(id, Dest::to)), id));
        }
        int i = 1;  // tour index
        while (!from_degs.empty()) {
            auto [deg, id] = from_degs.top();  // c++17以降で使えるunpack方法．c++14以前ならstd::tie()を使う
            from_degs.pop();
            best_job.tour[i] = Order(id, Dest::from);
            best_job.set_index_by_order(best_job.tour[i], i);
            i++;
        }
        assert(from_degs.empty());
        while (!to_degs.empty()) {
            auto [deg, id] = to_degs.top();
            to_degs.pop();
            best_job.tour[i] = Order(id, Dest::to);
            best_job.set_index_by_order(best_job.tour[i], i);
            i++;
        }
        assert(to_degs.empty());
        calc_obj(best_job);
        return;
    }

    void local_search() {
        // 局所探索
        hill_climbing();  // 山登り法; hill climbing
        // simulate_annealing(); // 焼きなまし法; Simulated Annealing
        // while (timer.check_time_limit()) {
        //     two_opt_search();
        //     or_opt_search();
        // }
        return;
    }

    void print() {
        std::vector<int> orders_id_;
        for (const int& id : best_job.orders_id) {
            orders_id_.emplace_back(id + 1);  // 1-indexed
        }
        std::vector<int> tour_;
        for (const Order& o : best_job.tour) {
            Point p = get_point_by_order(o);
            tour_.emplace_back(p.x);
            tour_.emplace_back(p.y);
        }
        if (LOCAL) {
            cout << "[result]" << endl;
            cout << "orders len:" << orders_id_.size() << endl;
            writer::printf(orders_id_, "orders id");
            cout << "tour len:" << tour_.size() / 2 << endl;
            writer::printf(tour_, "tour");
            cout << "objective value:" << best_job.obj << endl;
            cout << "score:" << get_current_score() << endl;
            cout << "process time:" << timer.get_elapsed_time() << "[msec]" << endl;
        } else {
            cout << orders_id_.size() << " ";
            writer::print(orders_id_, " ");
            cout << tour_.size() / 2 << " ";
            writer::print(tour_, " ");
        }
        return;
    }

    double get_current_score() {
        calc_obj(best_job);
        return (1e8 / (1000.0 + best_job.obj));
    }

private:
    int get_dist_by_point(const Point& p, const Point& q) const {
        // p(x,y)とq(x,y)のマンハッタン距離を返す
        return (abs(p.x - q.x) + abs(p.y - q.y));
    }

    Point get_point_by_order(const Order& order) const {
        // あるorderの座標(point)を返す
        if (order.dest == Dest::from) {
            return data.from[order.id];
        } else if (order.dest == Dest::to) {
            return data.to[order.id];
        } else {
            return data.OFFICE_POINT;
        }
    }

    double get_rad_by_order(const Order& order) {
        // orderの座標ベクトル(point)が中心(400, 400)を原点としてx軸となす角rad=[-pi,pi]を返す
        Point p = get_point_by_order(order);
        // 中心(400,400)を原点とみなすために全ての点を-400,-400する
        p.x -= 400, p.y -= 400;
        // atan2(y, x) = arc tan(y/x) = theta[rad]. theta=[-pi,pi]
        return atan2(double(p.y), double(p.x));
    }

    double get_rad_diff(const double& p_rad, const double& q_rad) {
        // TODO: 合ってる？
        return std::min(abs(p_rad - q_rad), 2 * PI - (p_rad - q_rad));
    }

    void calc_obj(Job& one_job) {
        // objective valueを計算する．
        one_job.obj = 0;
        for (int i = 0; i < data.TOUR_LEN; i++) {
            one_job.obj += get_dist_by_point(
                get_point_by_order(one_job.tour[i]),
                get_point_by_order(one_job.tour[(i + 1) % data.TOUR_LEN]));
        }
        return;
    }

    void hill_climbing() {
        // 山登り法
        // Dest::toのorderをDest::fromの間に挿入していく．つまり，一週目で配れる分は配っておく
        Job cur_job = Job();  // 山登りの結果格納用
        cur_job.init(data);
        cur_job.set_copy(best_job);
        const double RAD_EPS = PI / 4.0;  // 逆走の許容値．pi/4[rad]=45[deg]までなら逆走を許す
        const int threshold = data.MAP_SIZE / 4;
        for (const int& to_id : cur_job.orders_id) {
            Order to_order = Order(to_id, Dest::to);
            if ((cur_job.get_index_by_order(to_order) <= 0) || (data.TOUR_LEN - 2 <= cur_job.get_index_by_order(to_order))) {
                // 端っこ
                if (LOCAL) {
                    cout << "reach border" << endl;
                }
                continue;
            }
            Point to_point = get_point_by_order(to_order);
            int opt_delta = 0;  // score_deltaの絶対値が大きく，符号が負であれば良い
            int opt_delta_id = -1;
            for (int id = 1; id < data.ORDER_NUM; id++) {
                if (cur_job.get_index_by_order(cur_job.tour[id]) < cur_job.get_index_by_order(Order(to_id, Dest::from))) {
                    // 順序関係を満たしているかチェックする
                    // 常にfrom_pos < to_posが成立しないとダメ
                    if (LOCAL) {
                        cout << "illigal order" << endl;
                    }
                    continue;
                }
                if ((RAD_EPS < get_rad_diff(get_rad_by_order(cur_job.tour[id]), get_rad_by_order(to_order))) ||
                    (RAD_EPS < get_rad_diff(get_rad_by_order(to_order), get_rad_by_order(cur_job.tour[id + 1])))) {
                    // 離れすぎるとよくない
                    if (LOCAL) {
                        cout << "rad too large" << endl;
                    }
                    continue;
                }
                if ((threshold < get_dist_by_point(get_point_by_order(cur_job.tour[id]), get_point_by_order(to_order))) ||
                    (threshold < get_dist_by_point(get_point_by_order(to_order), get_point_by_order(cur_job.tour[id + 1])))) {
                    // 離れすぎるとよくない
                    if (LOCAL) {
                        cout << "dist too large" << endl;
                    }
                    continue;
                }
                Point to_prev_point = get_point_by_order(
                    cur_job.tour[cur_job.get_index_by_order(to_order) - 1]);
                Point to_next_point = get_point_by_order(
                    cur_job.tour[cur_job.get_index_by_order(to_order) + 1]);
                Point from_p = get_point_by_order(cur_job.tour[id]);
                Point from_q = get_point_by_order(cur_job.tour[id + 1]);
                int cur_dist = get_dist_by_point(to_prev_point, to_point) +
                               get_dist_by_point(to_point, to_next_point) +
                               get_dist_by_point(from_p, from_q);
                int new_dist = get_dist_by_point(from_p, to_point) +
                               get_dist_by_point(to_point, from_q) +
                               get_dist_by_point(to_prev_point, to_next_point);
                int cur_delta = new_dist - cur_dist;
                if (cur_delta < opt_delta) {
                    opt_delta = cur_delta;
                    opt_delta_id = id;
                }
            }
            if (0 < opt_delta_id) {
                // 良い感じに挿入できる場所があった
                if (LOCAL) {
                    cout << to_id << " is insert to " << opt_delta_id + 1 << endl;
                }
                cur_job.insert(opt_delta_id + 1, to_order);
                cur_job.set_index_by_order(to_order, opt_delta_id + 1);
            }
        }
        calc_obj(cur_job);
        if (LOCAL) {
            cout << "cur_job.obj:" << cur_job.obj << endl;
            cout << "best_job.obj:" << best_job.obj << endl;
        }
        if (cur_job.obj < best_job.obj) {
            // 山登りしたらscoreが改善(objが減少)した場合
            // best_job.set_copy(cur_job);
        }
        // best_job.set_copy(cur_job);  // debug
        return;
    }

    void simulate_annealing() {}
    void two_opt_search() {}
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
    solver.print();
    solver.local_search();
    solver.print();
    std::cerr << solver.get_current_score() << endl;  // pythonスクリプトによるスコア計算用
    return 0;
}
