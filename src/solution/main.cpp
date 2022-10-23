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

// RNG; rondom number generator
struct XorShift {
private:
    std::uint64_t x = 123456789;
    std::uint64_t y = 362436069;
    std::uint64_t z = 521288629;
    std::uint64_t w = 88675123;
    std::uint32_t _max = std::numeric_limits<std::uint32_t>::max();
    std::uint32_t _min = std::numeric_limits<std::uint32_t>::min();

public:
    XorShift() {}

    double rand() {
        // [0,1]の(double)numberを返す
        std::uint32_t t = (x ^ (x << 11));
        x = y;
        y = z;
        z = w;
        w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
        return double(w - _min) / (_max - _min);  // min-max normalization
    }

    double randrange(const int& l, const int& r) {
        // [l,r]の(double)numberを返す
        return l + rand() * (r - l);
    }
};

// 実行時間管理
struct Timer {
private:
    // const std::int64_t TIME_LIMIT = 5000;  // DEBUG:
    const std::int64_t TIME_LIMIT = 1950;  // [msec]
    chrono::system_clock::time_point start_time, cur_time;

public:
    Timer() {
        start_time = chrono::system_clock::now();
    }

    std::int64_t get_time_limit() const {
        return TIME_LIMIT;
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
    // あるidのorderのtour内での順番(index)
    std::vector<int> from_id_index;
    std::vector<int> to_id_index;

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
        from_id_index = std::vector<int>(data.ALL_ORDER_NUM, -1);
        to_id_index = std::vector<int>(data.ALL_ORDER_NUM, -1);
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
        for (int i = 0; i < int(from_id_index.size()); i++) {
            from_id_index[i] = one_job.from_id_index[i];
            to_id_index[i] = one_job.to_id_index[i];
        }
        return;
    }

    void set_order(const int& index, const Order& order) {
        // tour[index]=orderに設定
        assert(0 <= index && index < int(tour.size()));
        tour[index] = order;
        set_index_by_order(index, order);
        return;
    }

    int get_index_by_order(const Order& order) const {
        // あるorderのtour内の順番(index)を返す
        if (order.dest == Dest::from) {
            return from_id_index[order.id];
        } else if (order.dest == Dest::to) {
            return to_id_index[order.id];
        } else {
            return 0;  // office
        }
    }

    void insert(const int& index, const Order& order) {
        // イメージ
        // t=[a,b,c,d,e,f,g], t[1]=eの場合
        // shift(1,4)
        // t=[a,b,b,c,d,f,g]
        // insert(e)
        // t=[a,e,b,c,d,f,g]
        shift_tour(index, get_index_by_order(order));
        tour[index] = order;
        set_index_by_order(index, tour[index]);
        return;
    }

    void swap(const int& i, const int& j) {
        // tour[i]とtour[j]をswapする
        // Order tmp = tour[i];
        // tour[i] = tour[j];
        // tour[j] = tmp;
        std::swap(tour[i], tour[j]);  // MEMO: バグってる説ある
        set_index_by_order(i, tour[i]);
        set_index_by_order(j, tour[j]);
        return;
    }

private:
    void set_index_by_order(const int& index, const Order& order) {
        // あるorderのtour内の順番(index)をセットする
        assert(0 <= index && index < int(tour.size()));
        if (order.dest == Dest::from) {
            from_id_index[order.id] = index;
        } else if (order.dest == Dest::to) {
            to_id_index[order.id] = index;
        } else {
            return;
        }
    }

    void shift_tour(const int& top_index, const int& tail_index) {
        // tourの[top,tail]を1つずつ後ろにずらす．tour[tail]を捨てる
        for (int i = tail_index; i > top_index; i--) {
            // 後ろから前にずらす
            tour[i] = tour[i - 1];
            set_index_by_order(i, tour[i]);
        }
        return;
    }
};

// インターフェース
class Solver {
private:
    Timer timer = Timer();            // 実行時間管理
    Data data = Data();               // data set
    Job best_job = Job();             // 暫定解(今まで調べた中で一番良かったもの．最適解とは限らない)
    Job cur_job = Job();              // 仮の解
    XorShift xs_random = XorShift();  // 乱数生成器
    const double PI = std::acos(-1.0);

public:
    /*
    public function below
    */

    Solver() {
        data.read_data();
        best_job.init(data);
        cur_job.init(data);
    }

    void select_order() {
        // officeからのマンハッタン距離が近そうな上位50件を受け付ける．
        // fromだけ,toだけが近いものを排除するために，遠すぎる座標に対してペナルティを課している
        const int threshold = 250;
        // manhattan_dist[id] -> (from+to dist, id)
        std::vector<std::pair<int, int>> manhattan_dist;
        for (int id = 0; id < data.ALL_ORDER_NUM; id++) {
            int from_dist = get_dist_by_point(data.OFFICE_POINT, data.from[id]);
            int to_dist = get_dist_by_point(data.OFFICE_POINT, data.to[id]);
            // 遠すぎる(閾値よりも遠い)やつにペナルティ
            if (threshold < from_dist) from_dist *= 5;
            if (threshold < to_dist) to_dist *= 5;
            manhattan_dist.emplace_back(std::make_pair(
                from_dist + to_dist,
                id));
        }
        std::sort(manhattan_dist.begin(), manhattan_dist.end());
        for (int i = 0; i < data.ORDER_NUM; i++) {
            best_job.orders_id[i] = manhattan_dist[i].second;
        }
        return;
    }

    void build_init_solution() {
        // 初期解の作成
        // 中心からの角度[rad]が小さい順にfrom->toの順に回る．
        // pop(degs) -> min({deg, id})
        std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>> from_degs, to_degs;
        for (int i = 0; i < data.ORDER_NUM; i++) {
            int id = best_job.orders_id[i];
            from_degs.push(std::make_pair(
                get_rad_by_order(Order(id, Dest::from)),
                id));
            to_degs.push(std::make_pair(
                get_rad_by_order(Order(id, Dest::to)),
                id));
        }
        int tour_index = 1;  // tour index
        while (!from_degs.empty()) {
            auto [deg, id] = from_degs.top();  // c++17以降で使えるunpack方法．c++14以前ならstd::tie()を使う
            from_degs.pop();
            best_job.set_order(tour_index++, Order(id, Dest::from));
        }
        while (!to_degs.empty()) {
            auto [deg, id] = to_degs.top();
            to_degs.pop();
            best_job.set_order(tour_index++, Order(id, Dest::to));
        }
        best_job.obj = calc_obj(best_job);
        return;
    }

    void simulated_annealing() {
        // 焼きなまし法
        // 焼きなましのパラメータ → 温度;temperature
        double max_temp;  // TODO: 求める
        double min_temp;
        cur_job.set_copy(best_job);  // init
        while (timer.check_time_limit()) {
            // 幾何冷却スケジューリング
            // double _time = double(timer.get_elapsed_time()) / timer.get_time_limit();
            // double temp = pow(max_temp, 1.0 - _time) * pow(min_temp, _time);
            double temp = 0.0;  // DEBUG:
            nbhd_search(temp);
            if (cur_job.obj < best_job.obj) {
                best_job.set_copy(cur_job);
                print();
            }
        }
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
        return (1e8 / (1000.0 + calc_obj(best_job)));
    }

private:
    /*
    private function below
    */

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

    double get_rad_diff(const double& p_rad, const double& q_rad) const {
        // p[rad]とq[rad]のなす角[0,pi]を返す
        // p,q共に値域は[-pi,pi]
        // TODO: 合ってる？
        return std::min(abs(p_rad - q_rad), 2 * PI - (p_rad - q_rad));
    }

    int calc_obj(const Job& one_job) {
        // objective valueを計算して返す
        int obj = 0;
        for (int i = 0; i < data.TOUR_LEN; i++) {
            obj += get_dist_by_point(
                get_point_by_order(one_job.tour[i]),
                get_point_by_order(one_job.tour[(i + 1) % data.TOUR_LEN]));
        }
        return obj;
    }

    void nbhd_search(const double& temp) {
        // 近傍探索 (neighborhood search)
        /*
        近傍を「insert処理またはswap処理を一回実施した状態」と定義する
        insert処理 <=def=> tour[i]をtour[j]にinsertする．(j<i)
        swap処理 <=def=> tour[i]とtour[j]をswapする
        insert処理，swap処理共に，順序関係が壊れない場合にのみ実施する
        ただし，スコアが改悪した場合でも一定確率で採用する(焼きなましによって多様性を出したい)
        */

        std::vector<std::pair<Order, Order>> order_pairs;  // 近傍を全探索する時に使う
        for (const int& id_i : cur_job.orders_id) {
            for (const int& id_j : cur_job.orders_id) {
                if (id_i == id_j) {
                    continue;
                }
                order_pairs.emplace_back(std::make_pair(
                    Order(id_i, Dest::from),
                    Order(id_j, Dest::from)));
                order_pairs.emplace_back(std::make_pair(
                    Order(id_i, Dest::from),
                    Order(id_j, Dest::to)));
                order_pairs.emplace_back(std::make_pair(
                    Order(id_i, Dest::to),
                    Order(id_j, Dest::from)));
                order_pairs.emplace_back(std::make_pair(
                    Order(id_i, Dest::to),
                    Order(id_j, Dest::to)));
            }
        }

        const double RAD_EPS = PI / 6.0;   // 逆走の許容値．pi/6[rad]=30[deg]までなら逆走を許す
        int insert_cnt = 0, swap_cnt = 0;  // DEBUG:処理回数のメモ

        const int T = int(order_pairs.size());
        for (int t = 0; t < T; t++) {
            // for (const auto& [order_i, order_j] : order_pairs) {
            auto [order_i, order_j] = order_pairs[xs_random.randrange(0, T)];
            assert(order_i.id != order_j.id);
            int i = cur_job.get_index_by_order(order_i);
            int j = cur_job.get_index_by_order(order_j);
            if (i < j) {
                // iとjの順序関係．j<iを想定しているのでi<jはダメ
                continue;
            }
            if ((i <= 1) || (data.TOUR_LEN - 2 <= i) ||
                (j <= 1) || (data.TOUR_LEN - 2 <= j)) {
                // 端っこはダメ
                continue;
            }
            // if (false) {  // DEBUG:
            if (xs_random.rand() < 0.95) {  // 確率でinsert．swapを選択する
                // insert近傍
                if ((order_i.dest == Dest::to) &&
                    (j < cur_job.get_index_by_order(Order(order_i.id, Dest::from)))) {
                    // insertして順序関係を満たさなくなるならダメ
                    continue;
                }
                if ((RAD_EPS < get_rad_diff(get_rad_by_order(order_j), get_rad_by_order(order_i))) ||
                    (RAD_EPS < get_rad_diff(get_rad_by_order(order_i), get_rad_by_order(cur_job.tour[j + 1])))) {
                    // 角度が離れすぎると良くない?
                    continue;
                }
                Point j_point = get_point_by_order(order_j);
                Point next_j_point = get_point_by_order(cur_job.tour[j + 1]);
                Point prev_i_point = get_point_by_order(cur_job.tour[i - 1]);
                Point i_point = get_point_by_order(order_i);
                Point next_i_point = get_point_by_order(cur_job.tour[i + 1]);
                int cur_dist = get_dist_by_point(j_point, next_j_point) +
                               get_dist_by_point(prev_i_point, i_point) +
                               get_dist_by_point(i_point, next_i_point);
                int new_dist = get_dist_by_point(j_point, i_point) +
                               get_dist_by_point(i_point, prev_i_point) +
                               get_dist_by_point(prev_i_point, next_i_point);
                int delta = new_dist - cur_dist;
                if (delta < 0) {
                    // if ((delta < 0) || (xs_random.rand() < exp(-delta / temp))) {
                    // スコアが改善する or 一定確率 でこの近傍に移動
                    cur_job.insert(j + 1, order_i);
                    insert_cnt++;
                    // cur_job.obj = calc_obj(cur_job);
                    // return;
                }
            } else {
                // swap近傍
                if ((order_i.dest == Dest::to) &&
                    (j < cur_job.get_index_by_order(Order(order_i.id, Dest::from)))) {
                    // iは前に移動する
                    // swapしてi.fromとi.toが順序関係を満たさなくなるならダメ
                    continue;
                }
                if ((order_j.dest == Dest::from) &&
                    (cur_job.get_index_by_order(Order(order_j.id, Dest::to)) < i)) {
                    // jは後ろに移動する
                    // swapしてj.fromとj.toが順序関係を満たさなくなるならダメ
                    continue;
                }
                Point prev_j_point = get_point_by_order(cur_job.tour[j - 1]);
                Point j_point = get_point_by_order(order_j);
                Point next_j_point = get_point_by_order(cur_job.tour[j + 1]);
                Point prev_i_point = get_point_by_order(cur_job.tour[i - 1]);
                Point i_point = get_point_by_order(order_i);
                Point next_i_point = get_point_by_order(cur_job.tour[i + 1]);
                int cur_dist = get_dist_by_point(prev_j_point, j_point) +
                               get_dist_by_point(j_point, next_j_point) +
                               get_dist_by_point(prev_i_point, i_point) +
                               get_dist_by_point(i_point, next_i_point);
                int new_dist = get_dist_by_point(prev_j_point, i_point) +
                               get_dist_by_point(i_point, next_j_point) +
                               get_dist_by_point(prev_i_point, j_point) +
                               get_dist_by_point(j_point, next_i_point);
                int delta = new_dist - cur_dist;
                if (delta < 0) {  // DEBUG:
                    // if ((delta < 0) || (xs_random.rand() < exp(-delta / temp))) {
                    // スコアが改善する or 一定確率 でこの近傍に移動
                    cur_job.swap(i, j);
                    swap_cnt++;
                    // cur_job.obj = calc_obj(cur_job);
                    // return;
                }
            }
        }
        cur_job.obj = calc_obj(cur_job);
        if (LOCAL) {
            cout << "insert/swap cnt:" << insert_cnt << "/" << swap_cnt << endl;
            cout << "best/cur obj:" << best_job.obj << "/" << cur_job.obj << endl;
        }
        return;
    }
};

int main() {
    // std I/Oの高速化
    std::cin.tie(nullptr);
    std::ios::sync_with_stdio(false);

    Solver solver;
    solver.select_order();
    solver.build_init_solution();
    solver.print();
    solver.simulated_annealing();
    solver.print();
    std::cerr << solver.get_current_score() << endl;  // pythonスクリプトによるスコア計算用
    return 0;
}
