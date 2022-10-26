// #include <bits/stdc++.h>

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

// RNG; random number generator; 乱数生成器
struct RNG {
private:
    std::uint64_t x = 123456789;
    std::uint64_t y = 362436069;
    std::uint64_t z = 521288629;
    std::uint64_t w = 88675123;
    std::uint32_t max_ = std::numeric_limits<std::uint32_t>::max();
    std::uint32_t min_ = std::numeric_limits<std::uint32_t>::min();

public:
    RNG() {}

    double rand() {
        // [0,1)の(double)numberを返す
        // min-max normalization(正規化)
        return double(xorshift() - min_) / (max_ - min_);
    }

    double randrange(const int& l, const int& r) {
        // [l,r)の(double)numberを返す
        return l + rand() * (r - l);
    }

private:
    std::uint32_t xorshift() {
        // XorShift Algorithm
        // [0,??]の整数を返す
        std::uint32_t t = (x ^ (x << 11));
        x = y;
        y = z;
        z = w;
        return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    }
};

// 実行時間管理
struct Timer {
public:
    // const std::int64_t TIME_LIMIT = 5000;  // DEBUG:
    const std::int64_t TIME_LIMIT = 1950;  // [msec]
private:
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
            return -1;  // office
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
        set_order(index, order);
        return;
    }

    void swap(const int& i, const int& j) {
        // tour[i]とtour[j]をswapする
        std::swap(tour[i], tour[j]);
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
            assert(false);
        }
        return;
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
    Timer timer = Timer();  // 実行時間管理
    Data data = Data();     // data set
    Job best_job = Job();   // 暫定解(今まで調べた中で一番良かったもの．最適解とは限らない)
    Job cur_job = Job();    // 仮の解
    RNG rng = RNG();        // 乱数生成器; random number generator

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
        // 初期解の作成 - 最近近傍法(from/to混合)
        // 後の処理のため，一旦適当にtourを埋める
        int tour_index = 1;
        for (const int& id : best_job.orders_id) {
            best_job.set_order(tour_index, Order(id, Dest::from));
            best_job.set_order(tour_index + 1, Order(id, Dest::to));
            tour_index += 2;
        }
        std::set<int> from_done_id;
        // tour[i] (1<=i<=TOUR_LEN-2) に最も近いtour[j] (2<=j<=TOUR_LEN-1, i<j) を見つけて，tour[i+1]=tour[j]とする
        for (int i = 1; i < data.TOUR_LEN - 2; i++) {
            Order order_i = best_job.tour[i];
            int min_dist = std::numeric_limits<int>::max();
            int min_dist_j = i + 1;
            for (int j = i + 1; j < data.TOUR_LEN - 1; j++) {
                Order order_j = best_job.tour[j];
                if ((order_j.dest == Dest::to) &&
                    (from_done_id.find(order_j.id) == from_done_id.end())) {
                    continue;
                }
                int dist_ = get_dist_by_point(
                    get_point_by_order(order_i),
                    get_point_by_order(order_j));
                if (dist_ < min_dist) {
                    min_dist = dist_;
                    min_dist_j = j;
                }
            }
            best_job.swap(i + 1, min_dist_j);
            if (best_job.tour[i + 1].dest == Dest::from) {
                from_done_id.insert(best_job.tour[i + 1].id);
            }
        }
        best_job.obj = calc_obj(best_job);
        print();  // DEBUG:
        return;
    }

    void simulated_annealing() {
        // 焼きなまし法
        cur_job.set_copy(best_job);  // init
        // 温度;temperature パラメータの設定
        // insertとswap searchを一回ずつ回してdeltaのminとaveを求めたのをハードコーディングしてる
        double min_temp = 1.0, max_temp = 50.0;
        while (timer.check_time_limit()) {
            // 温度の更新 (幾何冷却スケジューリング)
            // time_ = [0,1]
            double time_ = double(timer.get_elapsed_time()) / (timer.TIME_LIMIT);
            // temp = [min_temp,max_temp]
            double temp = pow(max_temp, 1.0 - time_) * pow(min_temp, time_);
            /*
            近傍探索 (neighborhood search)
            近傍を「insert処理またはswap処理を実施した状態」と定義する
            insert処理，swap処理共に，順序関係が壊れない場合にのみ実施する
            ただし，スコアが改悪した場合でも一定確率で採用する(焼きなましによって多様性を出したい)
            */
            // TODO: insertとswapの比率どうしよう？これも温度によって変えても良いかも
            // FIXME: nbhdらへんでバグってる(AtCoderに投げると20caseくらいWAが出た)．多分どこかの添え字ミス
            if (rng.rand() < 0.50) {
                insert_nbhd_search(temp);
            } else {
                swap_nbhd_search(temp);
            }
            if (cur_job.obj < best_job.obj) {
                best_job.set_copy(cur_job);
                print();  // DEBUG:
            } else {
                cur_job.set_copy(best_job);  // reset
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

    double get_current_score() {
        return (1e8 / (1000.0 + calc_obj(best_job)));
    }

    void insert_nbhd_search(const double& temp) {
        // insert近傍探索 (insert neighborhood search)
        // insert処理 <=def=> tour[i]をtour[j+1]の場所にinsertする．(j<i)

        int insert_cnt = 0;  // DEBUG:処理回数のメモ
        const int T = 100;   // insertの試行回数
        for (int t = 0; t < T; t++) {
            const int j = rng.randrange(2, data.TOUR_LEN - 3);
            const int i = rng.randrange(j + 2, data.TOUR_LEN - 3);
            Order order_i = cur_job.tour[i];
            Order order_j = cur_job.tour[j];
            if ((order_i.dest == Dest::to) &&
                (j < cur_job.get_index_by_order(Order(order_i.id, Dest::from)))) {
                // insertして順序関係を満たさなくなるならダメ
                continue;
            }
            Point j_point = get_point_by_order(order_j);
            Point next_j_point = get_point_by_order(cur_job.tour[j + 1]);
            Point prev_i_point = get_point_by_order(cur_job.tour[i - 1]);
            Point i_point = get_point_by_order(order_i);
            Point next_i_point = get_point_by_order(cur_job.tour[i + 1]);
            const int cur_dist = get_dist_by_point(j_point, next_j_point) +
                                 get_dist_by_point(prev_i_point, i_point) +
                                 get_dist_by_point(i_point, next_i_point);
            const int new_dist = get_dist_by_point(j_point, i_point) +
                                 get_dist_by_point(i_point, next_j_point) +
                                 get_dist_by_point(prev_i_point, next_i_point);
            const int delta = new_dist - cur_dist;
            if ((delta < 0) || (rng.rand() < exp(-delta / temp))) {
                // スコアが改善する or 一定確率 でこの近傍に移動
                cur_job.insert(j + 1, order_i);
                // cur_job.obj += delta;
                insert_cnt++;
            }
        }
        cur_job.obj = calc_obj(cur_job);
        if (LOCAL) {
            cout << "temp:" << temp << ",insert cnt:" << insert_cnt << endl;
            cout << "\tbest/cur obj:" << best_job.obj << "/" << cur_job.obj << endl;
        }
        return;
    }

    void swap_nbhd_search(const double& temp) {
        // swap近傍探索 (swap neighborhood search)
        // swap処理 <=def=> tour[i]とtour[j]をswapする

        // tabu list (将棋の千日手的なやつを防ぐ)
        // swapした(j,i)をtabuリストに追加する
        // swapする前に(i,j)がtabuリストにあれば処理をskip
        std::set<std::pair<int, int>> tabu;

        int swap_i_cnt = 0, swap_ij_cnt = 0;  // DEBUG:処理回数のメモ
        const int T = 100;                    // swapの試行回数
        for (int t = 0; t < T; t++) {
            {  // 以下tour[i]とtour[i+1]のswap
                bool conditions_passed = true;
                const int i = rng.randrange(2, data.TOUR_LEN - 2);
                if (tabu.find(std::make_pair(i, i + 1)) != tabu.end()) {
                    // すでに(i,i+1)のswapは確認済み
                    conditions_passed = false;
                }
                if ((cur_job.tour[i].id == cur_job.tour[i + 1].id) &&
                    (cur_job.tour[i].dest == Dest::from) && (cur_job.tour[i + 1].dest == Dest::to)) {
                    // swapしてi.fromとi.toが順序関係を満たさなくなるならダメ
                    conditions_passed = false;
                }
                Point prev_i_point = get_point_by_order(cur_job.tour[i - 1]);
                Point i_point = get_point_by_order(cur_job.tour[i]);
                Point next_i_point = get_point_by_order(cur_job.tour[i + 1]);
                Point next_next_i_point = get_point_by_order(cur_job.tour[i + 2]);
                const int cur_dist = get_dist_by_point(prev_i_point, i_point) +
                                     get_dist_by_point(i_point, next_i_point) +
                                     get_dist_by_point(next_i_point, next_next_i_point);
                const int new_dist = get_dist_by_point(prev_i_point, next_i_point) +
                                     get_dist_by_point(next_i_point, i_point) +
                                     get_dist_by_point(i_point, next_next_i_point);
                const int delta = new_dist - cur_dist;
                if (conditions_passed && ((delta < 0) || (rng.rand() < exp(-delta / temp)))) {
                    // スコアが改善する or 一定確率 でこの近傍に移動
                    cur_job.swap(i, i + 1);
                    tabu.insert(std::make_pair(i, i + 1));
                    // cur_job.obj += delta;
                    swap_i_cnt++;
                }
            }
            {  // 以下tour[i]とtour[j]のswap (i<j)
                const int i = rng.randrange(2, data.TOUR_LEN - 2);
                const int j = rng.randrange(i + 2, data.TOUR_LEN - 2);
                if (tabu.find(std::make_pair(i, j)) != tabu.end()) {
                    // すでに(i,j)のswapは確認済み
                    continue;
                }
                Order order_i = cur_job.tour[i];
                Order order_j = cur_job.tour[j];
                if ((order_i.dest == Dest::from) &&
                    (cur_job.get_index_by_order(Order(order_i.id, Dest::to)) < j)) {
                    // iは後ろに移動する
                    // swapしてi.fromとi.toが順序関係を満たさなくなるならダメ
                    continue;
                }
                if ((order_j.dest == Dest::to) &&
                    (i < cur_job.get_index_by_order(Order(order_j.id, Dest::from)))) {
                    // jは前に移動する
                    // swapしてj.fromとj.toが順序関係を満たさなくなるならダメ
                    continue;
                }
                Point prev_i_point = get_point_by_order(cur_job.tour[i - 1]);
                Point i_point = get_point_by_order(order_i);
                Point next_i_point = get_point_by_order(cur_job.tour[i + 1]);
                Point prev_j_point = get_point_by_order(cur_job.tour[j - 1]);
                Point j_point = get_point_by_order(order_j);
                Point next_j_point = get_point_by_order(cur_job.tour[j + 1]);
                const int cur_dist = get_dist_by_point(prev_i_point, i_point) +
                                     get_dist_by_point(i_point, next_i_point) +
                                     get_dist_by_point(prev_j_point, j_point) +
                                     get_dist_by_point(j_point, next_j_point);
                const int new_dist = get_dist_by_point(prev_i_point, j_point) +
                                     get_dist_by_point(j_point, next_i_point) +
                                     get_dist_by_point(prev_j_point, i_point) +
                                     get_dist_by_point(i_point, next_j_point);
                const int delta = new_dist - cur_dist;
                if ((delta < 0) || (rng.rand() < exp(-delta / temp))) {
                    // スコアが改善する or 一定確率 でこの近傍に移動
                    cur_job.swap(i, j);
                    tabu.insert(std::make_pair(i, j));
                    // cur_job.obj += delta;
                    swap_ij_cnt++;
                }
            }
        }
        cur_job.obj = calc_obj(cur_job);
        if (LOCAL) {
            cout << "temp:" << temp << ",swap cnt:" << swap_i_cnt << "," << swap_ij_cnt << endl;
            cout << "\tbest/cur obj:" << best_job.obj << "/" << cur_job.obj << endl;
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
    solver.simulated_annealing();
    solver.print();
    return 0;
}
