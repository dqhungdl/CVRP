#include <bits/stdc++.h>

#include "TSP.hpp"
#include "ls.hpp"
#include "utilities.hpp"

using namespace std;

const double INF = 1e15;
const double eps = 1e-6;
const int MAXN = 1000 + 10;

int demand[MAXN];
problem instance;
int n_trucks, truck_capacity;
long my_seed = 1;
double init_cost;

bool check_tour_capacity(long int* tour, int n) {
    int s = 0;
    for(int i = 0; i < n; ++i) s += demand[tour[i]];
    return s <= truck_capacity;
}

double tour_cost(long int* tour, int n) {
    double res = 0.0;
    for(int i = 0; i < n; ++i) res += instance.distance[tour[i]][tour[i + 1]];
    return res;
}

template<typename T>
T tour_cost(long int* tour, int n, vector< vector<T> > &distance) {
    T res = 0;
    for(int i = 0; i < n; ++i) res += distance[tour[i]][tour[i + 1]];
    return res;
}

template<typename T>
T tour_cost(vector<int> tour, int n, vector< vector<T> > &distance) {
    double res = 0.0;
    for(int i = 0; i < n; ++i) res += distance[tour[i]][tour[i + 1]];
    return res;
}

problem extract_problem(long int* tour, int n, problem &instance) {
    problem res;
    res.n = n; res.n_near = n - 1; nn_ls = n - 1;
    res.distance = vector< vector<long int> > (n, vector<long int>(n));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j) {
            res.distance[i][j] = instance.distance[tour[i]][tour[j]];
        }
    res.nn_list = vector< vector<long int> > (n, vector<long int>(n - 1));
    for(int i = 0; i < n; ++i) {
        vector< pair<long int, int> > nei;
        for(int j = 0; j < n; ++j)
            if (j != i) nei.push_back(make_pair(res.distance[i][j], j));
        sort(nei.begin(), nei.end());
        for(int j = 0; j < n - 1; ++j)
            res.nn_list[i][j] = nei[j].second;
    }
    return res;
}

//tour includes depot (0)
template<typename T>
double ls(long int* tour, int n, vector< vector<T> > &distance) {
    //cout << "ls begin\n";
    tour[n] = tour[0];
    //return tour_cost(tour, n);

    double init_cost = tour_cost(tour, n, distance);

    //for(int i = 0; i < n; ++i) cout << tour[i] << " "; cout << "\n";
    long int* new_tour = (long int*)malloc((n + 1) * sizeof(long int));
    for(int i = 0; i < n; ++i) new_tour[i] = i; new_tour[n] = new_tour[0];
    //for(int i = 0; i <= n; ++i) cout << new_tour[i] << " "; cout << "\n";

    problem new_instance = extract_problem(tour, n, instance);

    vector< vector<T> > extract_distance = vector< vector<T> > (n, vector<T>(n));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            extract_distance[i][j] = distance[tour[i]][tour[j]];
    //cout << "call ls()\n";
    two_opt_first<T>(new_tour, n, new_instance, extract_distance, &my_seed);
    //cout << "call done\n";
    assert(new_tour[0] == new_tour[n]);

    //convert
    vector<int> init_tour(n);
    for(int i = 0; i < n; ++i) init_tour[i] = tour[i];
    for(int i = 0; i < n; ++i) tour[i] = init_tour[new_tour[i]];
    free(new_tour);

    //for(int i = 0; i < n; ++i) cout << tour[i] << " "; cout << "\n";
    for(int i = 0; i < n; ++i) init_tour[i] = tour[i];
    for(int i = 0; i < n; ++i)
        if (tour[i] == 0) {
            for(int j = i; j < n; ++j)
                tour[j - i] = init_tour[j];
            for(int j = 0; j < i; ++j)
                tour[n - 1 - (i - 1) + j] = init_tour[j];
            break;
        }
    //for(int i = 0; i < n; ++i) cout << tour[i] << " "; cout << "\n";
    assert(tour[0] == 0);
    tour[n] = 0;

    //cout << "ls done\n";
    T new_cost = tour_cost(tour, n, distance);
    assert((new_cost <= init_cost) || (abs(new_cost - init_cost) <= 1e-9));
    return new_cost;
}

template<typename T>
double ls_vector(vector<int> &a, vector< vector<T> > &distance) {
    assert(a[0] == 0); assert(a.back() == 0);
    long int* tour = (long int*)malloc(a.size() * sizeof(long int));
    for(int i = 0; i < a.size(); ++i) tour[i] = a[i];
    double res = ls(tour, a.size() - 1, distance);
    for(int i = 0; i < a.size(); ++i) a[i] = tour[i];
    free(tour);
    return res;
}

//tour includes depot (0)
double ls_1tour(long int* tour, int n, vector<int> &L) {
    assert(check_tour_capacity(tour, n) == true);
    double res = ls(tour, n, instance.distance);
    L.clear();
    for(int i = 1; i < n; ++i) L.push_back(tour[i]);
    return res;
}

double ls_split(long int* tour, int n, bool include_depot, int &best_split) {
    //cout << "ls_split() begin...\n";
    //for(int i = 0; i < n; ++i) cout << tour[i] << " "; cout << "\n";
    double res = INF;

    int tot_demand = 0;
    for(int i = 0; i < n; ++i) tot_demand += demand[tour[i]];

    int s = 0;
    long int* L = (long int*)malloc((n + 2) * sizeof(long int));
    long int* R = (long int*)malloc((n + 2) * sizeof(long int));

    int start = include_depot ? 1 : 0;
    int valid_split_cnt = 0;
    for(int split_point = start; split_point + 1 < n; ++split_point) {
        s += demand[tour[split_point]];
        if ((s <= truck_capacity) && (tot_demand - s <= truck_capacity)) {
            ++valid_split_cnt;

            //cout << "yeah\n";
            L[0] = 0;
            for(int i = start; i <= split_point; ++i) L[i - start + 1] = tour[i];
            R[0] = 0;
            for(int i = split_point + 1; i < n; ++i) R[i - (split_point + 1) + 1] = tour[i];
            //cout << "do it\n";
            //for(int i = 0; i < split_point - start + 1 + 1; ++i) cout << L[i] << " "; cout << "\n";
            double ls_L = ls(L, split_point - start + 1 + 1, instance.distance);
            //for(int i = 0; i < n - (split_point - start + 1) + 1; ++i) cout << R[i] << " "; cout << "\n";
            double ls_R = ls(R, n - (start == 1) - (split_point - start + 1) + 1, instance.distance);
            if (res > ls_L + ls_R) {
                res = ls_L + ls_R;
                best_split = split_point;
            }
            //cout << "ls() done\n";
        }
    }

    free(L); free(R);
    //cout << "n valid split = " << valid_split_cnt << " / " << n - 1 - start + 1 << "\n";
    //cout << "ls_split() done\n";
    return res;
}

void move_element(long int* tour, int n, int from, int to) {
    if (from == to) return;

    int e = tour[from];
    int dir = (from < to) ? 1 : -1;
    for(int i = from; i != to; i += dir) tour[i] = tour[i + dir];
    tour[to] = e;
}

//tour includes depot (0)
double ls_2tour(long int* tour, int n, vector<int> &L, vector<int> &R) {
    bool improved = true;
    int best_split = -1;
    //cout << "dd\n";
    //for(int i = 0; i < n; ++i) cout << tour[i] << " "; cout << "\n";
    double best_cost = ls_split(tour, n, true, best_split);
    //cout << "ee\n";

    while (improved) {
        improved = false;
        //cout << "in while\n";
        //swap
        double best_swap_cost = INF;
        int best_swap_i = -1, best_swap_j = -1;
        int best_swap_split = -1;
        bool ok = false;

        for(int i = 1; i < n; ++i) {
            for(int j = i + 1; j < n; ++j) {
                swap(tour[i], tour[j]);
                int split = -1;
                //cout << "11\n";
                double swap_cost = ls_split(tour, n, true, split);
                //cout << "22\n";
                swap(tour[i], tour[j]);
                if (best_swap_cost > swap_cost) {
                    best_swap_cost = swap_cost;
                    best_swap_i = i; best_swap_j = j;
                    best_swap_split = split;

                    if (best_swap_cost < best_cost) {
                        ok = true;
                    }
                }
            }
        }

        if ((abs(best_cost - best_swap_cost) > eps) && (best_swap_cost < best_cost)) {
            improved = true;
            best_cost = best_swap_cost;
            swap(tour[best_swap_i], tour[best_swap_j]);
            best_split = best_swap_split;
        }

        //move
        double best_move_cost = INF;
        int best_move_from = -1, best_move_to = -1;
        int best_move_split = -1;

        for(int i = 1; i < n; ++i) {
            for(int j = 1; j < n; ++j)
                if (i != j) {
                    move_element(tour, n, i, j);
                    int split = -1;
                    //cout << "1111\n";
                    double move_cost = ls_split(tour, n, true, split);
                    //cout << "2222\n";
                    move_element(tour, n, j, i);
                    if (move_cost < best_move_cost) {
                        best_move_cost = move_cost;
                        best_move_from = i; best_move_to = j;
                        best_move_split = split;
                    }
                }
        }

        if ((abs(best_cost - best_move_cost) > eps) && (best_move_cost < best_cost)) {
            improved = true;
            best_cost = best_move_cost;
            move_element(tour, n, best_move_from, best_move_to);
            best_split = best_move_split;
        }

    }

    //cout << "??\n";
    L.clear(); R.clear();
    if (best_split >= 0) {
        for(int i = 1; i <= best_split; ++i) L.push_back(tour[i]);
        for(int i = best_split + 1; i < n; ++i) R.push_back(tour[i]);
    }

    return best_cost;
}

vector<int> min_trucks(vector<int> &per) {
    vector<int> res;
    for(int i = 0; i < per.size(); ++i) {
        vector<int> a;
        for(int j = i; j < per.size(); ++j)
            a.push_back(demand[ per[j] ]);

        sort(a.begin(), a.end(), greater<int>());

        int x = 0;
        vector<bool> used(a.size(), false);
        while (true) {
            int s = 0;
            for(int j = 0; j < a.size(); ++j)
                if ((!used[j]) && (s + a[j] <= truck_capacity)) {
                    s += a[j];
                    used[j] = true;
                }
            if (s == 0) break;
            x++;
        }

        res.push_back(x);
    }
    res.push_back(0);
    return res;
}

void dp(vector<int> &per) {
    int n = per.size();
    assert(n == instance.n - 1);

    vector< vector<double> > f(n + 1, vector<double> (n_trucks + 1, INF));
    vector< vector<double> > g1(n + 1, vector<double> (n + 1, -1.0));
    vector< vector< vector<int> > > g1_L(n + 1, vector< vector<int> > (n + 1));
    vector< vector<double> > g2(n + 1, vector<double> (n + 1, -1.0));
    vector< vector< vector<int> > > g2_L(n + 1, vector< vector<int> > (n + 1));
    vector< vector< vector<int> > > g2_R(n + 1, vector< vector<int> > (n + 1));

    vector< vector<int> > trace(n + 1, vector<int> (n_trucks + 1));
    vector< vector< vector< vector<int> > > > trace_tours(n + 1, vector< vector< vector<int> > > (n_trucks + 1));
    f[0][0] = 0.0;

    long int* tour = (long int*)malloc((n + 10) * sizeof(long int));
    vector<int> L, R;
    vector<int> mt = min_trucks(per);
    for(int k = 1; k <= n_trucks; ++k) {
        cout << "k = " << k << " / " << n_trucks << "\n";
        for(int i = k; i <= n; ++i) {
            //cout << "k = " << k << ", i = " << i << "\n";
            if (mt[i] > n_trucks - k) continue;

            trace_tours[i][k].clear();

            int s = 0;
            for(int j = i; j >= max(1, k - 1); --j) {
                s += demand[per[j - 1]];
                if (s > truck_capacity * 2) break;

                //cout << "k = " << k << ", i = " << i << ", j = " << j << "\n";

                //cout << "j-i 1 tour\n";

                if ((j >= k) && (s <= truck_capacity) && (f[j - 1][k - 1] + 10 < INF)) {
                    if (g1[j][i] < 0) {
                        tour[0] = 0;
                        for(int p = j - 1; p <= i - 1; ++p) tour[p - (j - 1) + 1] = per[p];
                        g1[j][i] = ls_1tour(tour, i - j + 1 + 1, g1_L[j][i]);
                    }

                    double ls_1 = g1[j][i];
                    if (f[i][k] > f[j - 1][k - 1] + ls_1) {
                        trace[i][k] = j;
                        trace_tours[i][k].clear();
                        trace_tours[i][k].push_back(g1_L[j][i]);
                        f[i][k] = f[j - 1][k - 1] + ls_1;
                    }
                }


                //cout << "j-i 2 tour\n";
                //tour may have been modified

                if ((i - j + 1 >= 2) && (k >= 2) && (f[j - 1][k - 2] + 10.0 < INF)) {
                    if (g2[j][i] < 0) {
                        tour[0] = 0;
                        for(int p = j - 1; p <= i - 1; ++p) tour[p - (j - 1) + 1] = per[p];
                        g2[j][i] = ls_2tour(tour, i - j + 1 + 1, g2_L[j][i], g2_R[j][i]);
                    }
                    double ls_2 = g2[j][i];
                    if (f[i][k] > f[j - 1][k - 2] + ls_2) {
                        trace[i][k] = j;
                        trace_tours[i][k].clear();
                        trace_tours[i][k].push_back(g2_L[j][i]);
                        trace_tours[i][k].push_back(g2_R[j][i]);
                        f[i][k] = f[j - 1][k - 2] + ls_2;
                    }
                }
                //cout << "2 tour done\n";
            }

            //if (f[i][k] + 10 < INF) cout << "f(i = " << i << ", k = " << k << ") = " << f[i][k] << "\n";
        }
    }

    free(tour);

    double res = INF;
    int j = -1;
    for(int k = 1; k <= n_trucks; ++k)
        if (res > f[n][k]) {
            j = k;
            res = f[n][k];
        }

    cout << "best = " << res << "\n";
    int i = n;
    vector< vector<int> > sol;
    while (i > 0) {
        cout << i << " " << j << " " << f[i][j] << " " << " " << trace[i][j] << " " << trace_tours[i][j].size() << "\n";
        assert(trace_tours[i][j].size() > 0);
        for(vector<int> &t : trace_tours[i][j])
            sol.push_back(t);
        int ii = i;
        i = trace[ii][j] - 1;
        j -= trace_tours[ii][j].size();
    }
    assert(j == 0);

    double sc = 0;
    vector<bool> used(instance.n + 1, false);
    int nn = 0;
    for(vector<int> &t : sol) {
        cout << "0 "; for(int v : t) cout << v << " "; cout << "\n";
        for(int v : t) {
            if (!used[v]) {
                used[v] = true;
                nn++;
            }
            else {
                assert(false);
            }
        }
        int d = 0;
        for(int v : t) d += demand[v];
        assert(d <= truck_capacity);

        double c = instance.distance[0][t[0]] + instance.distance[t[(int)(t.size()) - 1]][0];
        for(int i = 1; i < t.size(); ++i) c += instance.distance[t[i - 1]][t[i]];
        cout << "cost = " << c << "\n";
        sc += c;
    }
    cout << "init cost = " << init_cost << "\n";
    cout << "total cost = " << sc << "\n";
    cout << "best = " << res << "\n";
    cout << "nn = " << n << " / " << instance.n << "\n";
}

//round distance
long int distance(point &a, point &b) {
    double xd = a.x - b.x;
    double yd = a.y - b.y;
    double r  = sqrt(xd*xd + yd*yd) + 0.5;
    return (long int) r;
}
double distance2(point &a, point &b) {
    double xd = a.x - b.x;
    double yd = a.y - b.y;
    double r  = sqrt(xd*xd + yd*yd);
    return r;
}

string inp_file, sol_file, log_file, instance_file;
int kind = 2;

void read_input() {
    ifstream f(inp_file);
    string tmp;
    getline(f, tmp);
    getline(f, tmp);
    getline(f, tmp);
    f >> tmp >> tmp;
    int n; f >> n;
    getline(f, tmp);
    getline(f, tmp);
    f >> tmp >> tmp;
    f >> truck_capacity;
    getline(f, tmp);
    getline(f, tmp);

    cout << "n = " << n << ", capacity = " << truck_capacity << "\n";
    instance.n = n; instance.n_near = n - 1;
    instance.distance = vector< vector<long int> > (n, vector<long int>(n, 0));
    instance.distance2 = vector< vector<double> > (n, vector<double>(n, 0));
    instance.nn_list = vector< vector<long int> > (n, vector<long int>(n - 1, 0));
    instance.points = vector<point>(n);

    for(int i = 0; i < n; ++i) {
        int id;
        f >> id >> instance.points[i].x >> instance.points[i].y;
        assert(id == i + 1);
    }

    f >> tmp;
    for(int i = 0; i < n; ++i) {
        int id;
        f >> id >> demand[i];
        assert(id == i + 1);
        if (i == 0) assert(demand[i] == 0);
    }

    f >> tmp;
    int depot; f >> depot;
    assert(depot == 1);
    int v; f >> v;
    assert(v == -1);

    f.close();

    //compute distance between points
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j) {
            instance.distance[i][j] = distance(instance.points[i], instance.points[j]);
        }

    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j) {
            instance.distance2[i][j] = distance2(instance.points[i], instance.points[j]);
        }

    //compute nn_list
    for(int i = 0; i < n; ++i) {
        vector< pair<long int, int> > nei;
        for(int j = 0; j < n; ++j)
            if (j != i) {
                nei.push_back(make_pair(instance.distance[i][j], j));
            }
        sort(nei.begin(), nei.end());
        for(int j = 0; j < n - 1; ++j)
            instance.nn_list[i][j] = nei[j].second;
    }
    nn_ls = n;
}

vector< vector<int> > tours;
vector<int> read_sol() {
    ifstream f(sol_file);
    long int v;
    f >> v;
    f >> n_trucks;
    f >> v;
    f >> v;
    cout << "n tours = " << n_trucks << "\n";

    tours;
    double tot_cost = 0.0;
    tours.clear();
    for(int i = 0; i < n_trucks; ++i) {
        f >> v >> v >> v >> v >> v;
        int s; f >> s;
        vector<int> tour(s, 0);
        for(int j = 0; j < s; ++j) {
            f >> tour[j];
        }
        tours.push_back(tour);
        double tc = tour_cost(tour, s - 1, instance.distance);

        int d = 0;
        for(int j = 0; j < s; ++j) d += demand[tour[j]];
        assert(d <= truck_capacity);

        for(int v : tour) cout << v << " "; cout << "\n";
        cout << "tour " << i + 1 << " length = " << tc << "\n";

        tot_cost += tc;
    }
    f.close();

    cout << "initial cost = " << tot_cost << "\n";
    init_cost = tot_cost;
    vector<int> per;
    for(auto &t : tours) {
        for(int i = 1; i + 1 < t.size(); ++i)
            per.push_back(t[i]);
    }

    vector<bool> used(instance.n + 1, false);
    for(int v : per) {
        if (used[v]) {
            cout << v << " appeared > 1 times\n";
            assert(false);
        }
        else {
            used[v] = true;
        }
    }
    cout << "per length = " << per.size() << " / " << instance.n - 1 << "\n";
    //sort(per.begin(), per.end()); for(int v : per) cout << v << " "; cout << "\n";
    return per;
}

vector<int> read_sol_2() {
    ifstream f(sol_file);
    long int v;
    n_trucks = 0;
    double tot_cost = 0.0;
    tours.clear();
    string tmp;
    f >> tmp;
    int Count=0;
    while (tmp == "Route") {
        Count++;
        n_trucks++;
        f >> tmp;
        vector<int> tour;
        tour.push_back(0);
        f >> tmp;
        while ((tmp != "Route") && (tmp != "Cost")) {
            int x = atoi(tmp.c_str());
            tour.push_back(x);
            f >> tmp;
        }
        tour.push_back(0);
        int s = tour.size();
        tours.push_back(tour);
        double tc = tour_cost(tour, s - 1, instance.distance);

        int d = 0;
        for(int j = 0; j < s; ++j) {
            d += demand[tour[j]];
        }
        assert(d <= truck_capacity);

        for(int v : tour) cout << v << " "; cout << "\n";
        cout << "tour " << n_trucks << " length = " << tc << "\n";

        tot_cost += tc;
    }
    assert(tmp == "Cost");
    int report_cost; f >> report_cost;
    f.close();

    cout << "initial cost = " << tot_cost << " \n";
    init_cost = tot_cost;
    vector<int> per;
    for(auto &t : tours) {
        for(int i = 1; i + 1 < t.size(); ++i)
            per.push_back(t[i]);
    }

    vector<bool> used(instance.n + 1, false);
    for(int v : per) {
        if (used[v]) {
            cout << v << " appeared > 1 times\n";
            assert(false);
        }
        else {
            used[v] = true;
        }
    }
    cout << "per length = " << per.size() << " / " << instance.n - 1 << "\n";
    //sort(per.begin(), per.end()); for(int v : per) cout << v << " "; cout << "\n";
    return per;
}

template<typename T>
bool change_2tour(vector<int> &a, vector<int> &b, int len_cut, vector< vector<T> > distance) {
    int m = a.size(), n = b.size();
    assert(a[0] == 0); assert(a[m - 1] == 0);
    assert(b[0] == 0); assert(b[n - 1] == 0);

    vector< vector< vector<double> > > f(m + 1, vector< vector<double> >(n + 1, vector<double>(truck_capacity + 1, INF)));
    vector< vector< vector< pair<int, int> > > > trace(m + 1, vector< vector< pair<int, int> > >(n + 1, vector< pair<int, int> >(truck_capacity + 1)));
    vector< vector< vector<int> > > trace_c1(m + 1, vector< vector<int> >(n + 1, vector<int>(truck_capacity + 1, -1)));

    f[0][0][0] = 0.0;
    int sda = 0;
    for(int i = 0; i < m; ++i) {
        sda += demand[a[i]];
        int sdb = 0;
        for(int j = 0; j < n; ++j) {
            sdb += demand[b[j]];
            for(int c1 = 0; c1 <= truck_capacity; ++c1) {
                int c2 = sda + sdb - c1;
                if (!(f[i][j][c1] + 10.0 < INF)) continue;
                assert(c2 <= truck_capacity);
                //cout << i << " " << j << " " << c1 << " " << f[i][j][c1] << "\n";

                //i -> i + 1
                if ((i + 1 < m) && (c1 + demand[a[i + 1]] <= truck_capacity)) {
                    if (f[i + 1][j][c1 + demand[a[i + 1]]] > f[i][j][c1] + distance[a[i]][a[i + 1]]) {
                        f[i + 1][j][c1 + demand[a[i + 1]]] = f[i][j][c1] + distance[a[i]][a[i + 1]];
                        trace[i + 1][j][c1 + demand[a[i + 1]]] = make_pair(i, j);
                        trace_c1[i + 1][j][c1 + demand[a[i + 1]]] = c1;
                    }

                    double dist_cut = 0.0;
                    int next_c1 = c1 + demand[a[i + 1]];
                    for(int l = 1; l <= len_cut; ++l) {
                        if (j + l + 1 >= n) break;
                        next_c1 += demand[b[j + l]];
                        if (next_c1 > truck_capacity) break;
                        if (l > 1) dist_cut += distance[b[j + l - 1]][b[j + l]];

                        double add_cost = dist_cut + distance[a[i]][b[j + 1]] + distance[b[j + l]][a[i + 1]] + distance[b[j]][b[j + l + 1]];
                        if (c2 + demand[b[j + l + 1]] <= truck_capacity) {
                            if (f[i + 1][j + l + 1][next_c1] > f[i][j][c1] + add_cost) {
                                f[i + 1][j + l + 1][next_c1] = f[i][j][c1] + add_cost;
                                trace[i + 1][j + l + 1][next_c1] = make_pair(i, j);
                                trace_c1[i + 1][j + l + 1][next_c1] = c1;
                            }
                        }
                    }
                }

                //j -> j + 1
                if ((j + 1 < n) && (c2 + demand[b[j + 1]] <= truck_capacity)) {
                    if (f[i][j + 1][c1] > f[i][j][c1] + distance[b[j]][b[j + 1]]) {
                        f[i][j + 1][c1] = f[i][j][c1] + distance[b[j]][b[j + 1]];
                        trace[i][j + 1][c1] = make_pair(i, j);
                        trace_c1[i][j + 1][c1] = c1;
                    }

                    double dist_cut = 0.0;
                    int next_c2 = c2 + demand[b[j + 1]];
                    for(int l = 1; l <= len_cut; ++l) {
                        if (i + l + 1 >= m) break;
                        next_c2 += demand[a[i + l]];
                        if (next_c2 > truck_capacity) break;
                        if (l > 1) dist_cut += distance[a[i + l - 1]][a[i + l]];

                        double add_cost = dist_cut + distance[b[j]][a[i + 1]] + distance[a[i + l]][b[j + 1]] + distance[a[i]][a[i + l + 1]];
                        if (c1 + demand[a[i + l + 1]] <= truck_capacity) {
                            if (f[i + l + 1][j + 1][c1 + demand[a[i + l + 1]]] > f[i][j][c1] + add_cost) {
                                f[i + l + 1][j + 1][c1 + demand[a[i + l + 1]]] = f[i][j][c1] + add_cost;
                                trace[i + l + 1][j + 1][c1 + demand[a[i + l + 1]]] = make_pair(i, j);
                                trace_c1[i + l + 1][j + 1][c1 + demand[a[i + l + 1]]] = c1;
                            }
                        }
                    }
                }

                //i -> j + 1, j->i + 1
                if ((i + 1 < m) && (j + 1 < n))
                    if ((c1 + demand[b[j + 1]] <= truck_capacity) && (c2 + demand[a[i + 1]] <= truck_capacity)) {
                        if (f[i + 1][j + 1][c2 + demand[a[i + 1]]] > f[i][j][c1] + distance[a[i]][b[j + 1]] + distance[b[j]][a[i + 1]]) {
                            f[i + 1][j + 1][c2 + demand[a[i + 1]]] = f[i][j][c1] + distance[a[i]][b[j + 1]] + distance[b[j]][a[i + 1]];
                            trace[i + 1][j + 1][c2 + demand[a[i + 1]]] = make_pair(i, j);
                            trace_c1[i + 1][j + 1][c2 + demand[a[i + 1]]] = c1;
                        }
                    }
            }
        }
    }

    double res = INF;
    int best_c1 = -1;
    for(int c1 = 0; c1 <= truck_capacity; c1++)
        if (f[m - 1][n - 1][c1] < res) {
            res = f[m - 1][n - 1][c1];
            best_c1 = c1;
        }

    cout << "--------------------------------------------------------\n";
    T init_cost = tour_cost(a, m - 1, distance) + tour_cost(b, n - 1, distance);
    cout << "init cost = " << init_cost << "\n";
    cout << "res dp = " << res << "\n";

    cout << m << ": "; for(int v : a) cout << v << " "; cout << "\n";
    cout << n << ": "; for(int v : b) cout << v << " "; cout << "\n";

    double eps = 1e-9;
    if ((abs(res - init_cost) > eps) && (res > init_cost)) {
        cout << "dp absolutely wrongggggg!!\n";
        assert(false);
    }

    //trace and check
    int i = m - 1, j = n - 1, c1 = best_c1;
    vector< vector<int> > res_tour(2, vector<int>(1, 0));
    int ti = 0, tj = 1;
    while ((i > 0) || (j > 0)) {
        cout << i << " " << j << "\n";
        int ni = trace[i][j][c1].first, nj = trace[i][j][c1].second;
        c1 = trace_c1[i][j][c1];

        if (i == ni) {
            assert(nj == j - 1);
            res_tour[tj].push_back(b[j - 1]);
        }
        else if (j == nj) {
            assert(ni == i - 1);
            res_tour[ti].push_back(a[i - 1]);
        }
        else if (ni == i - 1) {
            if (nj == j - 1) {
                res_tour[tj].push_back(a[i - 1]);
                res_tour[ti].push_back(b[j - 1]);
                swap(ti, tj);
            }
            else {
                for(int k = j - 1; k > nj; --k)
                    res_tour[ti].push_back(b[k]);
                res_tour[ti].push_back(a[i - 1]);
                res_tour[tj].push_back(b[nj]);
            }
        }
        else {
            assert(nj == j - 1);
            for(int k = i - 1; k > ni; --k)
                res_tour[tj].push_back(a[k]);
            res_tour[tj].push_back(b[j - 1]);
            res_tour[ti].push_back(a[ni]);
        }
        i = ni; j = nj;
    }

    reverse(res_tour[0].begin(), res_tour[0].end());
    reverse(res_tour[1].begin(), res_tour[1].end());

    for(int v : res_tour[0]) cout << v << " "; cout << "\n";
    for(int v : res_tour[1]) cout << v << " "; cout << "\n";

    assert(res_tour[0][0] == 0); assert(res_tour[0].back() == 0);
    assert(res_tour[1][0] == 0); assert(res_tour[1].back() == 0);

    set<int> s;
    for(int k = 1; k + 1 < a.size(); ++k) {
        assert(s.count(a[k]) == 0);
        s.insert(a[k]);
    }
    for(int k = 1; k + 1 < b.size(); ++k) {
        assert(s.count(b[k]) == 0);
        s.insert(b[k]);
    }
    for(int k = 1; k + 1 < res_tour[0].size(); ++k) {
        assert(s.count(res_tour[0][k]) == 1);
        s.erase(res_tour[0][k]);
    }
    for(int k = 1; k + 1 < res_tour[1].size(); ++k) {
        assert(s.count(res_tour[1][k]) == 1);
        s.erase(res_tour[1][k]);
    }
    assert(s.size() == 0);

    double check_cost = tour_cost(res_tour[0], res_tour[0].size() - 1, distance) + tour_cost(res_tour[1], res_tour[1].size() - 1, distance);
    assert(abs(check_cost - res) <= eps);
    int cap_0 = 0;
    for(int v : res_tour[0]) cap_0 += demand[v];
    assert(cap_0 <= truck_capacity);
    int cap_1 = 0;
    for(int v : res_tour[0]) cap_1 += demand[v];
    assert(cap_1 <= truck_capacity);
    assert((cap_0 == best_c1) || (cap_1 == best_c1));

    double ls_cost = ls_vector(res_tour[0], distance) + ls_vector(res_tour[1], distance);
    if ((ls_cost < res) && (abs(ls_cost - res) > eps)) {
        cout << "ls success!, reduced: " << res - ls_cost << "\n";
        res = ls_cost;
        cout << "new res = " << res << "\n";
        cout << "new tours after ls:\n";
        for(int v : res_tour[0]) cout << v << " "; cout << "\n";
        for(int v : res_tour[1]) cout << v << " "; cout << "\n";
    }

    if ((res < init_cost) && (abs(res - init_cost) > eps)) {
        cout << "found "<<res<<" "<<init_cost<<"\n";
        //assert(false);
        return true;
    }
    return false;
}

void parse_commandline(int argc, char *argv[]) {

     int i;

     for (i=1; i<argc; i=i+2) {
     
         // -i : ten file du lieu vao    
         if (strcmp(argv[i],"-i")==0) {

			instance_file = argv[i+1];
		}

		 // -k : ten file du lieu vao    
         if (strcmp(argv[i],"-k")==0) {

			kind = atoi(argv[i+1]);
		}

     } 
}

int main(int argc, char *argv[])
{
	//parse_commandline(argc, argv);
    
	
    inp_file = "input.txt";
    sol_file = "output.txt";
	log_file = "log.txt";
	cerr<<inp_file<<endl;
	cerr<<sol_file<<endl;
	cerr<<log_file<<endl;

    freopen(log_file.c_str(), "w", stdout);

    read_input();
    //if (kind == 2) vector<int> per = read_sol_2();
    //else vector<int> per = read_sol();
    vector<int> per = read_sol_2();

    //return 0;
    int n_tour = tours.size();
    //n_tour = 2;
    int n_round_dist_improved = 0, n_real_dist_improved = 0, n_both_dist_improved = 0, n_pair = 0;
    for(int i = 0; i < n_tour; ++i) {
        for(int j = i + 1; j < n_tour; ++j) {
            int round_dist = change_2tour(tours[i], tours[j], 3, instance.distance);
            int real_dist = change_2tour(tours[i], tours[j], 3, instance.distance2);
            n_pair++;
            n_round_dist_improved += round_dist;
            n_real_dist_improved += real_dist;
            n_both_dist_improved += ((round_dist == 1) && (real_dist == 1));
        }
    }
    for(int i=0;i<n_tour;i++)
        debug(tours[i]);
    cout << "done, no assertion failed\n";

    cout << "n_round_dist_improved = " << n_round_dist_improved << " / " << n_pair << "\n";
    cout << "n_real_dist_improved = " << n_real_dist_improved << " / " << n_pair << "\n";
    cout << "n_both_dist_improved = " << n_both_dist_improved << " / " << n_pair << "\n";

    cerr << "n_round_dist_improved = " << n_round_dist_improved << " / " << n_pair << "\n";
    cerr << "n_real_dist_improved = " << n_real_dist_improved << " / " << n_pair << "\n";
    cerr << "n_both_dist_improved = " << n_both_dist_improved << " / " << n_pair << "\n";
    PAUSE();
    //dp(per);
}
