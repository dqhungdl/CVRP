#ifndef ACO_HPP_INCLUDED
#define ACO_HPP_INCLUDED

#include <bits/stdc++.h>

#include "structs.hpp"
#include "Utilities.hpp"
#include "BipartiteMatching.hpp"
#include "LS.hpp"

class ACO {
public:
    //input:
    Problem* instance;

    //other variables:
    double INF = 1e15;
    int n, n_trucks, truck_capacity;
    double T_min, T_max; //trail min/max
    double alpha = 0, beta = 1.0, rho = 0.05;
    std::vector< std::vector<double> > trail;
    std::vector< std::vector<double> > heuristic;
    double penalty = 10000.0;
    std::vector< std::vector<int> > freq;

    //methods:
    void init(Problem* _instance);
    ACO();
    ACO(Problem* _instance);
    virtual ~ACO();

    void init_ACO();
    std::vector<int> choose_trucks(std::vector<Tour> &tours, std::vector<int> &available_trucks, int n_trucks);
    std::vector<int> visit_customers(std::vector<Tour*> truck_tours, std::vector<int> &customers);
    void fix_solution(Solution &sol);
    Solution ant_construct_solution();
    int best_solution(std::vector<Solution> &sols);
    std::vector<Solution> generate_solutions(int n_ants);
    void update_trails(Solution &sol);
    void update_trails(std::vector<Solution> &sols);
    void run();
};

void ACO::init(Problem* _instance) {
    instance = _instance;
    n = instance->n;
    n_trucks = instance->n_trucks;
    truck_capacity = instance->truck_capacity;

    freq = std::vector< std::vector<int> >(n, std::vector<int>(n, 0));
}

ACO::ACO() {
}

ACO::ACO(Problem* _instance) {
    init(_instance);
}

ACO::~ACO() {
}

void ACO::init_ACO() {
    T_max = 1.0; T_min = T_max / n;

    trail = std::vector< std::vector<double> >(n, std::vector<double>(n, T_min));
    heuristic = std::vector< std::vector<double> >(n, std::vector<double>(n, 0.0));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            if (i != j)
                heuristic[i][j] = instance->distance[i][j];
}

//side effect: changes 'available_trucks'
std::vector<int> ACO::choose_trucks(std::vector<Tour> &tours, std::vector<int> &available_trucks, int n_trucks) {
    n_trucks = std::min(n_trucks, (int)available_trucks.size());
    std::vector<int> truck_id_chosen(n_trucks);
    /*for(int i = 0; i < n_trucks; ++i) {
        int j = Utilities::rnd(0, available_trucks.size() - 1 - i);
        truck_id_chosen[i] = available_trucks[j];
        std::swap(available_trucks[j], available_trucks[available_trucks.size() - 1 - i]);
    }*/
    std::vector<double> w(available_trucks.size());
    for(int i = 0; i < w.size(); ++i) w[i] = 1.0 / (tours[ available_trucks[i] ].tot_demand + 1);
    for(int i = 0; i < n_trucks; ++i) {
        int j = Utilities::choice(w);
        truck_id_chosen[i] = available_trucks[j];
        std::swap(available_trucks[j], available_trucks[available_trucks.size() - 1 - i]);
        std::swap(w[j], w[w.size() - 1]);
        w.pop_back();
    }
    return truck_id_chosen;
}

std::vector<int> ACO::visit_customers(std::vector<Tour*> tours, std::vector<int> &customers) {
    int n_trucks = tours.size();
    int n_customers = customers.size();
    int n = std::max(n_trucks, n_customers);

    std::vector< std::vector<double> > cost(n, std::vector<double>(n, 1e8));
    for(int i = 0; i < n_trucks; ++i) {
        for(int j = 0; j < n_customers; ++j) {
            if (customers[j] >= 0) {
                //if (tours[i]->tot_demand + instance->demand[ customers[j] ] > truck_capacity) continue;
                double p = trail[ tours[i]->back() ][ customers[j] ];
                double h = heuristic[ tours[i]->back() ][ customers[j] ];
                cost[i][j] = pow(p, alpha) * pow(h, beta);
                if (tours[i]->tot_demand + instance->demand[ customers[j] ] > truck_capacity) cost[i][j] += penalty;

                //cost[i][j] = 1.0 / (cost[i][j] + 1);
            }
            else {
                assert(customers[j] == -1);
                cost[i][j] = 0.0; // virtual customer
            }
        }
    }

    /*for(int i = 0; i < n_trucks; ++i) {
        for(int j = 0; j < n_customers; ++j) {
            std::cout << cost[i][j] << " ";
        }
        std::cout << "\n";
    }*/

    std::vector<int> Lmate, Rmate;
    double matching = BipartiteMatching::min_cost_matching(cost, Lmate, Rmate);
    //std::cout << "Lmate = \n"; for(int i = 0; i < n_trucks; ++i) std::cout << Lmate[i] << " "; std::cout << "\n";
    //std::cout << "Rmate = \n"; for(int i = 0; i < n_customers; ++i) std::cout << Rmate[i] << " "; std::cout << "\n";

    std::vector<int> customers_visited;
    for(int i = 0; i < n_trucks; ++i) {
        /*std::vector< std::pair<int, int> > pp;
        for(int k = 0; k < n_customers; ++k) pp.push_back(std::make_pair(cost[i][k], k));
        sort(pp.begin(), pp.end());
        int top = std::min((int)pp.size(), 20);
        std::vector<double> w(top); for(int k = 0; k < top; ++k) w[k] = pp[k].first;
        int k = Utilities::choice(w);
        int j = pp[k].second;
        j = customers[j];*/
        if (Lmate[i] < n_customers) {
            int j = customers[ Lmate[i] ];
            //j = customers[ Utilities::choice(cost[i]) ];
            if (j < 0) continue;
            freq[ tours[i]->back() ][ j ]++;
            tours[i]->add(j, instance->demand[ j ]);
            customers_visited.push_back(j);
        }
    }
    return customers_visited;
}

void ACO::fix_solution(Solution &sol) {
    for(Tour &t : sol.tours) {
        LS::ls_tour(t, *instance);
    }

    std::cout << "begin dp...\n";
    int n_while = 500;
    while (true) {
        if (n_while <= 0) break;
        --n_while;

        bool improved = false;
        for(int i = 0; i < sol.tours.size(); ++i)
            for(int j = i + 1; j < sol.tours.size(); ++j) {
                bool ok = LS::exchange_2tour(sol.tours[i].order, sol.tours[j].order, 5, *instance);
                if (ok) {
                    sol.tours[i].calc_tot_demand(instance->demand);
                    sol.tours[j].calc_tot_demand(instance->demand);
                    improved = true;
                }
            }

        if (!improved) break;
    }
    std::cout << "dp end\n";
}

Solution ACO::ant_construct_solution() {
    std::cout << "single ant begin...\n";
    std::vector<Tour> tours(n_trucks);
    for(int i = 0; i < n_trucks; ++i) {
        tours[i].id = i; tours[i].add(0, instance->demand[0]);
    }

    std::vector<bool> full_truck(n_trucks, false);
    std::vector<bool> visited_customer(n, false);
    visited_customer[0] = true; //0 is depot

    while (true) {
        std::vector<int> available_trucks;
        for(int i = 0; i < n_trucks; ++i)
            if (!full_truck[i])
                available_trucks.push_back(i);
        //std::cout << "avail trucks = " << available_trucks.size() << "\n";
        if (available_trucks.size() == 0) break;

        int n_truck_continue = (available_trucks.size() * 2) / 2;
        std::vector<int> continue_trucks_id = choose_trucks(tours, available_trucks, n_truck_continue);
        std::vector<Tour*> continue_trucks;
        //std::cout << "choose trucks to continue done:\n"; for(int v : continue_trucks_id) std::cout << v << " "; std::cout << "\n";
        for(int v : continue_trucks_id) continue_trucks.push_back(&tours[v]);

        std::vector<int> customer_left;
        for(int i = 1; i < n; ++i)
            if (!visited_customer[i])
                customer_left.push_back(i);
        if (customer_left.size() == 0) break;
        int n_virtual = (int)(Utilities::rnd(0.1, 0.5) * available_trucks.size());
        n_virtual = std::max(1, n_virtual);
        for(int i = 0; i < n_virtual; ++i) customer_left.push_back(-1);
        //std::cout << "customers left = " << customer_left.size() << "\n"; for(int v : customer_left) std::cout << v << " "; std::cout << "\n";

        std::vector<int> customer_visited = visit_customers(continue_trucks, customer_left);
        //std::cout << "customers visited this turn:\n"; for(int v : customer_visited) std::cout << v << " "; std::cout << "\n";
        for(int v : customer_visited) visited_customer[v] = true;

        for(int i = 0; i < n_trucks; ++i) {
            /*bool ok = false;
            for(int j = 0; j < n; ++j)
                if ((!visited_customer[j]) && (tours[i].tot_demand + instance->demand[j] <= truck_capacity)) {
                    ok = true;
                }
            if (!ok) full_truck[i] = true;*/
            if ((tours[i].size() > 1) && (tours[i].back() == 0))
                full_truck[i] = true;
        }
    }

    for(int i = 0; i < n_trucks; ++i) tours[i].add(0, instance->demand[0]);
    std::cout << "single ant done\n";
    Solution res(tours);
    std::cout << "before fix:\n";
    res.print(); res.valid(*instance);
    std::cout << "tot len = " << res.get_tot_len(instance->distance) << "\n";

    fix_solution(res);
    std::cout << "after fix:\n";
    res.print(); res.valid(*instance);
    std::cout << "tot len = " << res.get_tot_len(instance->distance) << "\n";

    return res;
}

std::vector<Solution> ACO::generate_solutions(int n_ants) {
    std::vector<Solution> sols;
    for(int i = 0; i < n_ants; ++i) {
        Solution x = ant_construct_solution();
        sols.push_back(x);
    }
    return sols;
}

int ACO::best_solution(std::vector<Solution> &sols) {
    int sol_idx = 0;
    double min_cost = sols[0].get_cost(instance->distance, truck_capacity);
    for(int i = 1; i < sols.size(); ++i) {
        double cost = sols[i].get_cost(instance->distance, truck_capacity);
        if (cost < min_cost) {
            min_cost= cost;
            sol_idx = i;
        }
    }
    return sol_idx;
}

void ACO::update_trails(Solution &sol) {
    //std::cout << "update trail single sol begin...\n";
    sol.print();
    std::set< std::pair<int, int> > good_edges;
    for(int i = 0; i < sol.tours.size(); ++i)
        for(int j = 0; j + 1 < sol.tours[i].size(); ++j) {
            int u = sol.tours[i][j], v = sol.tours[i][j + 1];
            //std::cout << u << " " << v << "\n";
            good_edges.insert(std::make_pair(u, v));
            good_edges.insert(std::make_pair(v, u));
            trail[u][v] = (1 - rho) * trail[v][u] + rho * T_min;
            trail[v][u] = (1 - rho) * trail[v][u] + rho * T_min;
        }
    //std::cout << "hi\n";
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            if (good_edges.count(std::make_pair(i, j)) == 0) {
                trail[i][j] = (1 - rho) * trail[i][j] + rho * T_max;
            }

    double s = 0, mint = 1e15, maxt = -1.0;
    for(int i = 0; i < n; ++i)
        for(int j = i + 1; j < n; ++j) {
            s += trail[i][j];
            mint = std::min(mint, trail[i][j]);
            maxt = std::max(maxt, trail[i][j]);
        }
    std::cerr << "sum trail = " << s << ", min = " << mint << ", max = " << maxt << "\n";
}

void ACO::update_trails(std::vector<Solution> &sols) {
    int sol_idx = best_solution(sols);
    update_trails(sols[sol_idx]);
}

void ACO::run() {
    auto start_time = clock();

    init_ACO();
    int n_ants = 3;
    int n_rounds = 1000;
    Solution best_sol;

    while (true) {
        std::cout << n_rounds << " round(s) left\n";
        std::cerr << n_rounds << " round(s) left\n";
        if (n_rounds <= 0) break;
        --n_rounds;
        std::vector<Solution> sols = generate_solutions(n_ants);
        //std::cout << "all ants generate solutions done\n";
        int sol_idx = best_solution(sols);
        //std::cout << "find this turn best solution = " << sol_idx << "\n";
        std::cerr << "cur len = " << sols[sol_idx].get_tot_len(instance->distance) << "\n";
        std::cerr << "n tour exceed = " << sols[sol_idx].count_exceed_tour(truck_capacity) << "\n";
        if ((best_sol.tours.size() == 0) || (sols[sol_idx].get_tot_len(instance->distance) < best_sol.get_tot_len(instance->distance))) {
            best_sol = sols[sol_idx];
        }

        std::cerr << "best len = " << best_sol.get_tot_len(instance->distance) <<" at "<<(clock() - start_time)/CLOCKS_PER_SEC<<"\n";

        update_trails(sols);
        //std::cout << "update trails done\n";
    }

    std::cout << "ACO done\n";
    std::cout << "best solution:\n";
    best_sol.print();
    best_sol.valid(*instance);
    std::cout << "tot len " << best_sol.get_tot_len(instance->distance) << "\n";
    std::cout << "tot cost " << best_sol.get_cost(instance->distance, truck_capacity) << "\n";
    std::cout << "freq = \n"; Utilities::print_vector_2D(freq);

    auto end_time = clock();
    double running_time = (double)(end_time - start_time)/CLOCKS_PER_SEC;
    std::cout << "running time = " << running_time << "\n";
}

#endif // ACO_HPP_INCLUDED
