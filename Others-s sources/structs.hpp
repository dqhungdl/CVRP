#ifndef STRUCTS_HPP_INCLUDED
#define STRUCTS_HPP_INCLUDED

#include <bits/stdc++.h>

struct Point {
    double x, y;

    Point() {}
    Point(int _x, int _y) {
        x = _x; y = _y;
    }

    int round_Euclid_distance(Point &other) {
        double v = sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y));
        return (int)(v + 0.5);
    }
};

struct Problem {
    int n; //number of depot + customers, labeled 0->n-1. Depot is labeled 0
    std::vector<Point> points;
    std::vector< std::vector<double> > distance;
    std::vector< std::vector<int> > nn_list; //nearest neighbor list
    int n_trucks, truck_capacity;
    std::vector<int> demand;

    void init(int _n) {
        n = _n;
        points.resize(n);
        distance.resize(n);
        for(int i = 0; i < n; ++i) distance.resize(n);
        nn_list.resize(n);
        for(int i = 0; i < n; ++i) nn_list.resize(n - 1);
        demand.resize(n);
    }

    Problem() {}
    Problem(int _n) {
        init(_n);
    }
};

struct Tour {
    int id;
    std::vector<int> order;
    int tot_demand = 0;

    Tour() {}

    int& operator[] (int idx) {
        return order[idx];
    }

    int size() {
        return order.size();
    }

    int back() {
        return order.back();
    }

    void add(int v, int demand) {
        order.push_back(v);
        tot_demand += demand;
    }

    int calc_tot_demand(std::vector<int> &demand) {
        tot_demand = 0;
        for(int v : order) tot_demand += demand[v];
        return tot_demand;
    }

    double get_len(std::vector< std::vector<double> > &distance) {
        double res = 0.0;
        for(int i = 0; i + 1 < order.size(); ++i) res += distance[ order[i] ][ order[i + 1] ];
        return res;
    }

    double get_cost(std::vector< std::vector<double> > &distance, int truck_capacity) {
        double res = get_len(distance);
        double penalty = 10000;
        res += penalty * std::max(0, tot_demand - truck_capacity);
        return res;
    }

    void print() {
        for(int v : order) std::cout << v << " "; std::cout << "\n";
    }
};

struct Solution {
    std::vector<Tour> tours;

    Solution() {}
    Solution(std::vector<Tour> &_tours) {
        tours = _tours;
    }

    double get_tot_len(std::vector< std::vector<double> > &distance) {
        double res = 0.0;
        for(int i = 0; i < tours.size(); ++i)
            res += tours[i].get_len(distance);
        return res;
    }

    double get_cost(std::vector< std::vector<double> > &distance, int truck_capacity) {
        double res = 0.0;
        for(int i = 0; i < tours.size(); ++i)
            res += tours[i].get_cost(distance, truck_capacity);
        return res;
    }

    void print() {
        for(Tour &t : tours) t.print();
    }

    int count_exceed_tour(int truck_capacity) {
        int res = 0;
        for(Tour& t : tours)
            if (t.tot_demand > truck_capacity)
                ++res;
        return res;
    }

    bool valid(Problem &instance) {
        std::vector<bool> used(instance.n, false);
        for(Tour &t : tours) {
            assert(t.size() >= 2);
            assert((t[0] == 0) && (t.back() == 0));
            for(int i = 1; i + 1 < t.size(); ++i) {
                assert(used[t[i]] == false);
                used[t[i]] = true;
            }
        }
        assert(used[0] == false);
        for(int i = 1; i < instance.n; ++i) assert(used[i] == true);

        //check demand
        int n_exceed = 0;
        for(Tour &t : tours) {
            if (t.tot_demand > instance.truck_capacity) {
                n_exceed++;
                std::cout << "------------\n";
                t.print();
                std::cout << "Exceed allowed capacity: " << t.tot_demand << "/" << instance.truck_capacity << "\n";
                std::cout << "------------\n";
            }
        }
        std::cout << "n tours exceed capacity = " << n_exceed << "/" << tours.size() << "\n";
        return true;
    }
};

#endif // STRUCTS_HPP_INCLUDED
