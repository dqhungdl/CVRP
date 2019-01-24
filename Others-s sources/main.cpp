#include <bits/stdc++.h>

#include "structs.hpp"
#include "ACO.hpp"
#include "Utilities.hpp"

using namespace std;

Problem instance;
int n_trucks, truck_capacity;
int demand[100000];

void read_input(string inp_file) {
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
    instance.n = n;
    instance.distance = vector< vector<double> > (n, vector<double>(n, 0));
    instance.nn_list = vector< vector<int> >(n, vector<int>(n - 1, 0));
    instance.points = vector<Point>(n);

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
    instance.demand.resize(n);
    for(int i = 0; i < n; ++i) instance.demand[i] = demand[i];

    f >> tmp;
    int depot; f >> depot;
    assert(depot == 1);
    int v; f >> v;
    assert(v == -1);

    f.close();

    //compute distance between points
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j) {
            instance.distance[i][j] = instance.points[i].round_Euclid_distance(instance.points[j]);
        }

    //compute nn_list
    for(int i = 0; i < n; ++i) {
        vector< pair<double, int> > nei;
        for(int j = 0; j < n; ++j)
            if (j != i) {
                nei.push_back(make_pair(instance.distance[i][j], j));
            }
        sort(nei.begin(), nei.end());
        for(int j = 0; j < n - 1; ++j)
            instance.nn_list[i][j] = nei[j].second;
    }
}

int main() {
    freopen("log.txt", "w", stdout);

    string input_file = "X-n101-k25.vrp";
    read_input(input_file);
    cout << "read done\n";
    instance.n_trucks = 26;
    instance.truck_capacity = truck_capacity;

    ACO solver(&instance);
    solver.run();

    return 0;
}
