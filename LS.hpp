#ifndef LS_HPP_INCLUDED
#define LS_HPP_INCLUDED

#include <bits/stdc++.h>

#include "structs.hpp"
#include "Utilities.hpp"

namespace LS {
    long int dlb_flag = 1;  /* flag indicating whether don't look bits are used. I recommend
			      to always use it if local search is applied */

    long int* generate_random_permutation(int n) {
        std::vector<int> a(n);
        for(int i = 0; i < n; ++i) a[i] = i;

        long int* res = (long int*)malloc(n * sizeof(long int));
        for(int i = 0; i < n; ++i) {
            int j = Utilities::rnd(0, n - 1 - i);
            res[i] = a[j];
            std::swap(a[j], a[n - 1 - i]);
        }
        return res;
    }

    void two_opt_first( long int *tour, int n, Problem &instance )
    /*
          FUNCTION:       2-opt a tour
          INPUT:          pointer to the tour that undergoes local optimization
          OUTPUT:         none
          (SIDE)EFFECTS:  tour is 2-opt
          COMMENTS:       the neighbourhood is scanned in random order (this need
                          not be the best possible choice). Concerning the speed-ups used
                          here consult, for example, Chapter 8 of
                          Holger H. Hoos and Thomas Stuetzle,
                          Stochastic Local Search---Foundations and Applications,
                          Morgan Kaufmann Publishers, 2004.
                          or some of the papers online available from David S. Johnson.
    */
    {
        int nn_ls = n - 1; /* maximal depth of nearest neighbour lists used in the local search */

        long int c1, c2;             /* cities considered for an exchange */
        long int s_c1, s_c2;         /* successor cities of c1 and c2     */
        long int p_c1, p_c2;         /* predecessor cities of c1 and c2   */
        long int pos_c1, pos_c2;     /* positions of cities c1, c2        */
        long int i, j, h, l;
        long int improvement_flag, help, n_improves = 0, n_exchanges=0;
        long int h1=0, h2=0, h3=0, h4=0;
        double radius;             /* radius of nn-search */
        double gain = 0;
        long int *random_vector;
        long int *pos;               /* positions of cities in tour */
        long int *dlb;               /* vector containing don't look bits */

        double eps = 1e-9;

        pos = (long int*)malloc(n * sizeof(long int));
        dlb = (long int*)malloc(n * sizeof(long int));
        for ( i = 0 ; i < n ; i++ ) {
            pos[tour[i]] = i;
            dlb[i] = false;
        }

        improvement_flag = true;
        random_vector = generate_random_permutation( n );

        while ( improvement_flag ) {
            //std::cout << "in while...\n";
            improvement_flag = false;

            for (l = 0 ; l < n; l++) {
                //std::cout << l << "/" << n << "\n";
                c1 = random_vector[l];
                //std::cout << "c1 = " << c1 << "\n";
                assert ( c1 < n && c1 >= 0);
                if ( dlb_flag && dlb[c1] )
                    continue;
                pos_c1 = pos[c1];
                s_c1 = tour[pos_c1+1];
                radius = instance.distance[c1][s_c1];

                //std::cout << "mark 1\n";
                //std::cout << "s_c1 = " << s_c1 << "\n";
                /* First search for c1's nearest neighbours, use successor of c1 */
                for ( h = 0 ; h < nn_ls ; h++ ) {
                    //std::cout << "h = " << h << "/" << nn_ls << "\n";
                    c2 = instance.nn_list[c1][h]; /* exchange partner, determine its position */
                    //std::cout << "c1 = " << c1 << ", c2 = " << c2 << "\n";
                    if (( radius > instance.distance[c1][c2] ) && (abs(radius - instance.distance[c1][c2]) > eps)) {
                        //std::cout << "pos[c2]+1 = " << pos[c2] << " " << pos[c2] + 1 << "\n";
                        s_c2 = tour[pos[c2]+1];
                        gain =  - radius + instance.distance[c1][c2] +
                                instance.distance[s_c1][s_c2] - instance.distance[c2][s_c2];
                        if (( gain < 0 ) && (abs(gain) > eps)) {
                            h1 = c1; h2 = s_c1; h3 = c2; h4 = s_c2;
                            goto exchange2opt;
                        }
                    }
                    else
                        break;
                }

                //std::cout << "mark 2\n";
                /* Search one for next c1's h-nearest neighbours, use predecessor c1 */
                //std::cout << "pos_c1 = " << pos_c1 << "\n";
                if (pos_c1 > 0)
                    p_c1 = tour[pos_c1-1];
                else
                    p_c1 = tour[n-1];
                //std::cout << "mark 2.1\n";
                //std::cout << p_c1 << " " << c1 << "\n";
                radius = instance.distance[p_c1][c1];
                //std::cout << "mark 2.2\n";
                for ( h = 0 ; h < nn_ls ; h++ ) {
                    //std::cout << "h = " << h << "/" << nn_ls << "\n";
                    c2 = instance.nn_list[c1][h];  /* exchange partner, determine its position */
                    if (( radius > instance.distance[c1][c2] ) && (abs(radius - instance.distance[c1][c2]) > eps)) {
                        pos_c2 = pos[c2];
                        if (pos_c2 > 0)
                            p_c2 = tour[pos_c2-1];
                        else
                            p_c2 = tour[n-1];
                        if ( p_c2 == c1 )
                            continue;
                        if ( p_c1 == c2 )
                            continue;
                        gain =  - radius + instance.distance[c1][c2] +
                            instance.distance[p_c1][p_c2] - instance.distance[p_c2][c2];
                        if (( gain < 0 ) && (abs(gain) > eps)) {
                            h1 = p_c1; h2 = c1; h3 = p_c2; h4 = c2;
                            goto exchange2opt;
                        }
                    }
                    else
                        break;
                }
                //std::cout << "mark 2.5\n";
                    /* No exchange */
                    dlb[c1] = 1;
                    continue;

                exchange2opt:
                //std::cout << "go to exchange2opt\n";
                n_exchanges++;
                improvement_flag = 1;
                dlb[h1] = 0; dlb[h2] = 0;
                dlb[h3] = 0; dlb[h4] = 0;
                //std::cout << "mark 3\n";
                /* Now perform move */
                if ( pos[h3] < pos[h1] ) {
                    help = h1; h1 = h3; h3 = help;
                    help = h2; h2 = h4; h4 = help;
                }
                if ( pos[h3] - pos[h2] < n / 2 + 1) {
                    /* reverse inner part from pos[h2] to pos[h3] */
                    i = pos[h2]; j = pos[h3];
                    while (i < j) {
                        c1 = tour[i];
                        c2 = tour[j];
                        tour[i] = c2;
                        tour[j] = c1;
                        pos[c1] = j;
                        pos[c2] = i;
                        i++; j--;
                    }
                }
                else {
                    /* reverse outer part from pos[h4] to pos[h1] */
                    i = pos[h1]; j = pos[h4];
                    if ( j > i )
                        help = n - (j - i) + 1;
                    else
                        help = (i - j) + 1;
                    help = help / 2;
                    for ( h = 0 ; h < help ; h++ ) {
                        c1 = tour[i];
                        c2 = tour[j];
                        tour[i] = c2;
                        tour[j] = c1;
                        pos[c1] = j;
                        pos[c2] = i;
                        i--; j++;
                        if ( i < 0 )
                            i = n-1;
                        if ( j >= n )
                            j = 0;
                    }
                    tour[n] = tour[0];
                }
            }

            if ( improvement_flag ) {
                n_improves++;
            }
        }
        free( random_vector );
        free( dlb );
        free( pos );
    }

    double tour_cost(std::vector<int> &tour, int n, std::vector< std::vector<double> > &distance) {
        double res = 0.0;
        for(int i = 0; i < n; ++i) 
            res += distance[tour[i]][tour[i + 1]];
        return res;
    }

    Problem extract_problem(long int* tour, int n, Problem &instance) {
        Problem res;
        res.n = n;
        res.distance = std::vector< std::vector<double> > (n, std::vector<double>(n));
        for(int i = 0; i < n; ++i)
            for(int j = 0; j < n; ++j) {
                res.distance[i][j] = instance.distance[tour[i]][tour[j]];
            }
        res.nn_list = std::vector< std::vector<int> > (n, std::vector<int>(n - 1));
        for(int i = 0; i < n; ++i) {
            std::vector< std::pair<double, int> > nei;
            for(int j = 0; j < n; ++j)
                if (j != i) nei.push_back(std::make_pair(res.distance[i][j], j));
            sort(nei.begin(), nei.end());
            for(int j = 0; j < n - 1; ++j)
                res.nn_list[i][j] = nei[j].second;
        }
        return res;
    }

    //tour includes depot (0)
    void ls(long int* tour, int n, Problem &instance) {
        //cout << "ls begin\n";
        tour[n] = tour[0];
        //return tour_cost(tour, n);

        //for(int i = 0; i < n; ++i) cout << tour[i] << " "; cout << "\n";
        long int* new_tour = (long int*)malloc((n + 1) * sizeof(long int));
        for(int i = 0; i < n; ++i) new_tour[i] = i; new_tour[n] = new_tour[0];
        //for(int i = 0; i <= n; ++i) cout << new_tour[i] << " "; cout << "\n";

        Problem new_instance = extract_problem(tour, n, instance);

        //cout << "call ls()\n";
        two_opt_first(new_tour, n, new_instance);
        //cout << "call done\n";
        assert(new_tour[0] == new_tour[n]);

        //convert
        std::vector<int> init_tour(n);
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
    }

    double ls_vector(std::vector<int> &a, Problem &instance) {
        assert(a[0] == 0); assert(a.back() == 0);
        double init_cost = tour_cost(a, a.size() - 1, instance.distance);
        long int* tour = (long int*)malloc(a.size() * sizeof(long int));
        for(int i = 0; i < int(a.size()); ++i) tour[i] = a[i];
        ls(tour, a.size() - 1, instance);
        for(int i = 0; i < int(a.size()); ++i) a[i] = tour[i];
        free(tour);

        double new_cost = tour_cost(a, a.size() - 1, instance.distance);
        assert((new_cost <= init_cost) || (fabs(new_cost - init_cost) <= 1e-9));
        return new_cost;
    }

    double ls_tour(Tour &a, Problem &instance) {
        assert(a[0] == 0); assert(a.back() == 0);
        double init_cost = a.get_len(instance.distance);
        long int* tour = (long int*)malloc(a.size() * sizeof(long int));
        for(int i = 0; i < a.size(); ++i) tour[i] = a[i];

        ls(tour, a.size() - 1, instance);

        for(int i = 0; i < a.size(); ++i) a[i] = tour[i];
        free(tour);
        double new_cost = a.get_len(instance.distance);
        assert((new_cost <= init_cost) || (fabs(new_cost - init_cost) <= 1e-9));
        return new_cost;
    }

    bool exchange_2tour(std::vector<int> &a, std::vector<int> &b, int len_cut, Problem &instance) {
        int truck_capacity = instance.truck_capacity;
        std::vector<int>& demand = instance.demand;
        std::vector< std::vector<double> >& distance = instance.distance;

        double INF = 1e15;

        int m = a.size(), n = b.size();
        assert(a[0] == 0); assert(a[m - 1] == 0);
        assert(b[0] == 0); assert(b[n - 1] == 0);

        int tot_cap_a = 0;
        for(int v : a) tot_cap_a += demand[v];
        int tot_cap_b = 0;
        for(int v : b) tot_cap_b += demand[v];
        int tot_cap = tot_cap_a + tot_cap_b;
        //std::cout << "tot cap = " << tot_cap << "\n";
        int max_capacity = std::max(truck_capacity, std::max(tot_cap_a, tot_cap_b));

        std::vector< std::vector< std::vector<double> > > f(m + 1, std::vector< std::vector<double> >(n + 1, std::vector<double>(max_capacity + 1, INF)));
        std::vector< std::vector< std::vector< std::pair<int, int> > > > trace(m + 1, std::vector< std::vector< std::pair<int, int> > >(n + 1, std::vector< std::pair<int, int> >(max_capacity + 1)));
        std::vector< std::vector< std::vector<int> > > trace_c1(m + 1, std::vector< std::vector<int> >(n + 1, std::vector<int>(max_capacity + 1, -1)));

        f[0][0][0] = 0.0;
        int sda = 0;
        for(int i = 0; i < m; ++i) {
            sda += demand[a[i]];
            int sdb = 0;
            for(int j = 0; j < n; ++j) {
                sdb += demand[b[j]];
                for(int c1 = 0; c1 <= max_capacity; ++c1) {
                    int c2 = sda + sdb - c1;
                    if (!(f[i][j][c1] + 10.0 < INF)) continue;
                    assert(c2 <= max_capacity);
                    //cout << i << " " << j << " " << c1 << " " << f[i][j][c1] << "\n";

                    //i -> i + 1
                    if ((i + 1 < m) && (c1 + demand[a[i + 1]] <= max_capacity)) {
                        if (f[i + 1][j][c1 + demand[a[i + 1]]] > f[i][j][c1] + distance[a[i]][a[i + 1]]) {
                            f[i + 1][j][c1 + demand[a[i + 1]]] = f[i][j][c1] + distance[a[i]][a[i + 1]];
                            trace[i + 1][j][c1 + demand[a[i + 1]]] = std::make_pair(i, j);
                            trace_c1[i + 1][j][c1 + demand[a[i + 1]]] = c1;
                        }

                        double dist_cut = 0.0;
                        int next_c1 = c1 + demand[a[i + 1]];
                        for(int l = 1; l <= len_cut; ++l) {
                            if (j + l + 1 >= n) break;
                            next_c1 += demand[b[j + l]];
                            if (next_c1 > max_capacity) break;
                            if (l > 1) dist_cut += distance[b[j + l - 1]][b[j + l]];

                            double add_cost = dist_cut + distance[a[i]][b[j + 1]] + distance[b[j + l]][a[i + 1]] + distance[b[j]][b[j + l + 1]];
                            if (c2 + demand[b[j + l + 1]] <= max_capacity) {
                                if (f[i + 1][j + l + 1][next_c1] > f[i][j][c1] + add_cost) {
                                    f[i + 1][j + l + 1][next_c1] = f[i][j][c1] + add_cost;
                                    trace[i + 1][j + l + 1][next_c1] = std::make_pair(i, j);
                                    trace_c1[i + 1][j + l + 1][next_c1] = c1;
                                }
                            }
                        }
                    }

                    //j -> j + 1
                    if ((j + 1 < n) && (c2 + demand[b[j + 1]] <= max_capacity)) {
                        if (f[i][j + 1][c1] > f[i][j][c1] + distance[b[j]][b[j + 1]]) {
                            f[i][j + 1][c1] = f[i][j][c1] + distance[b[j]][b[j + 1]];
                            trace[i][j + 1][c1] = std::make_pair(i, j);
                            trace_c1[i][j + 1][c1] = c1;
                        }

                        double dist_cut = 0.0;
                        int next_c2 = c2 + demand[b[j + 1]];
                        for(int l = 1; l <= len_cut; ++l) {
                            if (i + l + 1 >= m) break;
                            next_c2 += demand[a[i + l]];
                            if (next_c2 > max_capacity) break;
                            if (l > 1) dist_cut += distance[a[i + l - 1]][a[i + l]];

                            double add_cost = dist_cut + distance[b[j]][a[i + 1]] + distance[a[i + l]][b[j + 1]] + distance[a[i]][a[i + l + 1]];
                            if (c1 + demand[a[i + l + 1]] <= max_capacity) {
                                if (f[i + l + 1][j + 1][c1 + demand[a[i + l + 1]]] > f[i][j][c1] + add_cost) {
                                    f[i + l + 1][j + 1][c1 + demand[a[i + l + 1]]] = f[i][j][c1] + add_cost;
                                    trace[i + l + 1][j + 1][c1 + demand[a[i + l + 1]]] = std::make_pair(i, j);
                                    trace_c1[i + l + 1][j + 1][c1 + demand[a[i + l + 1]]] = c1;
                                }
                            }
                        }
                    }

                    //i -> j + 1, j->i + 1
                    if ((i + 1 < m) && (j + 1 < n))
                        if ((c1 + demand[b[j + 1]] <= max_capacity) && (c2 + demand[a[i + 1]] <= max_capacity)) {
                            if (f[i + 1][j + 1][c2 + demand[a[i + 1]]] > f[i][j][c1] + distance[a[i]][b[j + 1]] + distance[b[j]][a[i + 1]]) {
                                f[i + 1][j + 1][c2 + demand[a[i + 1]]] = f[i][j][c1] + distance[a[i]][b[j + 1]] + distance[b[j]][a[i + 1]];
                                trace[i + 1][j + 1][c2 + demand[a[i + 1]]] = std::make_pair(i, j);
                                trace_c1[i + 1][j + 1][c2 + demand[a[i + 1]]] = c1;
                            }
                        }
                }
            }
        }

        double res = INF, penalty = 10000.0;
        int best_c1 = -1;
        for(int c1 = 0; c1 <= max_capacity; c1++) {
            int c2 = tot_cap - c1;
            double cost = f[m - 1][n - 1][c1] + penalty * ( std::max(0, c1 - truck_capacity) + std::max(0, c2 - truck_capacity) );
            if (cost < res) {
                res = cost;
                best_c1 = c1;
            }
        }
        //std::cout << "--------------------------------------------------------\n";
        double init_cost = tour_cost(a, m - 1, distance) + tour_cost(b, n - 1, distance);
        init_cost += penalty * ( std::max(0, tot_cap_a - truck_capacity) + std::max(0, tot_cap_b - truck_capacity) );
        //std::cout << "init cost = " << init_cost << "\n";
        //std::cout << "res dp = " << res << "\n";

        //std::cout << m << ": "; for(int v : a) std::cout << v << " "; std::cout << "\n";
        //std::cout << n << ": "; for(int v : b) std::cout << v << " "; std::cout << "\n";

        double eps = 1e-9;
        if ((fabs(res - init_cost) > eps) && (res > init_cost)) {
            std::cout << "dp absolutely wrongggggg!!\n";
            assert(false);
        }

        //trace and check
        int i = m - 1, j = n - 1, c1 = best_c1;
        std::vector< std::vector<int> > res_tour(2, std::vector<int>(1, 0));
        int ti = 0, tj = 1;
        while ((i > 0) || (j > 0)) {
            //std::cout << i << " " << j << "\n";
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
                    std::swap(ti, tj);
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

        //for(int v : res_tour[0]) std::cout << v << " "; std::cout << "\n";
        //for(int v : res_tour[1]) std::cout << v << " "; std::cout << "\n";

        assert(res_tour[0][0] == 0); assert(res_tour[0].back() == 0);
        assert(res_tour[1][0] == 0); assert(res_tour[1].back() == 0);

        std::set<int> s;
        for(int k = 1; k + 1 < int(a.size()); ++k) {
            assert(s.count(a[k]) == 0);
            s.insert(a[k]);
        }
        for(int k = 1; k + 1 < int(b.size()); ++k) {
            assert(s.count(b[k]) == 0);
            s.insert(b[k]);
        }
        for(int k = 1; k + 1 < int(res_tour[0].size()); ++k) {
            assert(s.count(res_tour[0][k]) == 1);
            s.erase(res_tour[0][k]);
        }
        for(int k = 1; k + 1 < int(res_tour[1].size()); ++k) {
            assert(s.count(res_tour[1][k]) == 1);
            s.erase(res_tour[1][k]);
        }
        assert(s.size() == 0);

        int cap_0 = 0;
        for(int v : res_tour[0]) cap_0 += demand[v];
        assert(cap_0 <= max_capacity);
        int cap_1 = 0;
        for(int v : res_tour[1]) cap_1 += demand[v];
        assert(cap_1 <= max_capacity);
        assert((cap_0 == best_c1) || (cap_1 == best_c1));
        assert(cap_0 + cap_1 == tot_cap);

        double check_cost = tour_cost(res_tour[0], res_tour[0].size() - 1, distance) + tour_cost(res_tour[1], res_tour[1].size() - 1, distance);
        check_cost += penalty * ( std::max(0, cap_0 - truck_capacity) + std::max(0, cap_1 - truck_capacity) );
        //std::cout << "check dp: " << check_cost << " " << res << "\n";
        assert(fabs(check_cost - res) <= eps);

        /*double ls_cost = ls_vector(res_tour[0], instance) + ls_vector(res_tour[1], instance);
        if ((ls_cost < res) && (fabs(ls_cost - res) > eps)) {
            //std::cout << "ls success!, reduced: " << res - ls_cost << "\n";
            res = ls_cost;
            //std::cout << "new res = " << res << "\n";
            //std::cout << "new tours after ls:\n";
            //for(int v : res_tour[0]) std::cout << v << " "; std::cout << "\n";
            //for(int v : res_tour[1]) std::cout << v << " "; std::cout << "\n";
        }*/

        if ((res < init_cost) && (fabs(res - init_cost) > eps)) {
            //std::cout << "found\n";
            a = res_tour[0]; b = res_tour[1];
            //assert(false);
            return true;
        }
        return false;
    }

}

#endif // LS_HPP_INCLUDED
