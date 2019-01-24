#ifndef UTILITIES_HPP_INCLUDED
#define UTILITIES_HPP_INCLUDED

#include <bits/stdc++.h>

namespace Utilities {

    //return a random number in range [0.0, 1.0]
    double rnd() {
        return (double)(std::rand()) / (RAND_MAX);
    }

    //return a random (int) number in range [u, v]
    int rnd(int u, int v) {
        return u + std::rand() % (v - u + 1);
    }

    //return a random (float) number in range [u, v]
    double rnd(double u, double v) {
        return u + (v - u) * rnd();
    }

    template<typename T>
    int choice(std::vector<T> &a) {
        T sum = std::accumulate(a.begin(), a.end(), (T)(0));
        double p = rnd();
        double sx = 0.0;
        for(int i = 0; i < int(a.size()); ++i) {
            double x = (double)(a[i]) / sum;
            sx += x;
            if (sx >= p) return i;
        }
        return (int)(a.size()) - 1;
    }

    bool equals(double a, double b, double eps) {
        return (fabs(a - b) <= eps);
    }

    const double eps = 1e-9;
    bool equals(double a, double b) {
        return equals(a, b, eps);
    }

    template<typename T>
    void print_vector_1D(std::vector<T> &a) {
        for(int v : a) std::cout << v << " "; std::cout << "\n";
    }

    template<typename T>
    void print_vector_2D(std::vector< std::vector<T> > &a) {
        for(std::vector<T> &b : a) {
            print_vector_1D(b);
        }
    }
}

#endif // UTILITIES_HPP_INCLUDED
