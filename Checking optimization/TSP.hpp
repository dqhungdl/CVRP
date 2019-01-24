#ifndef TSP_HPP_INCLUDED
#define TSP_HPP_INCLUDED

#include <vector>

struct point {
  double x;
  double y;
};

struct problem {
  long int      optimum;                /* optimal tour length if known, otherwise a bound */
  long int      n;                      /* number of cities */
  long int      n_near;                 /* number of nearest neighbors */
  std::vector<point> points;               /* array of structs containing coordinates of nodes */
  std::vector< std::vector<long int> > distance;	        /* distance matrix: distance[i][j] gives distance
					   between city i und j */
  std::vector< std::vector<double> > distance2;
  std::vector< std::vector< long int> > nn_list;              /* nearest neighbor list; contains for each node i a
                                           sorted list of n_near nearest neighbors */
};



#endif // TSP_HPP_INCLUDED
