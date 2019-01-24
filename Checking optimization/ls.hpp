#ifndef LS_HPP_INCLUDED
#define LS_HPP_INCLUDED

/*

       AAAA    CCCC   OOOO   TTTTTT   SSSSS  PPPPP
      AA  AA  CC     OO  OO    TT    SS      PP  PP
      AAAAAA  CC     OO  OO    TT     SSSS   PPPPP
      AA  AA  CC     OO  OO    TT        SS  PP
      AA  AA   CCCC   OOOO     TT    SSSSS   PP

######################################################
##########    ACO algorithms for the TSP    ##########
######################################################

      Version: 1.0
      File:    ls.h
      Author:  Thomas Stuetzle
      Purpose: header file for local search routines
      Check:   README and gpl.txt
      Copyright (C) 1999  Thomas Stuetzle
*/

/***************************************************************************

    Program's name: acotsp

    Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS, BWAS) for the
    symmetric TSP

    Copyright (C) 2004  Thomas Stuetzle

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    email: stuetzle no@spam ulb.ac.be
    mail address: Universite libre de Bruxelles
                  IRIDIA, CP 194/6
                  Av. F. Roosevelt 50
                  B-1050 Brussels
		  Belgium

***************************************************************************/

#include "TSP.hpp"
#include "utilities.hpp"
#include <vector>

long int ls_flag;

long int nn_ls;

long int dlb_flag;

//long int * generate_random_permutation( long int n, long *seed);

long int * generate_random_permutation( long int n, long *seed)
/*
      FUNCTION:       generate a random permutation of the integers 0 .. n-1
      INPUT:          length of the array
      OUTPUT:         pointer to the random permutation
      (SIDE)EFFECTS:  the array holding the random permutation is allocated in this
                      function. Don't forget to free again the memory!
      COMMENTS:       only needed by the local search procedures
*/
{
   long int  i, help, node, tot_assigned = 0;
   double    rnd;
   long int  *r;

   r = (long int*)malloc(n * sizeof(long int));

   for ( i = 0 ; i < n; i++)
     r[i] = i;

   for ( i = 0 ; i < n ; i++ ) {
     /* find (randomly) an index for a free unit */
     rnd  = ran01 ( seed );
     node = (long int) (rnd  * (n - tot_assigned));
     assert( i + node < n );
     help = r[i];
     r[i] = r[i+node];
     r[i+node] = help;
     tot_assigned++;
   }
   return r;
}

double ran01( long *idum )
/*
      FUNCTION:       generate a random number that is uniformly distributed in [0,1]
      INPUT:          pointer to variable with the current seed
      OUTPUT:         random number uniformly distributed in [0,1]
      (SIDE)EFFECTS:  random number seed is modified (important, this has to be done!)
      ORIGIN:         numerical recipes in C
*/
{
  long k;
  double ans;

  k =(*idum)/IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0 ) *idum += IM;
  ans = AM * (*idum);
  return ans;
}

//template<typename T>
//void two_opt_first( long int *tour, int n, problem &instance, std::vector< std::vector<T> > &distance, long *seed);

template<typename T>
void two_opt_first( long int *tour, int n, problem &instance, std::vector< std::vector<T> > &distance, long *seed)
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
    long int c1, c2;             /* cities considered for an exchange */
    long int s_c1, s_c2;         /* successor cities of c1 and c2     */
    long int p_c1, p_c2;         /* predecessor cities of c1 and c2   */
    long int pos_c1, pos_c2;     /* positions of cities c1, c2        */
    long int i, j, h, l;
    long int improvement_flag, help, n_improves = 0, n_exchanges=0;
    long int h1=0, h2=0, h3=0, h4=0;
    T radius;             /* radius of nn-search */
    T gain = 0;
    long int *random_vector;
    long int *pos;               /* positions of cities in tour */
    long int *dlb;               /* vector containing don't look bits */

    double eps = 1e-9;

    pos = (long int*)malloc(n * sizeof(long int));
    dlb = (long int*)malloc(n * sizeof(long int));
    for ( i = 0 ; i < n ; i++ ) {
        pos[tour[i]] = i;
        dlb[i] = FALSE;
    }

    improvement_flag = TRUE;
    random_vector = generate_random_permutation( n, seed);

    while ( improvement_flag ) {
        //std::cout << "in while...\n";
        improvement_flag = FALSE;

        for (l = 0 ; l < n; l++) {
            //std::cout << l << "/" << n << "\n";
            c1 = random_vector[l];
            //std::cout << "c1 = " << c1 << "\n";
            DEBUG ( assert ( c1 < n && c1 >= 0); )
            if ( dlb_flag && dlb[c1] )
                continue;
            pos_c1 = pos[c1];
            s_c1 = tour[pos_c1+1];
            radius = distance[c1][s_c1];

            //std::cout << "mark 1\n";
            //std::cout << "s_c1 = " << s_c1 << "\n";
            /* First search for c1's nearest neighbours, use successor of c1 */
            for ( h = 0 ; h < nn_ls ; h++ ) {
                //std::cout << "h = " << h << "/" << nn_ls << "\n";
                c2 = instance.nn_list[c1][h]; /* exchange partner, determine its position */
                //std::cout << "c1 = " << c1 << ", c2 = " << c2 << "\n";
                if (( radius > distance[c1][c2] ) && (abs(radius - distance[c1][c2]) > eps)) {
                    //std::cout << "pos[c2]+1 = " << pos[c2] << " " << pos[c2] + 1 << "\n";
                    s_c2 = tour[pos[c2]+1];
                    gain =  - radius + distance[c1][c2] +
                    distance[s_c1][s_c2] - distance[c2][s_c2];
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
            radius = distance[p_c1][c1];
            //std::cout << "mark 2.2\n";
            for ( h = 0 ; h < nn_ls ; h++ ) {
                //std::cout << "h = " << h << "/" << nn_ls << "\n";
                c2 = instance.nn_list[c1][h];  /* exchange partner, determine its position */
                if (( radius > distance[c1][c2] ) && (abs(radius - distance[c1][c2]) > eps)) {
                    pos_c2 = pos[c2];
                    if (pos_c2 > 0)
                    p_c2 = tour[pos_c2-1];
                    else
                    p_c2 = tour[n-1];
                    if ( p_c2 == c1 )
                    continue;
                    if ( p_c1 == c2 )
                    continue;
                    gain =  - radius + distance[c1][c2] +
                    distance[p_c1][p_c2] - distance[p_c2][c2];
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
                dlb[c1] = TRUE;
                continue;

            exchange2opt:
            //std::cout << "go to exchange2opt\n";
            n_exchanges++;
            improvement_flag = TRUE;
            dlb[h1] = FALSE; dlb[h2] = FALSE;
            dlb[h3] = FALSE; dlb[h4] = FALSE;
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

void two_h_opt_first( long int *tour, int n, problem &instance, long *seed);

void three_opt_first( long int *tour, int n, problem &instance, long *seed);


#endif // LS_HPP_INCLUDED
