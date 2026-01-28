The code computes archetypes of  a data matrix X and provides r archetypes.
The default value is r=4. 

The code requires essntially the data matrix X. X could be 2-dimensional,
say  n*m,  where n is the sample size, and m is the number of variables. 
It could also be 3-dimensional for example, n*la*lo, for a field, in which
n again is the sample size, and la is the latitudinal dimension and lo is 
the longitudinal dimension

Example: [A B Za Zb]=my_aaals_4github(randn(100,2));

A and B are weight matrices Zb are the archetypes, Za are similar to archetypes
 but have smaller amplitudes (referred to in the ms as duals).
 
In some cases the matrix may be scaled to get smoother convergence.
