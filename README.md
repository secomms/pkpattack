# pkpattack
Source code for simulating and estimating the complexity of algorithms for solving the permuted kernel problem.

The script test_attack.sage simulates our attack; the performances of the algorithm (such as lists sizes and number of collisions) are compared with the theoretical estimates. To launch the script, type

load('test_attack.sage')

The pkp parameters (such as n, q and ell) and attack paramaters need to be specified in the script. The files list_sorting.sage and ISD.sage contains useful functions to simulate the attack.

The script kmp_cost.sage computes the minimum cost of KMP attack. 

The script our_cost.sage computes the minimum cost of our attack.

The script compare_cost.sage can be used to compare KMP and our algorithms, for fixed values of n, R and increasing values of ell; the value of q is computed accordingly.
