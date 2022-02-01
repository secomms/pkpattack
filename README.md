# pkpattack
Source code for simulating and estimating the complexity of new approaches for solving the permuted kernel problem

The script test_attack.sage simulates the full attack using low-weight codewords; the performances of the algorithm (such as lists sizes and number of collisions) are compared with the theoretical estimates. To launch the script, type

load('test_attack.sage')

The script cost.sage computed the minimum cost of our attack. To launch the script, type

load('cost.sage')
