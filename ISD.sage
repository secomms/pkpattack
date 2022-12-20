import numpy as np

def gauss_binomial(m,r,q):

    x = 1.;
    for i in range(r):
        x = x*(1.-q^(m-i*1.))/(1.-q^(1.+i));

    return x;
##########################################################

##Lee and Brickell ISD algorithm to find subcodes
##Input:
# - Fq: finite field;
# - G: generator matrix
# - n: code length;
# - k: code dimension;
# - w: searched support size
# - d: dimension of subcodes

#The algorithm outputs three values:
# - first value: 0 if RREF is not possible, 1 if it has been computed
# - second value: 0 if no subcode has been found, 1 if ISD has found a subcode with the desired support size
# - third value: found subcode, or 0 (in case ISD has failed)

def subcodes_lee_brickell_isd(Fq,G,n,k,w,d):

    Sn = Permutations(range(n));
   
    ##Pick a random permutation
    P = Sn.random_element();
    Gp = G[:,P];
    Gp_0 = Gp[:,0:k];
   
    #Verify rank of left submatrix
    if rank(Gp_0) < k:
        del Gp;
        del Gp_0;
        return 0,0,0;
    else:

        #Go on and put matrix in systematic form
        Gp_0_inv = Gp_0^-1;
        Gp_sys = Gp_0_inv*Gp;

        #Go through all combinations of d rows, check weight and stop as soon as a subcode is found
        Cd = Combinations(range(k),d);

        for rows_pos in Cd:

            subcode_generator = Gp_sys[rows_pos,:];

            #Check support of matrix
            this_w = 0;
            for i in range(k,n):
                g_i = subcode_generator[:,i];
                num_zeros = g_i.list().count(0);
                if num_zeros < d:
                    this_w += 1;

            if this_w == (w-d):

                subcode_solution = matrix(Fq,d,n);
                for i in range(d):
                    subcode_solution[i, rows_pos[i]] = 1;

                for j in range(n):
                    subcode_solution[:,P[j]] = subcode_generator[:,j];
                ok = 1;
                del Gp;
                del Gp_0;
                del Gp_0_inv;
                del Gp_sys;
                return 1,1,subcode_solution;

        del Gp;
        del Gp_0;
        return 1,0,0;
    
########################################################################   


#Uncomment the following lines to run a quick test

#Code and algorithm parameters
#n = 40; #code length
#k = 20; #code dimension
#q = 5; #finite field size
#w = 18; #codewords weight
#d = 3;

#Fq = GF(q);
#G = random_matrix(Fq, k, n);
#H = codes.LinearCode(G).parity_check_matrix();

#N_w = binomial(n,w)*(q^d-1)^(w-d)*gauss_binomial(k,d,q)/#gauss_binomial(n,d,q);

#print("Expected num of subcodes = ",N_w*1.);


#ok2 = 0;
#num_test = 0;
#while ok2 == 0:
#    ok1, ok2, C = subcodes_lee_brickell_isd(Fq,G,n,k,w,d);
#    num_test += 1;
#    print(num_test);
