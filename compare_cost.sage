reset();

load('cost_isd.sage');
load('kmp_cost.sage');
load('our_cost.sage');


###################################################
#Params for KMP
n = 106;
r = 48;
q = 4093;

min_ell = 1; #min value of ell to consider
max_ell = 10; #max value of ell to consider


R = 1-r/n;


P = Primes(); #set of prime numbers


#Keep n and r constant, vary ell and compute corresponding value of q
for ell in range(min_ell, max_ell+1):

    #Compute costs of KMP and our algorithm
    best_u1, best_u2, best_L1, best_L1, best_num_coll, kmp_cost = kmp_cost_numerical(n,r,ell,q);
    best_d, best_w, best_w1, best_u, best_isd, new_cost = compute_new_cost(n,r,q,ell);

    #Print values
    str_params = "$("+str(n)+", "+str(r)+", "+str(q)+", "+str(ell)+")$ &";
    str_kmp = "$2^{"+str(kmp_cost)+"}$ &";
    str_our = "$("+str(best_d)+", "+str(best_w)+", "+str(best_w1)+", "+str(best_u)+")$ & $2^{"+str(new_cost)+"}$\\\\"
    print(str_params+str_kmp+str_our);


    #Find value of q for next instance (i.e., the one with ell = ell+1)
    q = round((n/exp(1.))^(1/((ell+1)*(1-R)))*(1+1/(2*ell*(1-R)*n)*log2(6.28*n)));
    if (q in P)==False:
        q = P.next(q);

    while (factorial(n)>(q^((ell+1)*r))):
        q = P.next(q);
