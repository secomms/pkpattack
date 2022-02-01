reset();

load('cost_isd.sage');

###################################################

def gauss_binomial(m,r,q):
   
    x = 1.;
    for i in range(r):
        x = x*(1.-q^(m-i*1.))/(1.-q^(1.+i));
   
    return x;

##################################à

def new_cost(n,m,q):
   
   
    best_cost = 10000000000000000;
   
    best_d = 0;
    best_w = 0;
    best_w1 = 0;
    best_ell = 0;
    
    best_isd = 0;
    
    for w in range(1,n):
        
        for d in range(1,w):
            
            N_w = binomial(n,w)*(q^d-1)^(w-d)*gauss_binomial(m,d,q)/gauss_binomial(n,d,q);
            
            if N_w>1:
                    
                c_isd = cost_isd(q,n,m,d,w,N_w);
                if c_isd.imag()==0:
                    for w1 in range(1,w):

                        w2 = w-w1;
                        K1 = log2(binomial(n,w1)*factorial(w1));
                        K2 = log2(binomial(n,w2)*factorial(w2));

                        if (K1<best_cost)&(K2<best_cost):

                            K1_merge_K2 = K1+K2-d*log2(q);
                            
                            L2 = log2(factorial(n)/factorial(n-w))-d*log2(q);

                            for ell in range(1,m):

                                pos1 = n-m+ell-w;
                            
                                L1 = log2(factorial(n)/factorial(n-pos1))
                                L1_merge_L2 = L1 + L2 - (ell-d)*log2(q);
                                L = log2(factorial(n)/factorial(m-ell)/(q^ell));

                                cost = log2(2^c_isd + 2^K1 + 2^K2 + 2^K1_merge_K2 + 2^L1 + 2^L2+2^L1_merge_L2 +2^L);

                                if cost < best_cost:

                                    best_cost = cost;

                                    best_d = d;
                                    best_w = w;
                                    best_w1 = w1;
                                    best_ell = ell;
                                    best_isd = c_isd;
                                    
                                 #   print("cost = ",cost);
                                 #   print("[d, w, w1, ell] = ",[d,w,w1,ell]);
                                 #   print("[Nw, C_isd] = ",[log2(N_w),c_isd]);
                                 #   print("[K1, K2, merge] = ",[K1,K2,K1_merge_K2]);
                                 #   print("[L1, L2, merge] = ",[L1,L2,L1_merge_L2]);
                                 #   print("[L] = ",[L]);
                                 #   print("-----------------------");
                                    
    return best_d, best_w, best_w1, best_ell, best_isd, best_cost;

   
#######################################   
   

n = 69;
m = 41;
q = 251;

best_d, best_w, best_w1, best_ell, c_isd, c1 = new_cost(n,m+1,q);

print("[n,m,q]",[n,m,q]);
print( "ATTACK: [d, w, w1, ell, c_isd, cost] = "+str([best_d, best_w, best_w1, best_ell, c_isd, c1]) );
print("-----------------------------------------------------");


n = 94;
m = 54;
q = 509;

best_d, best_w, best_w1, best_ell, c_isd, c1 = new_cost(n,m+1,q);

print("[n,m,q]",[n,m,q]);
print( "ATTACK: [d, w, w1, ell, c_isd, cost] = "+str([best_d, best_w, best_w1, best_ell, c_isd, c1]) );
print("-----------------------------------------------------");

n = 106;
m = 47;
q = 4093;

best_d, best_w, best_w1, best_ell, c_isd, c1 = new_cost(n,m+1,q);

print("[n,m,q]",[n,m,q]);
print( "ATTACK: [d, w, w1, ell, c_isd, cost] = "+str([best_d, best_w, best_w1, best_ell, c_isd, c1]) );
print("-----------------------------------------------------");
