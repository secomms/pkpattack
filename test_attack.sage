reset();

load("list_sorting.sage");

def log2(x):
    return log(x*1.)/log(2.);
               
############################################



##Lee and Brickell ISD algorithm
##Input:
# - Fq: finite field;
# - G: generator matrix
# - n: code length;
# - k: code dimension;
# - w: searched weight
# - Perm_set: set of length-n permutations
# - op_mode: 0 if you want only codewords with weight w, 1 if you want weight <= w

#The algorithm outputs three values:
# - first value: 0 if RREF is not possible, 1 if it has been computed
# - second value: 0 if no codeword has been found, 1 if ISD has found a codeword with the desired weight
# - third value: found codeword, or 0 (in case ISD has failed)

def lee_brickell_isd_Fq(Fq,G,n,k,w,Perm_set,op_mode):


    Fq_list = Fq.list();
    q = len(Fq_list);
    P = Perm_set.random_element(); #Random permutation
    Gp = G[:,P]; #Apply permutation
    Gp_0 = Gp[:,0:k]; 
    
    #If Gp_0 is singular, report failure
    if rank(Gp_0) < k:
        del Gp;
        del Gp_0;
        return 0,0,0;
    else:
        
        #Gp_0 is non singular, compute RREF
        Gp_0_inv = Gp_0^-1;
        new_G_1 = Gp_0_inv*Gp[:,k:n];
        
        flag_weight = 0; #flag_weight = 1 when we find codewor with desired weight
        
        c = vector(GF(2),n);
        
        #Go through all vectors with weight 2 and first set entry equal to 1
        #We directly compute the associated linear combination
        for i0 in range(0,k-1):
            
            g0 = new_G_1[i0,:];
            for i1 in range(i0+1,k):
                g1 = new_G_1[i1,:];
                
                #Compute all codewords in the form g0 + u*g1 (with u \neq 0)
                for i in range(1,q):
                    u = Fq_list[i];
                    cw = g0 + u*g1;
                    
                    this_w = len(cw.support())+2; #Hamming weight of the found codeword
                    if this_w <= w:
                        
                        
                        if op_mode == 0:
                            if this_w == w:
                                
                                #Compute full codeword and reverse permutation
                                full_cw = vector(Fq,n);
                                full_cw[i0] = 1;
                                full_cw[i1] = u;
                                
                                for j in range(n-k):                                    
                                    full_cw[k+j] = cw[0,j];
                                output_cw = vector(Fq,n);
                                for j in range(n):
                                    output_cw[P[j]] = full_cw[j];
                                    
                                    
                                ok = 1; #ok = 1, ISD has been successful
                                del Gp;
                                del Gp_0;
                                del Gp_0_inv;
                                del new_G_1;
                                return 1,1,output_cw;
                        else:
                            
                            #Compute full codeword and reverse permutation
                            full_cw = vector(Fq,n);
                            full_cw[i0] = 1;
                            full_cw[i1] = u;
                            
                            for j in range(n-k):
                                full_cw[k+j] = cw[0,j];
                            output_cw = vector(Fq,n);
                            for j in range(n):
                                output_cw[P[j]] = full_cw[j];

                            del Gp;
                            del Gp_0;
                            del Gp_0_inv;
                            del new_G_1;
                            return 1,1,output_cw;    
                        
        #ISD has failed: no valid codeword has been found                
        del Gp;
        del Gp_0;
        return 1,0,0;


##################################################################

#Prepare lists with entries picked from all_entries; the number of considered elements is specified by num_dispositions
def prepare_lists(Zq, all_entries, num_dispositions, sub_Htr, target_s, sign_val):

    list_values = [];
    list_sums = [];
    
    entries_list = Combinations(all_entries.list(),num_dispositions);
    
    for values in entries_list:
        all_dispositions = Permutations(values);
        for u in all_dispositions:
            u = matrix(vector(u));
            x = target_s + sign_val*u*sub_Htr;
            list_sums.append(x);
            list_values.append(u);

    return list_values, list_sums;

################################################

#Merge two lists, using the indexes of collisions provided in list_indexes
def merge_lists(Zq, list_a, list_b, list_indexes, num_distinct):
    
    
    #Construct merged lists and remove duplicates
    list_values = [];
    list_sums = [];
    
    for pos in list_indexes:
        i_a = pos[0];
        i_b = pos[1];
        c_a = list_a[i_a];
        c_b = list_b[i_b];        
        
        
        #check that there are not repeated entries    
        values = Set(c_a.list()).union(Set(c_b.list())).difference(Set([0]));

        if len(values)==num_distinct:
            
            new_c = block_matrix(Zq, 1,2, [c_a, c_b]);
            list_values.append(new_c);
    
    return list_values;


##################################################

##Generate a random pkp instance over Zq. The output public vector has all distinct elements

def random_pkp_instance(q,n,m):
    
    Zq = Integers(q);
    
    #Generate random matrix A with m rows, n columns and full rank
    rank_A = 0;
    while rank_A < m:
        A = random_matrix(Zq,m,n);
        rank_A = rank(A);
        
    
    #Find kernel element with all distinct elements
    #Use A as parity-check matrix of code C, compute generator for C and find codeword with distinct entries
    #Exploit information set to speed-up the search
    
    k = n-m;
    
    Zq_list = Zq.list();
    
    C_dual = codes.LinearCode(A); #code generated by A
    G_A = C_dual.parity_check_matrix(); #generator matrix for code whose parity-check is A
    
    #Check if there G_A can be put into systematic form, otherwise apply random permutation
    while rank(G_A[:,0:k])< k:
        P = Permutations(range(1,n+1)).random_element().to_matrix();
        A = A*P;
        G_A = G_A*P;
    
    
    G_A_sys = G_A[:,0:k]^-1*G_A;
    
    flag_c = 0; #flag_c = 1 when codeword with distinct entries is found
    while flag_c == 0:
        
        #Generate u with k distinct entries
        P = Permutations(range(n)).random_element();
        u = matrix(Zq,1,k); 
        for i in range(k):
            u[0,i] = Zq_list[P[i]];
        
        #Compute codeword corresponding to u
        c = u*G_A_sys;
        
        #Check if u has all distinct entries
        unique_values = uniq(c.list());
        if len(unique_values) == n:
            flag_c = 1;
            
    ##Apply random permutations
    secret_P = Permutations(range(n)).random_element();
    c_perm = vector(Zq,n);
    for i in range(n):
        c_perm[i] = c[0,secret_P[i]];
    
    
    return A, c, c_perm, secret_P;
    
    
#################################################

##Parameters for the PKP instance
n = 27;
m = 22;
q = 127;


#Parameters for the attack
w = 4;
w1 = 3;
w2 = 1;
ell = 2;


#Sanity-check on the parameters; if the expected number of solutions is higher than 1, raise a warning
num_perm = factorial(n);
pr_kernel = q^(-m);
if num_perm*pr_kernel>1:
    print("Wait, you may have more than one solution...");


###########################
##Generate random PKP instance

Zq = Integers(q); #ring of integers mod q

A, c, challenge_c, secret_P = random_pkp_instance(q,n,m);  #random instance; c is the un-permuted vector, challenge_c is the permuted one; secret_P is the secret permutation

#Enrich A with the row formed by all ones and compute s
r = m+1; 
H = block_matrix(Zq,2,1,[A,ones_matrix(Zq,1,n)]);

s = matrix(Zq,1,r);
s[0,r-1] = sum(vector(challenge_c)); #set sum of entries as the last coordinate of s

    


#Call ISD to find codeword with weight w
print("Calling ISD to find codeword with weight "+str(w)+"...");
ok2 = 0;
num_attempts = 0;
while ok2 == 0:
    num_attempts+=1;
#    print("num ISD calls = "+str(num_attempts));
    ok1, ok2, cw = lee_brickell_isd_Fq(Zq,H,n,r,w,Permutations(range(n)),1);    

#Apply first transformation (permutation and change of basis) on H and s 
first_pos = cw.support();
rem_pos = Set(range(n)).difference(Set(first_pos)).list();
perm = rem_pos[0:n-r+ell-w]+first_pos+rem_pos[n-r+ell-w:n-w];

H_perm = matrix(Zq,r,n);
for i in range(n):
    H_perm[:,i] = H[:,perm[i]];
    
cw_perm = matrix(Zq,1,n);
for i in range(n):
    cw_perm[0,i] = cw[perm[i]];

c_perm = matrix(Zq,1,n);
for i in range(n):
    c_perm[0,i] = c[0,perm[i]];
    
B = matrix(cw_perm)[0,0:r]*H_perm[:,0:r]^-1;
    
target_s = s*B.transpose();

#Start preparing lists K1 and K2; every list is represented by a pair of values and corresponding sub-syndromes through H1tr and H2tr
cw_for_lists = matrix(cw)[0,first_pos];
sub_H1tr = cw_for_lists[:, 0:w1].transpose();
sub_H2tr = cw_for_lists[0, w1:w].transpose();

print("---> Preparing initial lists K1 and K2...");
list_values_1, list_sums_1 = prepare_lists(Zq, challenge_c, w1, sub_H1tr, matrix(Zq,1,1), 1);
list_values_2, list_sums_2 = prepare_lists(Zq, challenge_c, w2, sub_H2tr, target_s, -1);


#Find collisions 
indexes_collisions = colliding_indexes(list_sums_1,list_sums_2); #indexes of collisions for lists 1 and 2


##Check sizes of the obtained lists, and number of collisions (compare with theoretical values)
emp_list_sizes = [log2(len(list_values_1)), log2(len(list_values_2))];
th_list_sizes = [log2( factorial(n)/factorial(n-w1)), log2(factorial(n)/factorial(n-w2) )];

print("Emp initial sizes (in log2) = "+str(emp_list_sizes));
print("Th theoretical sizes  (in log2) = "+str(th_list_sizes));

print("Emp num_collisions (in log2) = "+str(log2(len(indexes_collisions))));
print("Th num collisions (in log2) = "+str(th_list_sizes[0]+th_list_sizes[1]-log2(q)));

#Merge the lists to remove badly constructed elements
list_values_1 = merge_lists(Zq, list_values_1, list_values_2, indexes_collisions, w1+w2);

print("Emp size of produced list (in log2) = "+str(log2(len(list_values_1))));
print("Th size of produced list  (in log2) = "+str(log2( factorial(n)/factorial(n-w)/q )));


#Put H_perm in systematic form, apply same transformation to s

B = H_perm[:,n-r:n]^-1;
new_H_perm = B*H_perm;
new_s = s*B.transpose();
new_target_s = new_s[0,1:ell];

#Prepare lists L1 and L2
sub_H1tr = new_H_perm[1:ell, 0:n-r+ell-w].transpose();
sub_H2tr = new_H_perm[1:ell, n-r+ell-w:n-r+ell].transpose();

print("---> Preparing lists L1 and L2...");

#Use the already found values to prepare the second list
list_values_2 = copy(list_values_1);
list_sums_2 = [];
for x in list_values_2:
    list_sums_2.append(new_target_s - x*sub_H2tr);

list_values_1, list_sums_1 = prepare_lists(Zq, challenge_c, n-r+ell-w, sub_H1tr, matrix(Zq,1,ell-1), 1);


#Find collisions and then merge
indexes_collisions = colliding_indexes(list_sums_1,list_sums_2); #indexes of collisions for lists 1 and 2


##Check sizes of the obtained lists, and number of collisions (compare with theoretical values)
emp_list_sizes = [log2(len(list_values_1)), log2(len(list_values_2))];
th_list_sizes = [log2( factorial(n)/factorial(n-(n-r+ell-w))), log2( factorial(n)/factorial(n-w)/q )];

print("Emp sizes (in log2) = "+str(emp_list_sizes));
print("Th theoretical sizes  (in log2) = "+str(th_list_sizes));

print("Emp num_collisions (in log2) = "+str(log2(len(indexes_collisions))));
print("Th num collisions (in log2) = "+str(th_list_sizes[0]+th_list_sizes[1]-(ell-1)*log2(q)));


list_values = merge_lists(Zq, list_values_1, list_values_2, indexes_collisions, n-r+ell);

print("Emp size of produced list (in log2) = "+str(log2(len(list_values))));
print("Th size of produced list  (in log2) = "+str(log2( factorial(n)/factorial(n-(n-r+ell))/q^ell )));

#Check validity of the found candidates
print("---> Checking final candidates...");

#G = codes.LinearCode(new_H_perm).parity_check_matrix(); 
flag_sys = 0;
while flag_sys == 0:
    info_set = Combinations(range(n-r+ell),n-r).random_element();
    J = Set(range(n)).difference(Set(info_set)).list();
    sub_H = new_H_perm[:,J];
    if rank(sub_H)==r:
        flag_sys = 1;

inv_sub_H = sub_H^-1;        
H_sys = inv_sub_H*new_H_perm;        
target_entries = uniq(challenge_c.list());
check_s = new_s*inv_sub_H.transpose();

for u in list_values:
    
    remaining_entries = check_s - u[:,info_set]*H_sys[:,info_set].transpose();
    values = uniq(u.list()+remaining_entries.list());
   
    if values == target_entries:
        print("!!!! We found a solution !!!!");
    

