def log2(x):
    return log(x*1.)/log(2.);


###################################

q = 251
n = 198;
k = 94;
w = 70;

def peters_isd(q,n,k,w):

    x = floor(k/2);


    log2q=log2(q);
    mincost=10000000;
    bestp=0;
    bestl=0;
    for p in range(1,11):
        Anum=binomial(x,p);
        Bnum=binomial(k-x,p);
        #for(l=1,floor(log(Anum)/log(q)+p*log(q-1)/log(q))+10,\
        for l in range(1,floor( log(Anum)/log(q)+p*log(q-1)/log(q))+10 +1):
        # if(q==2):
            ops=0.5*(n-k)^2*(n+k)+ ((0.5*k-p+1)+(Anum+Bnum)*(q-1)^p)*l+ q/(q-1.)*(w-2*p+1)*2*p*(1+(q-2)/(q-1.))*Anum*Bnum*(q-1)^(2*p)/q^l;
        #ops=(n-k)^2*(n+k)\
        #     + ((0.5*k-p+1)+(Anum+Bnum)*(q-1)^p)*l\
        #     + q/(q-1.)*(w-2*p+1)*2*p*(1+(q-2)/(q-1.))*\
        #       Anum*Bnum*(q-1)^(2*p)/q^l;
    
            prob=Anum*Bnum*binomial(n-k-l,w-2*p)/binomial(n,w);
            cost=log2(ops)+log2(log2q)-log2(prob);
            if cost<mincost:
                mincost=cost;
                bestp=p; bestl=l;
    
#    cost=mincost;
#    p=bestp;
#    l=bestl;
#    print("Given q=",q,", n=",n,", k=",k,", w=",w);
#    print("parameters p=",p, " and l=",l," yield 2^",cost," bit ops");

    return mincost;


#########################################

def beullens_isd(q,n,k,d,w):
    
    success_pr = binomial(w,d)*binomial(n-w,k-d)/binomial(n,k);
    c_iter = k^3+binomial(k,d);
    
    return log2(c_iter)-log2(success_pr);

#########################################

def cost_isd(q,n,k,d,w,Nw):
    
    if d == 1:
        
        c_isd = peters_isd(q,n,k,w);
        
        if c_isd.imag()!=0:
           
            pr_isd = binomial(w,1)*binomial(n-w,k-1)/binomial(n,k);
            cost = k^3+binomial(k,1);
            
            if (k-2)<(n-w):
                pr_isd+= binomial(w,2)*binomial(n-w,k-2)/binomial(n,k);
                cost += binomial(k,2)*(q-1);
            
            if pr_isd == 0:
                c_isd = 10000000000;
            else:
                c_isd = log2(cost/pr_isd);
    
    else:
        
        c_isd = beullens_isd(q,n,k,d,w);
    
    return max(0,c_isd-log2(Nw));
        
