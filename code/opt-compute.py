import math,random
from math import log,ceil,floor,sqrt

w_value=0

## Compute binom(n,k) in big integers
## using a precomputed table if possible
def binom(n,k):
    if binoms[n][k]==0:
        value = factorials[n]//(factorials[k]*factorials[n-k])
        binoms[n][k]=value
        return value
    else:
        return binoms[n][k]

## Compute the layer-@k of hypercube of dimension @n and range [0;@m]
## Equals the coefficient of x^k in the product (1+x+x^2+... +x^{m})^n
## Requires n>1 and k<=n*m
def nb(k,m,n):
    sum=0
    for s in range(0,1+floor(k/(m+1))):
        summand = binom(n,s)*binom(k-s*(m+1)+n-1,n-1)
        if(s&1==0):
            sum += summand
        else:
            sum -= summand
    return sum

def test_precompute():
    precompute_real(10,11)
    print(nb(50,10,10))
    print(layer_sizes[50])

## Precomputing the hypercube layer sizes. Hypercube consists of
## @dimension elements each taking values from 0 to val-1
def precompute_real(dimension,val):
    global powcc
    powcc=1   #compute the size of the hypercube
    for i in range(dimension):
        powcc*=val
    global max_distance
    max_distance = (val-1)*dimension  #distance between all-0 and all-val-1
    global factorials;  #array of factorials
    factorials = [0 for i in range(max_distance+dimension+1)] #some margin
    factorials[0]=1
    for i in range(1,max_distance+dimension):
        factorials[i] = factorials[i-1]*i
    global binoms;   #double array of binoms
    binoms = [[0 for i  in range(max_distance+dimension)] for j in range(max_distance+dimension)]
    global layer_sizes
    layer_sizes =[0 for i in range(max_distance+1)]
    for i in range(max_distance+1):
        layer_sizes[i] = nb(i,val-1,dimension)
    global w_value
    w_value = val

## Verification cost of layer d
## Should be modified if the checksum is counted
def compute_cost(d):
    return d

## checksum cost
def compute_checksum_cost(sum_value):
    total_cost=0
    while(sum_value !=0):
        global w_value
        total_cost = sum_value % w_value
        sum_value //=w_value
    return total_cost



#- \frac{\sum_d |L_d|C_d}{2^v}\pm\sqrt{ \frac{(\sum_d |L_d|C_d)^2}{2^{2v}}-
#  \frac{\mu (\sum_d |L_d|C_d)^2-\sum_d |L_d|C_d^2}{\mu 2^{2v}-2^v}}

# @mu_scaled = mu_(l2)*(powcc) ;
def compute_extrema(mu_scaled,d0):
    # 1. Computing sums
    sum_ld=0   # sum_d L_v(d)
    sum_ld_cd=0  # sum_d L_v(d)*C_d
    sum_ld_cd2=0 # sum_d L_v(d)*C_d^2
    for d in range(d0+1):
        cd = compute_cost(d)
        sum_ld += layer_sizes[d]
        sum_ld_cd += layer_sizes[d]*cd
        sum_ld_cd2 += layer_sizes[d]*cd*cd

    # 2. Computing lambda_2
    #  Discriminant = v1 - v2/v3
    v1 = (sum_ld_cd*sum_ld_cd)/(sum_ld*sum_ld)
    v2=(mu_scaled*sum_ld_cd*sum_ld_cd)-sum_ld_cd2*powcc #multiply right term by powcc since mu is scaled
    v3 = (mu_scaled*sum_ld- powcc)*sum_ld #multiply right term by powcc since mu is scaled
    if(v3==0):
        return(0,0)
    v4 = v2/v3
    discr = v1-v4
    if(discr<0):
        #print("Discriminant is negative")
        return (0,0)
    lambda_21 = sum_ld_cd/sum_ld+sqrt(discr)  # first lambda_2
    #lambda_22 = -sum_ld_cd/sum_ld-sqrt(discr)  # second lambda_2


    #t1 = root1*root1 + 2*root1*sum_ld_cd/sum_ld+v4   #for debug
    #t2 = root2*root2 + 2*root2*sum_ld_cd/sum_ld+v4   #for debug
    #print("Lambda 2 variants: ", root1, root2)

    # 3. Compute two variants for lambda_1.  We compute it downscaled by sum_ld= sum_d L_v(d)
    lambda_11_downscaled = lambda_21/2-sum_ld_cd/(2*sum_ld)
    if(lambda_11_downscaled ==0):
        return (0,0)
    #lambda_12_downscaled = -lambda_22/2-sum_ld_cd/(2*sum_ld)

    # 4. Compute two variants for mu_d , upscaled by sum_ld
    opt_mu1_scaled = [0 for i in range(max_distance)] #upscaled by sum_ld
    opt_mu2_scaled = [0 for i in range(max_distance)]
    full_cost1=0   #total cost, variant 1
    checksum_cost=0 #total cost with checksum into account
    full_cost2=0   #total cost, variant 2
    var1_neg = False  #flag that one of mu_d values (variant 1) is negative
    var2_neg = False  #flag that one of mu_d values (variant 2) is negative
    #test_v1=0  #checking summation
    #test_v2 = 0
    #test2_v1=0 #checking quadratic
    #test2_v2 = 0
    for i in range(d0+1):
        opt_mu1_scaled[i] = (lambda_21-compute_cost(i))/(2*lambda_11_downscaled)
        #test_v1 += layer_sizes[i]*opt_mu1_scaled[i]
        #test2_v1 += layer_sizes[i]*opt_mu1_scaled[i]*opt_mu1_scaled[i]
        if(opt_mu1_scaled[i]<0):
            var1_neg = True
        full_cost1 +=  opt_mu1_scaled[i]*(layer_sizes[i]*compute_cost(i))  #the cost is temporarily scaled up by sum_ld
        checksum_cost +=  opt_mu1_scaled[i]*layer_sizes[i]*(compute_cost(i)+compute_checksum_cost(d0-i))  #the cost is temporarily scaled up by sum_ld
        #opt_mu2_scaled[i] = -(lambda_22+compute_cost(i))/(2*lambda_12_downscaled)
        #test_v2 += layer_sizes[i]*opt_mu2_scaled[i]
        #test2_v2 += layer_sizes[i]*opt_mu2_scaled[i]*opt_mu2_scaled[i]
        #if(opt_mu2_scaled[i]<0):
        #    var2_neg = True
        #full_cost2 +=  opt_mu2_scaled[i]*layer_sizes[i]*compute_cost(i)
    full_cost1 /= sum_ld
    checksum_cost /= sum_ld
    #full_cost2 /= sum_ld
    #test_v1 /= sum_ld  #should be equal to 1
    #test_v2 /= sum_ld   #should be equal to 1
    #test2_v1 /= (sum_ld*sum_ld/powcc) #should be equal to mu_scaled = mu/2^l
    #test2_v2 /= (sum_ld*sum_ld/powcc) #should be equal to mu_scaled = mu/2^l
    #a1 = log(test2_v1,2)
    #a2 = log(test2_v2,2)
    #a3 = log(mu_scaled,2) #should be equal to a2 and a1
    if(var1_neg):
        full_cost1=0
    #if(var2_neg):
    #    full_cost2=0
    #if(full_cost1!=0 or full_cost2!=0):
    #    print("opt at 0: ",opt_mu2_scaled[0]*pow(2,l)/sum_ld, opt_mu2_scaled[0]*pow(2,l)/sum_ld/threshold)
    #    print("opt_scaled at 0: ",opt_mu2_scaled[0])
    return (full_cost1,checksum_cost)

## Computes the lower bound on the verification cost for L2 value @mu_scaled
## Hypercube layers must be precomputed with `precompute_real()`
def compute_L2_lower_bound(chains,chain_len,mu_scaled):
    precompute_real(chains,chain_len) #compute layer sizes
    min_cost=100000000
    min_cost_checksum=100000000
    opt_d=0
    opt_d_checksum=0
    for d0 in range(1,max_distance):
        (cost1,cost_checksum) = compute_extrema(mu_scaled,d0)
        if(cost1<min_cost and cost1>0):
            opt_d = d0
            min_cost = cost1
        if(cost_checksum<min_cost_checksum and cost_checksum>0 and cost1>0):
            opt_d_checksum = d0
            min_cost_checksum = cost_checksum
    return (min_cost,opt_d,min_cost_checksum,opt_d_checksum)

## Computes the lower bound on the verification cost for max entropy value @mu_scaled
## Hypercube layers must be precomputed with `precompute_real()`
# Algorithm: find max d0 such that sum_d<d_0 Ld*mu_scaled/val^dimension <1
# the cost_scaled (multiplied by val^dimension) is then sum_d<d_0 Ld*mu_scaled*Cd
def compute_max_entropy_lower_bound(chains,chain_len,mu_scaled):
    precompute_real(chains,chain_len) #compute layer sizes
    min_cost=100000
    sum_ld=0
    sum_ld_cd=0
    for d in range(max_distance):
        sum_ld += layer_sizes[d]
        sum_ld_cd += layer_sizes[d]*compute_cost(d)
        if((sum_ld+layer_sizes[d+1])*mu_scaled>powcc):  #we will go over in the next step, so d+1 = d_0
            mu_d0_scaled = (powcc-mu_scaled*(sum_ld))/layer_sizes[d+1];
            cost_scaled = sum_ld_cd*mu_scaled+(powcc-mu_scaled*(sum_ld))*compute_cost(d+1)
            cost = cost_scaled/powcc
            return (cost,d+1)
    return (min_cost,0)

def compute_slice_size(d0,slice_len):
    sum_ld=0
    for i in range(d0-slice_len+1,d0+1):
        sum_ld += layer_sizes[i]
    return sum_ld

def print_optimal_L2_values(sec_level,min_v,max_v,min_w,max_w):
    fname = "opt-L2-sec"+str(sec_level)+".txt"
    for v in range(min_v,max_v+1): #we need signatures between 0.5K and 5K assuming 28 bytes per chain
        chain_len_winternitz=ceil(pow(2,sec_level/v))
        wint_cost = v*(chain_len_winternitz-1)/2
        print("Chains=",v," Winternitz cost w/0 checksum=",wint_cost,file=open(fname,"a+"))
        best_cost = wint_cost
        best_w = min_w
        for w in range(min_w,max_w+1): #no need in too expensive verifiers
            global w_value
            w_value = w
            hcube_size_log2 = v*log(w,2)
            if(hcube_size_log2<sec_level):
                continue
            mu_scaled = ceil(pow(2,hcube_size_log2-sec_level)) #the loss we can tolerate
            (c,d0,cc,dc) = compute_L2_lower_bound(v,w,mu_scaled)  #compute lower bound and the d0 value
            if(c<best_cost-0.5):
                best_cost = c
                best_w = min_w
                print("w=",w, " L2 lower bound=","{0:0.2f}".format(c)," layer=",d0,
                      " mu_log2=","{0:0.2f}".format(hcube_size_log2-sec_level),file=open(fname,"a+"))#Th Ver cost w/o checksum:  683.23  d_0= 775  bit loss= 0.4 chain len= 54
                #Pract Ver cost with 1-chain checksum:  774.00  bit loss= 2.1

def checksum_cost(layer,max_layer,w):
    inverted_checksum_value = max_layer-layer
    cost =0
    while inverted_checksum_value>0:
        cost += inverted_checksum_value % w
        inverted_checksum_value //= w
    return cost

def compute_cost_winternitz(w,v):
    precompute_real(v,w)
    cost=0
    for d in range((w-1)*v+1):
        cost += (d+checksum_cost(d,(w-1)*v,w))*layer_sizes[d]/powcc;
    return cost

def find_best_wots(sec_level,v):
    v_base=v-1
    w=2
    while(v_base>0):
        w=ceil(pow(2,sec_level/v_base))
        checksum_len = ceil(log((w-1)*v_base+1,w))
        if(checksum_len<= v- v_base):
            break
        else:
            v_base -=1
    cost = compute_cost_winternitz(w,v_base)
    return (cost,v_base,w)

def find_best_wots_ts(sec_level,v):
    w=ceil(pow(2,sec_level/v))
    middle_layer = (w-1)*v/2
    return (middle_layer,w)

def find_best_wots_cs(sec_level,v):
    w=2
    while(True):
        precompute_real(v,w)
        middle_layer = (w-1)*v//2
        if(log(layer_sizes[middle_layer],2)>=sec_level):
            return(middle_layer,w)
        else:
            w+=1

def find_best_single_layer(sec_level,v,min_w,max_w):
    best_cost = 100000
    best_w = 0
    for w in range(min_w,max_w+1):
        precompute_real(v,w)
        for d0 in range((w-1)*v):
            if(log(layer_sizes[d0],2)>=sec_level):
                if(d0<best_cost):
                    best_cost = d0
                    best_w = w
    return (best_cost,best_w)


def find_best_small_top(sec_level,v,min_w,max_w):
    best_cost = 100000
    best_w = 0
    for w in range(min_w,max_w+1):
        precompute_real(v-1,w)
        sum_ld=0
        for d0 in range((w-1)*(v-1)):
            sum_ld += layer_sizes[d0]
            if(log(sum_ld,2)>=sec_level):
                if(d0<best_cost):
                    best_cost = d0
                    best_w = w
    return (best_cost,best_w)


def find_opt(sec_level,v,min_w,max_w):
    best_cost = 100000
    best_w = 0
    for w in range(min_w,max_w+1):
        hcube_size_log2 = v*log(w,2)
        if(hcube_size_log2<sec_level):
            continue
        mu_scaled = ceil(pow(2,hcube_size_log2-sec_level))
        precompute_real(v,w)
        (c,d0,cc,dc) = compute_L2_lower_bound(v,w,mu_scaled)  #compute lower bound and the d0 value
        if(c<best_cost-0.2):
            best_cost = c
            best_w = w
    return (best_cost,best_w)


def find_best_full_top(sec_level,v,min_w,max_w):
    best_cost = 100000
    best_w = 0
    v_base=v-1
    best_v = v-2
    while(v_base>v-5):
        for w in range(min_w,max_w+1):
            hcube_size_log2 = v_base*log(w,2)
            if(hcube_size_log2<sec_level):
                continue
            mu_scaled = ceil(pow(2,hcube_size_log2-sec_level))
            precompute_real(v_base,w)
            (c,d0,cc,dc) = compute_L2_lower_bound(v_base,w,mu_scaled)  #compute lower bound and the d0 value
            checksum_len = ceil(log(dc+1,w))
            if(checksum_len!= v- v_base):
                continue
            if(cc<best_cost):
                best_cost = cc
                best_w = w
                best_v = v_base
        v_base -=1
    return (best_cost,best_w, best_v)

def find_best_constructions(sec_level,v,min_w,max_w):
    print("Optimal signatures of ",v," chains")
    #1. WOTS
    (cost_wots,v_wots,w_wots)  = find_best_wots(sec_level,v)
    print("WOTS: cost=","{0:0.2f}".format(cost_wots), " v=",v_wots," w=",w_wots)
    #1.2 WOTS-CS
    (cost_cs,w_cs)  = find_best_wots_cs(sec_level,v)
    print("WOTS-CS: cost=","{0:0.2f}".format(cost_cs), " w=",w_cs)
    #1.3 WOTS-TS
    (cost_tsw,w_tsw)  = find_best_wots_ts(sec_level,v)
    print("WOTS-TS: cost=","{0:0.2f}".format(cost_tsw), " w=",w_tsw)
    # 2. TLFC
    (cost_ft,best_ft, best_vt) = find_best_full_top(sec_level,v,min_w,max_w)
    print("TLFC: cost=","{0:0.2f}".format(cost_ft), " w=",best_ft, " v=",best_vt)
    # 3. TL1C
    (cost_st,best_st) = find_best_small_top(sec_level,v,min_w,max_w)
    print("TL1C: cost=","{0:0.2f}".format(cost_st), " w=",best_st)
    # 4. TSL
    (cost_single,best_single) = find_best_single_layer(sec_level,v,min_w,max_w)
    print("TSL: cost=","{0:0.2f}".format(cost_single), " w=",best_single)
    # 5. Lower bound
    (opt,w_opt) = find_opt(sec_level,v,min_w,max_w)
    print("Optimal: cost=","{0:0.2f}".format(opt), " w=",w_opt)

## Compute cost of some variants
def find_constructions(sec_level,v,min_w,max_w):
    #1. Winternitz
    # compute checksum chains
    chain_len_winternitz=ceil(pow(2,sec_level/v))
    wint_cost = v*(chain_len_winternitz-1)/2
    wint_checksum = ceil(log(v*chain_len_winternitz+1,chain_len_winternitz))
    print("Chains=",v," Winternitz cost w/0 checksum=",wint_cost, " checksum chains=",wint_checksum)
    best_cost = wint_cost
    best_w = min_w
    for w in range(min_w,max_w+1): #no need in too expensive verifiers
        global w_value
        w_value=w
        hcube_size_log2 = v*log(w,2)
        if(hcube_size_log2<sec_level):
            continue
        mu_scaled = ceil(pow(2,hcube_size_log2-sec_level))
        (c,d0,cc,dc) = compute_L2_lower_bound(v,w,mu_scaled)  #compute lower bound and the d0 value
        if(c<best_cost-0.1):
            best_cost = c
            best_w = min_w
            print("w=",w, " L2 lower bound=","{0:0.2f}".format(c)," layer=",d0,
                    " mu_log2=","{0:0.2f}".format(hcube_size_log2-sec_level))
            # 1. Full checksum
            ## 1.1 Cost is sum_d mu_d*ell_d*(d+checksum_cost(d))
            print(" checksum cost=","{0:0.2f}".format(cc)," layer=",dc, " chains=", ceil(log(dc+1,w)))
            # 2. Uniformly on first d0 layers
            sum_ld =0
            for d0 in range((w-1)*v):
                sum_ld += layer_sizes[d0]
                if(log(sum_ld,2)>=sec_level):
                    print("1-chain checksum cost=",d0)
                    break
            # 3. On one layer
            for d0 in range((w-1)*v):
                if(log(layer_sizes[d0],2)>=sec_level):
                    print("no checksum cost=",d0)
                    break



## Example of an efficient instance: hypercube with dimension 76 and range 6 of size 2^(196.4571)
## When mapping message to the top 86 layers, which constitute the 2^{-41} fraction of whole hypercube
## We obtain the average cost of 85.6 hashes and size 2212 bytes
def test_tcr():
    chains = 76
    chain_len = 6
    loss = 41.4571
    (c,d0) = compute_max_entropy_lower_bound(chains,chain_len,pow(2,loss))  #compute lower bound and the d0 value
    print("Hypercube: ", chains, " chains of length ",chain_len)
    print("Cost min-entropy lower bound is ", c, " d0=",d0)
    print("Size w/o checksum ", (chains)*4*7, " bytes") #assuming 28 bytes per chain
    print("Signature size w 1-chain checksum ", (chains+1)*4*7, " bytes")
    print("Full signature cost: ", d0)
    img_size = compute_slice_size(d0,chain_len)
    print("Image size ", log(img_size,2))
    print("Layer sizes: ")
    print(d0+2, log(layer_sizes[d0+2  ],2))
    print(d0+1, log(layer_sizes[d0+1],2))
    print(d0, log(layer_sizes[d0],2))
    print(d0-1, log(layer_sizes[d0-1],2))
    print(d0-2, log(layer_sizes[d0-2],2))

def test_L2():
    chains = 76
    chain_len = 6
    loss = 41.4571
    (c,d0) = compute_L2_lower_bound(chains,chain_len,pow(2,loss))  #compute lower bound and the d0 value
    print("Hypercube: ", chains, " chains of length ",chain_len)
    print("Cost L2 lower bound is ", c, " d0=",d0)
    print("Size w/o checksum ", (chains)*4*7, " bytes") #assuming 28 bytes per chain
    print("Signature size w 1-chain checksum ", (chains+1)*4*7, " bytes")
    print("Full signature cost: ", d0)
    img_size = compute_slice_size(d0,chain_len)
    print("Image size ", log(img_size,2))
    print("Layer sizes: ")
    print(d0+2, log(layer_sizes[d0+2  ],2))
    print(d0+1, log(layer_sizes[d0+1],2))
    print(d0, log(layer_sizes[d0],2))
    print(d0-1, log(layer_sizes[d0-1],2))
    print(d0-2, log(layer_sizes[d0-2],2))

#find_best_constructions(128,64,2,16)
#find_best_constructions(128,32,2,32)

#find_best_constructions(128,136,2,8)
#find_best_constructions(128,132,2,8)
#find_best_constructions(128,128,2,8)


#find_best_constructions(128,86,2,12)
#find_best_constructions(128,84,2,12)
#find_best_constructions(128,81,2,12)

#find_best_constructions(128,68,2,16)
#find_best_constructions(128,67,2,16)
#find_best_constructions(128,64,2,16)

#find_best_constructions(128,55,2,32)
#find_best_constructions(128,50,2,32)
#find_best_constructions(128,45,2,32)
#find_best_constructions(128,40,2,32)
#find_best_constructions(128,35,2,32)
#find_best_constructions(128,30,2,48)
find_best_constructions(128,25,2,96)

#find_best_constructions(160,168,2,8)
#find_best_constructions(160,165,2,8)
#find_best_constructions(160,160,2,8)

#find_best_constructions(160,106,2,12)
#find_best_constructions(160,104,2,12)
#find_best_constructions(160,101,2,12)

#find_best_constructions(160,84,2,12)
#find_best_constructions(160,80,2,12)

#find_best_constructions(160,70,2,16)
#find_best_constructions(160,60,2,24)
#find_best_constructions(160,50,2,32)

#find_best_constructions(160,45,2,32)
#find_best_constructions(160,40,2,48)
#find_best_constructions(160,35,2,64)

#test_precompute()
