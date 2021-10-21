#########################################################
 # Projection for the Truncation Method
#########################################################

function  convertBasisE10(vE::Array{Int64,1}, E::Int64, N::Int64)
    #converts the vector vE into an integer p
    #vE is of length N, containing integers between 1 and E
    #p is the decimal represention of v (which is in basis E)
    #convention: vE(1) represents the idiosyncratic state in the present period and vE(N) represents the past idosyncratic state N-1 periods before

    if (length(vE)!=N)
        println("incompatible length for vE");
    end

    if minimum(vE)<1
        println("all elements of vE should be >= 1");
    end

    if maximum(vE)>E
        println("all elements of vE should be <= E");
    end

    p = 0;
    for k=1:N
        p = p+(vE[k]-1)*E^(N-k);
    end

    p=p+1; #otherwise between 0 and (E+1)^N - 1
    return p
end


function  convertBasis10E(p::Int64, E::Int64, N::Int64)
    #converts the integer p into a vector vE of length N, containing integers between 1 and E+1 representing p in basis E
    #convention: vE(1) represents the idiosyncratic  state in the present period and vE(N) represents the past idosyncratic state N-1 periods before

    vE=Array{Int64,1}

    if ((p<1)|(p>E^N))
        println("incompatible integer p");
    end

    ptemp = p-1;
    vE = zeros(Int64,1,N);r=0;
    for k=1:N
        r = mod(ptemp, E);
        vE[N-k+1] = round(Int,1+r);
        ptemp = floor(Int,(ptemp-r)/(E));
    end
    return vE
end


#############Truncation Method
function Projection_plan(
        N::Integer, # length of the truncation
        AA0::AiyagariModelEGM,
        polA_ss::Array{T,1},
        polC_ss::Array{T,1},
        upc::Array{T,1},
        D_ss::Array{T,1},
        Aggs_ss::AggVarsEGM{T,T},
        Mat::Array{T,2},
        A::T,
        B::T,
        vstan::Array{T,1},
        #Ls::Array{T,1},
        #Ulevel::Array{T,1},
        states::Array{T,1}) where {T <: Real}

    @unpack params,aGrid,na,dGrid,nd,ns,Transv,states = AA0
    @unpack γ,β = params

    R = Aggs_ss.R
    w = Aggs_ss.w
    Ltot = dot(vstan,states)
    Trans   = reshape(Transv,ns,ns)

    #############Computing CC Bins
    Vmat = Mat[1:na,:]
    aGrid = AA0.aGrid
    aGrid[1] = 0
    Va = repeat(aGrid,ns)'*Mat #policy rule in a a'(a)
    pol = reshape(Va,na,ns)
    Vind = zeros(Int,na*ns)
    Wind = zeros(Float64,na*ns)

    for i = 1:length(Va)
        indcc=findall(x->x<=Va[i],aGrid) #this is the index
        Vind[i] = maximum(indcc)
        Wind[i] = 1-(Va[i] - aGrid[Vind[i]])/(aGrid[Vind[i]+1]-aGrid[Vind[i]])
    end

    Ls  = fill(dot(vstan,states),ns) #vector of labor supplies, it depends on the state
    Lv = kron(Ls[1:ns],ones(na,1))
    ytype = kron(collect(StepRange(1, Int8(1), ns)),ones(Int64,na)) #on the initial grid this define the type of agents in this economy
    Cpol =  - Va' + Aggs_ss.R*repeat(aGrid,ns) + Aggs_ss.w*states[ytype].*Lv + Aggs_ss.TT*ones(Float64,na*ns)

    Xe = (Aggs_ss.w .*states)
    Xc = 0*Xe
    Xc = kron(Xc,ones(na,1))

    resE = zeros(Float64,na*ns) #this is to see whether the agents will be credit constrained in the model
    cimp = zeros(Float64,na*ns)
    for i=1:na
        for j=1:ns
            ind = i+(j-1)*na
            EUc = 0
            for jp=1:ns
                    EUc = EUc + Trans[jp,j] * (
                    (Wind[ind]*Cpol[Vind[ind]+(jp-1)*na] +
                    (1-Wind[ind])*Cpol[Vind[ind]+1+(jp-1)*na] -
                                Xc[Vind[ind]+(jp-1)*na])^-γ)
            end
            resE[ind] = (Cpol[ind] -Xc[ind] )^-γ - β*R*EUc
            cimp[ind] = (β*R*EUc)^(-1/γ)+Xc[ind]
        end
    end

    resi = upc - β*R*(Mat')*upc #this is to see whether the agents will be credit constrained. Notice the agents will not be credit constraint in case the value is equal to zero.



    Nidio = ns #this is to update the sizes of the bin according to the function we had defined previously
    Size = zeros(Float64,Nidio^N) #size of the bin
    compt=0
    Ntoi = zeros(Int,Nidio^N)
    for i=1:Nidio^N
        vE =  convertBasis10E(i, Nidio, N)
        s0 = vstan[vE[N]]
        Siz = s0
        for j=1:N-1
            Siz = Siz*Trans[vE[N-j],vE[N+1-j]]
        end
        if Siz>0
            Size[compt+1]=Siz
            Ntoi[compt+1]=i
            compt+=1
        end
    end
    Nbin = compt

    CCv = zeros(Float64,na) #to compute the fraction of credit-constrained bin
    CCv[1] = 1
    CCv[2] = 1 #(for accuracy reason)
    CCv = repeat(CCv,Nidio)  #vector of 1 only for the lowest wealth bin for each Nidio

    Ci = (Cpol .-Xc)
    Ci2 = Ci.^2
    UCi = Uvv(Ci)
    dGridl_t = repeat(dGrid,ns)
    conso   = zeros(Float64,Nbin) #consumption
    at   = zeros(Float64,Nbin) #initial of period assets
    ae   = zeros(Float64,Nbin) #this represents the policy function for total assets
    ae2   = zeros(Float64,Nbin) #this represents the policy function for total assets
    uc   = zeros(Float64,Nbin) #marginal utility of the bin
    Euc  = zeros(Float64,Nbin) #expected utility of the bin
    ECib  = zeros(Float64,Nbin) #expected utility of the bin
    ECi2b  = zeros(Float64,Nbin) #expected utility of the bin
    EUCib  = zeros(Float64,Nbin) #expected utility of the bin
    ls   = zeros(Float64,Nbin) #labor supply of the bin
    Uli   = zeros(Float64,Nbin) #utility level of the bin per capita
    Ul   = zeros(Float64,Nbin) #utility level of the bin per capita
    xsiU = zeros(Float64,Nbin)
    xsiUc= zeros(Float64,Nbin)
    inc  = zeros(Float64,Nbin) #income of the agents
    resp = zeros(Float64,Nbin) #fraction of agents credit constrained
    CC = zeros(Float64,Nbin) #share of credit constrained agents
    X = zeros(Float64,Nbin) #object within the utility function (c)
    Transu = spzeros(Nbin,Nbin); #general transition matrix agents
    Size_check = zeros(Float64,Nbin)

    for k=1:Nbin
        i = Ntoi[k] #number of the history
        vE =  convertBasis10E(i, Nidio, N) #this is to recover the vector of histories in which the Size is greater than zero
        yy = mod(vE[1] - 1 ,ns)+1 #productivity type between 1 and ns
        ls[k] = Ls[vE[1]] #labor supply
        inc[k] = ls[k]*w*states[yy] #per capita earnings of the agent. This is the income of the agents
        id0 = 1+ (vE[N]-1)*na  #index of the initial distribution in D_ss
        if0 =  (vE[N])*na #final index of the initial distribution in D_ss
        DD = 0 .*D_ss
        DD[id0:if0] .= D_ss[id0:if0] #this is the initial distribution of agents with the beginning of history 4 for example, considering the case in which we have k = 18 and in this case we have i = 28 and vE= [5 4]
        for j=1:N-1
                Mat_trans = 0 .*Mat
                id0 = 1+ (vE[N+1-j]-1)*na
                if0 =  (vE[N+1-j])*na
                id1 = 1+ (vE[N-j]-1)*na
                if1 =  (vE[N-j])*na
                Mat_trans[id1:if1,id0:if0] = Mat[id1:if1,id0:if0]
                DD = Mat_trans*DD
         end
         Size_check[k] = sum(DD) #this is to check whether the size we obtained above is the same as this one
         CC[k]         = sum(resE.*DD)/Size_check[k] #fraction of CC agents in k
         conso[k]      = sum(DD.*Cpol) #this is the consumption in the bin k
         at[k]         = sum(dGridl_t.*DD) #beginnning of period assets in the bin k
         ae2[k]        = sum(dGridl_t.*(Mat*DD)) #total assets in the bin k
         ae[k]         = at[k]*Aggs_ss.R + inc[k]*Size_check[k] + Aggs_ss.TT*Size_check[k] - conso[k] #this represents the total assets
         uc[k]         = sum(upc.*DD) #marginal utility of the bin
         Euc[k]        = sum(upc.*(Mat*DD)) #following utility
         ECib[k]        = sum(Ci.*(Mat*DD)) #Expected consumption in the bin
         ECi2b[k]       = sum(Ci2.*(Mat*DD)) #Expected consumption in the bin
         EUCib[k]       = sum(UCi.*(Mat*DD)) #Expected utility of consumption in the bin
         Euc[k]        = sum(upc.*(Mat*DD)) #following utility
         X[k]          = conso[k]/Size_check[k] #object inside the utility function
         Ul[k]         = Uvv(X[k]) #utility of the bin
         Uli[k]         = sum(DD.*Uvv(Cpol))/Size_check[k] #average utility of the agents in the bin in leve 
         resp[k]       = sum(resi.*DD)/ Size_check[k] #fracton of CC agents in k
     end

     ECi = ECib./Size_check #Expected consumption of agents in the bin
     EUCi = EUCib./Size_check #Expected utility of consumption of agents in the bin
     ECi2 = ECi2b./Size_check #Expected consumption to power 2 of agents in the bin

     Vyn = ECi2 .- X.^2

     #xsiyn = EUCi./Ul
     xsiyn = Uli./Ul

     #Ul         =Uvv(X)

     Ucccp = (γ).*(γ+1).*(X.^(γ-2)) #this is the third derivative
     Uccp = (-γ).*(X.^(-γ-1)) #this is the second derivative with respect to c
     UcpV = X.^-γ #this is the first derivative with respect to c

     xsiyn2 = 1 .+ ((Vyn)./2).*(Uccp./Ul)

     xsis = 1 .+ ((Vyn)./2).*(Ucccp./UcpV)

     #############Transition Matrix for agents and wealth
     #constructing the transition matrices on the projection
     Transp       = spzeros(Nbin,Nbin); # general transition matrix AGENTS

    for k  = 1:Nbin
        i  = Ntoi[k] #number of the history
        v  = convertBasis10E(i,Nidio,N)
        statei   = v[1]
        vf = v                     #to initialize
        vf[2:end] = v[1:end-1]     #translated vector
        for statej=1:Nidio
            if Trans[statej,statei]>0 #if possible continuation
                vf[1] = statej
                iff = convertBasisE10(vec(vf),Nidio,N)
                j = findall(x->x==iff,Ntoi) #index of the history
                if isempty(j)
                else j=j[1]
                     Transp[j,k] = Trans[statej,statei] # from k to j per capita
                end
            end
        end
    end

    #############Define the number of credit constrainted bins
    DDD = reshape(D_ss,na,ns) #disitribution of end-of-period wealth
    Sh = sum(DDD[1,:])
    XX = [CC Size_check]
    XX= XX[sortperm(XX[:, 1],rev=true), :]
    ca=cumsum(XX[:,2])
    
    nc= findall(x->x<Sh, ca)
    nc = nc[end] #number of credit constrained bins

    #nc= findall(x->x>=Sh, ca)
    #nc = nc[1] #number of credit constrained bins

    Th = sort(CC,rev=true)[nc] #the item nc of CC in decreasing order

    Id = sparse(I,Nbin,Nbin)
    indcc=findall(x->x>=Th,CC)        #index of bin where fraction of CC is above Th. Those will be the bins where thte agents are credit constrained
    indnc=findall(x->x<Th,CC)        #index of bin where fraction of CC is below Th. Those will be the bins where the agents are not credit constrained

    #indcc = argmax(proj.Res)
    PP = copy(Id) #this is a diagonal matrix having 1 on the diagonal only if the history is not credit constrained
    for ii=1:length(indcc)
        PP[indcc[ii],indcc[ii]]= 0
    end
    PPc = Id - PP #this is the definition of Pc in the paper

    Onev = ones(Float64,Nbin)
    Ucv = X.^-γ
    Mat2 = PP*(Id - β*R*(Transp'))*diagm(Ucv) + PPc #to obtain the values of the epsilon in the euler equation
    xsiUc = inv(Mat2)*PPc*Onev #For the credit constrained cases the values of the xsi will be equal to zero
    #xsiUc = inv(Mat2)*PPc*resp #For the credit constrained cases the values of the xsi will be equal to zero
    Res = ((xsiUc.*(X.^-γ)) -  β*R*(Transp')*(xsiUc.*(X.^-γ))) #we should check whether those values are close to zero

    #############Defining types of projection
    ytypeb = zeros(Int64,Nbin)

        for k=1:Nbin
            i = Ntoi[k] #number of the history
            vE =  convertBasis10E(i, Nidio, N)
            ytypeb[k] = mod(vE[1] - 1 ,ns)+1 #productivity type between 1 and ns

    end

    #############Computing the standard deviation of the distribution
    EΓ   = 0
    EΓc   = 0
    EΓ2 = 0

    for i=1:Nbin
            Vi = zeros(Float64,size(Transp)[1])
            Vi[i]=1
            ind = findall(x->x>10^-6,Transp*Vi) #vector of possible successors of history i
            for j=ind
                EΓc = EΓc +Size_check[i]*Transp[j,i]*(at[j]/Size_check[j]-ae[i]/Size_check[i]) #should be 0
                EΓ = EΓ +Size_check[i]*Transp[j,i]*((at[j]/Size_check[j]-ae[i]/Size_check[i]) /(at[j]*Aggs_ss.R/Size_check[j] + inc[j] ) ) #divided by income per capita
                EΓ2 = EΓ2 + Size_check[i]*Transp[j,i]*((at[j]/Size_check[j]-ae[i]/Size_check[i]) /(at[j]*Aggs_ss.R/Size_check[j] + inc[j] ) )^2
            end
    end

        stdΓ = (EΓ2 - EΓ^2)^0.5

    Proj = Projection(
    Aggs_ss,
    N,
    ns,        #number of states
    Nbin,         #number of bins
    R,
    w,
    A, #capital stock
    B,
    states,
    Transv,
    Size_check,
    at,
    ae,
    Ls,
    conso, #conso (not per capita)
    Ucv,
    Transp,
    uc, #average marginal utility in the bin
    Ul, #utility in level of the bin (not ajusted)
    Res,    #xsiU,#level
    xsiUc, #marginal
    Ntoi,
    CC,
    indnc,
    indcc,
    ytypeb,
    xsiyn,
    xsiyn2)

    return Proj,upc,resE, stdΓ, resi,ytype
end
