#########################################################
 # Below the functions to solve Aiyagari through EGM
#########################################################

function PricesEGM(K,TT,Z,params::AiyagariParametersEGM) #this function set initial prices assuming no taxes (as a initial guess)
    @unpack β,α,δ,γ = params
    L =dot(vstan,states)
    w = (1- α)*Z*(K/L)^(α)
    R = Z*α*(K/L)^(α-1.0) + 1.0 - δ
    return AggVarsEGM(R,w,TT)
end


function AiyagariEGM(
    K::T,
    ns::I,
    Transv::Array{T,1},
    states::Array{T,1},
    amin::T = 1e-9,
    amax::T = 200.0,
    curv::Real = 4,
    β::T = 0.98,
    α::T = 0.4,
    δ::T = 0.02,
    γ::T = 1.5,
    TT::T = 0.0,
    na::I = 201,
    nd::I = 201) where{T <: Real,I <: Integer}

    #############Parameters
    params = AiyagariParametersEGM(β,α,δ,γ)
    AggVars = PricesEGM(K,TT,1.0,params)
    @unpack R,w = AggVars

    #############Policy grid
    function grid_fun(a_min,a_max,na, pexp)
        x = range(a_min,step=0.5,length=na)
        grid = a_min .+ (a_max-a_min)*(x.^pexp/maximum(x.^pexp))
        return grid
    end
    aGrid = grid_fun(amin,amax,na,curv)

    #############Distribution grid
    #dGrid = collect(range(aGrid[1],stop = aGrid[end],length = nd))
    #aGrid = collect(range(amin,stop = amax,length = na))

    dGrid=aGrid

    guess = repeat(10 .+ aGrid,ns)
    aGridl_t = repeat(aGrid,ns)
    return AiyagariModelEGM(params,aGrid,aGridl_t,na,dGrid,nd,states,ns,Transv),guess,AggVars
end


function interpEGM(
    pol::AbstractArray,
    grid::AbstractArray,
    x::T,
    na::Integer) where{T <: Real}
    np = searchsortedlast(pol,x)

    #if you give x (current asset) (out of the policy pol), gives ap (saving decision) out of the grid (grid)
    (np > 0 && np < na) ? np = np : #adjust indices if assets fall out of bounds
        (np == na) ? np = na-1 :
            np = 1

    ap_l,ap_h = pol[np],pol[np+1]
    a_l,a_h = grid[np], grid[np+1]
    ap = a_l + (a_h-a_l)/(ap_h-ap_l)*(x-ap_l)

    above =  ap > 0.0
    return above*ap,np #return either 0 or ap
end


function get_cEGM(
    pol::AbstractArray,
    Aggs::AggVarsEGM,
    CurrentAssets::AbstractArray,
    AiyagariModel::AiyagariModelEGM,
    cpol::AbstractArray)

    @unpack params,aGrid,na,ns,states = AiyagariModel
    X = copy(states)
    X = ((Aggs.w .*states))

    pol = reshape(pol,na,ns)
    for si = 1:ns
        for ai = 1:na
            asi = (si - 1)*na + ai
             #cpol[asi] = Aggs.R*CurrentAssets[asi] + Aggs.w*states[si] - interpEGM(pol[:,si],aGrid,CurrentAssets[asi],na)[1]
             cpol[asi] = Aggs.R*CurrentAssets[asi] + Aggs.TT + X[si] - interpEGM(pol[:,si],aGrid,CurrentAssets[asi],na)[1]
        end
    end
    return cpol
end


function get_ZEGM(
    pol::AbstractArray,
    Aggs::AggVarsEGM,
    CurrentAssets::AbstractArray,
    AiyagariModel::AiyagariModelEGM,
    cpol::AbstractArray)

    @unpack params,aGrid,na,ns,states = AiyagariModel
    X = copy(states)
    X = ((Aggs.w .*states))
    Z = copy(cpol)
    pol = reshape(pol,na,ns)
    for si = 1:ns
        for ai = 1:na
            asi = (si - 1)*na + ai
            #cpol[asi] = Aggs.R*CurrentAssets[asi] + Aggs.w*states[si] - interpEGM(pol[:,si],aGrid,CurrentAssets[asi],na)[1]
            Z[asi] = Aggs.R*CurrentAssets[asi] + Aggs.TT + X[si] - interpEGM(pol[:,si],aGrid,CurrentAssets[asi],na)[1] #compact form
        end
    end
    return Z
end


function EulerBackEGM(
    pol::AbstractArray,
    Aggs::AggVarsEGM,
    Aggs_P::AggVarsEGM,
    AiyagariModel::AiyagariModelEGM,
    cpol::AbstractArray,
    apol::AbstractArray)

    @unpack params,na,nd,ns,aGridl,Transv,states = AiyagariModel
    @unpack γ,β= params

    Trans = reshape(Transv,ns,ns)

    R_P,w_P = Aggs_P.R,Aggs_P.w
    R,w,TT = Aggs.R,Aggs.w,Aggs.TT

    #cp = get_cEGM(pol,Aggs_P,aGridl,AiyagariModel,cpol) #vector ns*na
    Zp  = get_ZEGM(pol,Aggs_P,aGridl,AiyagariModel,cpol) #vector ns*na, future prices

    upcp = uPrime(Zp,γ)
    Eupcp = copy(cpol)
    #Eupcp_sp = 0.0

    X = ((Aggs.w .*states))

    Ls = dot(vstan,states) #labor supply

    for ai = 1:na
        for si = 1:ns
            asi = (si-1)*na + ai
            Eupcp_sp = 0.0
            for spi = 1:ns
                aspi = (spi-1)*na + ai
                #Eupcp_sp += Trans[spi,si]*upcp[aspi]*R_P
                Eupcp_sp += Trans[spi,si]*upcp[aspi]
            end
            Eupcp[asi] = Eupcp_sp
        end
    end

    upc = β*Eupcp*R_P
    c = copy(upc)
    Z = uPrimeInv(upc,γ)

    for ai = 1:na
        for si = 1:ns
            asi = (si-1)*na+ai
            c[asi]    = Z[asi]
            apol[asi] = (aGridl[asi] + c[asi] - TT- X[si])/R
            #apol[asi] = (aGridl[asi] + c[asi] - w*states[si])/R
        end
    end
    UU = Uvv.(Z)
    return apol,c,upc,Ls,UU
end


function SolveEGM(
    pol::AbstractArray,
    Aggs::AggVarsEGM,
    AiyagariModel::AiyagariModelEGM,
    cpol::AbstractArray,
    apol::AbstractArray,tol = 1e-10)
    @unpack ns,na = AiyagariModel

    for i = 1:10000
        a = EulerBackEGM(pol,Aggs,Aggs,AiyagariModel,cpol,apol)[1]
        if (i-1) % 50 == 0
            test = abs.(a - pol)/(abs.(a) + abs.(pol))
            #println("iteration: ",i," ",maximum(test))
            if maximum(test) < tol
                println("Solved in ",i," ","iterations")
                break
            end
        end
        pol = copy(a)
    end
    return pol
end


#########################################################
 # Below the functions  for the stationary distribution
#########################################################

function MakeTransMatEGM(pol,AiyagariModel,tmat)
    @unpack ns,na,aGrid,Transv = AiyagariModel
    pol = reshape(pol,na,ns)

    Trans = reshape(Transv,ns,ns)
    tmat = copy(tmat)
    tmat  .=0.0

    for a_i = 1:na
        for j = 1:ns
            x,i = interpEGM(pol[:,j],aGrid,aGrid[a_i],na)
            p = (aGrid[a_i] - pol[i,j])/(pol[i+1,j] - pol[i,j])
            p = min(max(p,0.0),1.0)
            sj = (j-1)*na
            for k = 1:ns
                sk = (k-1)*na
                tmat[sk+i+1,sj+a_i] = p * Trans[k,j]
                tmat[sk+i,sj+a_i] = (1.0-p) * Trans[k,j]
            end
        end
    end
    return tmat
end


function StationaryDistributionEGM(Mat,AiyagariModel)
    @unpack params,ns,nd,na = AiyagariModel
    x = zeros(Float64,(na*ns,1))
    λ, x = powm!(Mat, ones(ns*na)/(ns*na), maxiter = 100000,tol = 1e-8)
    x = x/sum(x)
    return x
end


function steady(
    K0::T,
    initialpol::AbstractArray,
    AiyagariModel::AiyagariModelEGM,
    R::T,
    w::T,
    TT::T =0.0,
    tol = 1e-10,maxn = 50)where{T <: Real}
    @unpack params,aGrid,na,dGrid,nd,ns,na = AiyagariModel

    tmat = zeros(eltype(initialpol),(na*ns,na*ns))
    cmat = zeros(eltype(initialpol),na*ns)
    pol = initialpol

    relax  = 5e-2
    uK,lK = 2*K0, 0.0

    pol = initialpol

    for kit = 1:maxn
        Aggs = PricesEGM(K0,TT,1.0,params) #steady state
        cmat .= 0.0
        pol = SolveEGM(pol,Aggs,AiyagariModel,cmat,cmat)

        #############Stationary transition
        tmat .= 0.0
        trans = MakeTransMatEGM(pol,AiyagariModel,tmat)
        D = StationaryDistributionEGM(trans,AiyagariModel)

        #############Aggregate savings
        dGridl_t = repeat(dGrid,ns)
        EA = dot(dGridl_t,D)

        if (EA > K0) #too little lending -> low r -> too much borrowing
            uK = min(EA,uK)
            lK = max(K0,lK)
            K0 = 1.0/2.0*(lK + uK)
        else #too much lending -> high r -> too little borrowing
            uK = min(K0,uK)
            lK = max(EA,lK)
            K0 = 1.0/2.0*(lK + uK)
        end

        println("It. ",kit," Int. rate: ",Aggs.R," Bond S.: ",EA, " Bond D.: ",K0)

        if abs(EA - K0) < tol
            println("Markets clear!")
            cmat .= 0.0
            polA,polC,upc, Ls, Ulevel = EulerBackEGM(pol,Aggs,Aggs,AiyagariModel,cmat,cmat)
            return polA,polC,D,EA,Aggs,trans,upc, Ls, Ulevel
            break
        end
    end
    return println("Markets did not clear")
end


#########################################################
 # Below the function to calculate the residuals in Euler
#########################################################

function Euler_err(
    N::Integer, #length of the truncation
    AA0::AiyagariModelEGM,
    polA_ss::Array{T,1},
    polC_ss::Array{T,1},
    upc::Array{T,1},
    D_ss::Array{T,1},
    Aggs_ss::AggVarsEGM{T,T},
    Mat::Array{T,2},
    A::T,
    vstan::Array{T,1},
    Ls::Array{T,1},
    Ulevel::Array{T,1},
    states::Array{T,1}) where {T <: Real}

    @unpack params,aGrid,na,dGrid,nd,ns,Transv,states = AA0
    @unpack γ,β = params

    R = Aggs_ss.R
    w = Aggs_ss.w
    TT = Aggs_ss.TT
    Ltot = dot(vstan,states)

    Va = repeat(aGrid,ns)'*Mat #policy rule in a a'(a)
    pol = reshape(Va,na,ns)

    γ =  AA0.params.γ
    β =  AA0.params.β

    Trans   = reshape(AA0.Transv,ns,ns)

    #############Computing CC Bins
    Vmat = Mat[1:na,:]
    aGrid = AA0.aGrid
    aGrid[1] = 0
    Va = repeat(aGrid,ns)'*Mat #policy rule in a a'(a)
    pol = reshape(Va,na,ns)
    nd = zeros(Int,na*ns)
    Wind = zeros(Float64,na*ns)
    for i = 1:length(Va)
        indcc=findall(x->x<=Va[i],aGrid)
        Vind[i] = maximum(indcc)
        Wind[i] = 1-(Va[i] - aGrid[Vind[i]])/(aGrid[Vind[i]+1]-aGrid[Vind[i]])
    end

    γ =  AA0.params.γ
    β =  AA0.params.β
    R = Aggs_ss.R
    Trans   = reshape(AA0.Transv,ns,ns)

    Ls  = fill(dot(vstan,states),ns) #vector of labor supplies
    Lv = kron(Ls[1:ns],ones(na,1))

    ytype = repeat(kron(collect(StepRange(1, Int8(1), ns)),ones(Int,na))) #types of agents
    Cpol =  - Va' + Aggs_ss.R*repeat(aGrid,ns) + Aggs_ss.w*states[ytype].*Lv + Aggs_ss.TT*ones(Float64,na*ns)

    Xe = (Aggs_ss.w .*states)
    Xc = 0*Xe
    Xc = kron(Xc,ones(na,1))

    resi = β*R*(Mat')*upc
    cimp = (resi).^(-1/γ) .+Xc

    #############Construct the vector of 0 for constrained bin and 1 otherwise, using the transition matrix
    gridv = repeat([1;zeros(na-1)],ns) #grid of 1 only if minimum wealt (ie. CC)
    matt = kron(gridv,ones(1,na*ns))
    CCind = findall(x->x<=amin,(gridv.*Mat*ones(na*ns)))
    Vc = 1.0 .-sum(matt.*Mat,dims=1)'

    EulerRes = Vc .*(Cpol./cimp .- 1.) #not population weighted
    return EulerRes
end
