#########################################################
 # MIT Functions
#########################################################

function PricesMIT(K,TT,Z,params::AiyagariParametersEGM) #K is K(-1)
      @unpack β,α,δ,γ = params
      L =dot(vstan,states)
      w = (1- α)*Z*(K/L)^(α)
      R = Z*α*(K/L)^(α-1.0) + 1.0 - δ
      return AggVarsEGM(R,w,TT)
end


function EulerBackMIT(
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
            apol[asi] = (aGridl[asi] + c[asi] + TT- X[si])/R
            #apol[asi] = (aGridl[asi] + c[asi] - w*states[si])/R
        end
    end
    UU = Uvv.(Z)
    return apol,c,upc,Ls,UU
end


function MakeTransMatMIT(pol,AiyagariModel,tmat)
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


function UpdateAggs(
    initialpolc::AbstractArray,
    initialpolW::AbstractArray,
    initialpol::AbstractArray,
    initialdis::AbstractArray,
    Kguess::AbstractArray,
    Zpath::AbstractArray,
    AiyagariModel::AiyagariModelEGM,
    TT::Real,
    σf::Real,
    tol = 1e-10,maxn = 50)

    @unpack params,aGridl,na,nd,ns = AiyagariModel

    TimePeriods = length(Kguess)
    tmat = zeros(eltype(initialpol),(na*ns,na*ns))
    cmat = zeros(eltype(initialpol),na*ns)
    cmat2 = zeros(eltype(initialpol),na*ns)

    apols = zeros(eltype(initialpol),na*ns,TimePeriods) #policy function for assets - transition

    cpols = zeros(eltype(initialpolc),na*ns,TimePeriods-1) #policy function for consumption - transition

    Utpols = zeros(eltype(initialpolW),na*ns,TimePeriods-1) #utility function - transition

    devol = zeros(eltype(initialpol),na*ns,TimePeriods)
    aggK  = zeros(TimePeriods)

    ##Find policies back through time
    pol = initialpol
    apols[:,TimePeriods] = pol

    polc = initialpolc
    cpols[:,TimePeriods-1] = polc

    Ut = initialpolW
    Utpols[:,TimePeriods-1] = Ut

    for i = TimePeriods:-1:2
        #println("time: ",i-1)
        K,Z_P = Kguess[i],Zpath[i] # K[i] is the productive capital in period i, it is K(t-1)
        K_m,Z = Kguess[i-1],Zpath[i-1]
        Pr_P  = PricesMIT(K,TT,Z_P,params) #tomorrow
        Pr    = PricesMIT(K_m,TT,Z,params) #today
        cmat  = zeros(eltype(initialpol),na*ns)
        cmat2 = zeros(eltype(initialpol),na*ns)

        polc = EulerBackMIT(pol,Pr,Pr_P,AiyagariModel,cmat,cmat2)[2]

        Ut = EulerBackMIT(pol,Pr,Pr_P,AiyagariModel,cmat,cmat2)[5]

        pol   = EulerBackMIT(pol,Pr,Pr_P,AiyagariModel,cmat,cmat2)[1]

        apols[:,i-1] = pol #apols[:,i] is thus between i and i+1

        cpols[:, i-1] = polc #cpols[:,i] is thus between i and i+1
        Utpols[:, i-1] = Ut #Utpols[:,i] is thus between i and i+1
    end

    ##Find Evolution of distribution
    dis = initialdis
    trans = MakeTransMatMIT(initialpol,AiyagariModel,tmat)
    #dis = trans*dis
    aggK[1] = dot(aGridl,dis)
    devol[:,1] = dis
    for i = 1:TimePeriods-1  #constructing the distribution at period i+1
        #println("time: ",i)
        pol = apols[:,i]
        polc = cpols[:,i]
        Ut = Utpols[:,i]
        tmat .= 0.0
        trans = MakeTransMatMIT(pol,AiyagariModel,tmat)
        dis = trans*dis
        aggK[i+1] = dot(aGridl,dis) #it is K(t)
        devol[:,i+1] = dis
    end
    return apols,cpols,devol, Utpols, aggK
end


function equilibrium(Kguess,Zpath,initialpolc,initialpolW,initialpol,initialdis,AiyagariModel,TTprime, σf,tol = 1e-10,maxn = 50)

    polc_ss,U_ss,pol_ss,dis_ss, TT =initialpolc,initialpolW, initialpol,initialdis, TTprime
    aggK = copy(Kguess)
    apols,cpols, devol,Utpols,aggK = UpdateAggs(polc_ss,U_ss ,pol_ss,dis_ss,Kguess,Zpath,AiyagariModel,TT,σf)

    for i = 1:100
        #@show Kguess
        apols,cpols,devol,Utpols,aggK = UpdateAggs(polc_ss,U_ss,pol_ss,dis_ss,Kguess,Zpath,AiyagariModel,TT,σf)
        @show dif = maximum(abs.(aggK - Kguess))
        if dif < 1e-6
            return apols,cpols,devol,Utpols,aggK
        end
        Kguess = 0.2*aggK + 0.8*Kguess
    end
    #return print("did not converge")
    return apols,cpols,devol,Utpols,aggK
end
