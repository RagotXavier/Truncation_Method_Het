#########################################################
 # Optim Function to find the Pareto Weights such that the weights are positive
#########################################################
using JuMP
using Ipopt

function Optim(eco::Projection,params::AiyagariParametersEGM)
    amin   = minimum(eco.aep./eco.Sp) #minimum value of asset per capita
    ns     = eco.ns #number of idiosyncratic states
    Nbin   = eco.Nbin #number of bins
    w      = eco.w #wage
    R      = eco.R #interest rate
    states = eco.states #idiosyncratic states
    ytype  = eco.ytype #individual types
    β      = params.β #discount factor
    γ      = params.γ #curvature of the utility function
    α      = params.α #capital-share
    δ      = params.δ #depreciation rate
    Ucp    = vec(eco.Ucp) #vector of the utility of the bin (not per capita)
    Sp     = vec(eco.Sp) #size of each bin
    Ap     = Ucp./Sp #utility per capita
    A      = eco.A #total capital in the economy
    tauc   = A/sum(eco.aep)-1 #tax
    Ly     = eco.lb[ytype].*states[ytype] #effective labor supply
    Lv     = eco.lb[ytype] #labor supply
    TT     = eco.Agg.TT #transfers
    PI     = eco.Matab #transition matrix for the bins

    Xp     = vec(eco.cp./eco.Sp)  #argument utility function
    θ      = 0
    uSecond(X) = (-γ).*(X.^(-γ-1))

    Uccp = (-γ).*(Xp.^(-γ-1)) #this is the second derivative with respect to c
    UcpV = Xp.^-γ #this is the first derivative with respect to c

    ####################
    #Constructing H1 and H2 and Finding the Pareto weights ω
    ####################
    Sp[findall(x->x<=10^-20,Sp)] .=10^-10 #this is the vector of steady state histories. This is the N-tot vector of steady-state history by stacking history-sizes for all histories

    Diagp = diagm(Sp) #this is the diagram of S DS
    PIt = transpose(PI) #PI is the transition matrix across histories such that S = ΠS
    iDiagp = diagm(1.0 ./Sp) #this is the value of D^(-1)S

    PIS = Diagp*transpose(PI)*iDiagp # per capita, whereas PIup is not. This is the Π^S

    Id = sparse(I,Nbin,Nbin) #this is the I (Ntot x Ntot) identity matrix
    Bp  = diagm(eco.xsip .*Uccp)*(R*PI - Id)  #this is the value of matrix B

    P = sparse(I,eco.Nbin,eco.Nbin); #P is the diagonal matrix having 1 on the diagonal at y^N if and only if the history y^N is not credit constrained and O otherwise
    for i=eco.indcc
        P[i,i]=0
    end

    Pc = Id - P

    Ltot = sum(Sp.*Ly)
    K = Ltot*(( (1/β - 1 + δ))/α ) ^(-1/(1-α)) #capital
    FL = (1-α)*(K/Ltot)^α #marginal productivity of labor
    FK = 1/β - 1 #marginal productivity of capital

    M  = Id - β*R*PIS +β*((R-1/β)/(Ltot*(1+0*(w-FL)/w)))*Sp*transpose(collect(Ly)) #M Matrix
    J  = -inv(Pc+P*M*Bp)*P*M #J Matrix
    NN  = Id - Bp*inv(Pc+P*M*Bp)*P*M #N Matrix

    abc = vec(eco.abp ./ eco.Sp) #beginning-of-period per capita
    abc[findall(x->x>10^30,abc)] .=0   #if size is very small

    Q = (1/(Ltot*(1+0*(w-FL)/w)))*transpose(Ly)*NN #Q matrix

    H1t= ones(1,Nbin)*(diagm(abc)*NN + diagm(UcpV.*eco.xsip)*PI*J)-A*Q #H1 tilde
    H2t = ones(1,Nbin)*NN - Q #H2 tilde

    H1 = H1t*diagm(Sp.*eco.xsip.*UcpV) #H1
    H2 = H2t*diagm(Sp.*eco.xsip.*UcpV) #H2

    ####################
    #Finding the weights
    ####################
    ω=ones(Nbin,1)
    transpose(H1)
    transpose(H2)
    eco.Sp

    C=ones(3,Nbin)
    C[1,:] = H1
    C[2,:] = H2
    C[3,:] = eco.Sp

    d=ones(3,1)
    d[1,1] = 0
    d[2,1] = 0
    d[3,1] =1

    A= 1* Matrix(I, Nbin, Nbin)
    y=ones(Nbin,1)

    f(x...) = sum((x[i]-1)^2 for i in 1:length(x))

    model = Model();
    vector_model = Model(Ipopt.Optimizer)


    @variable(vector_model, x[1:eco.Nbin] >= 0)

    f(x...) = sum((x[i]-1)^2 -1 for i in 1:length(x))
    register(vector_model, :f, eco.Nbin, f; autodiff = true)

    @constraint(vector_model, C * x .== d)
    @NLobjective(vector_model, Min, f(x...) )

    optimize!(vector_model)

    objective_value(vector_model)

    length(x)

    ω =ones(eco.Nbin,1)

    for i = 1:length(x)
        ω[i] =value(x[i])
    end

    ####################
    #Check
    ####################
    Uv  = diagm(Sp.*eco.xsip.*UcpV)*vec(ω)
    Uvc  = eco.xsip.*UcpV.*vec(ω)
    Sλ  = J*Uv
    Sψ  = NN*Uv
    μ   = Q*Uv

    ψh = μ*ones(Nbin,1)- Sψ./Sp
    λ = Sλ./Sp
    λt  = iDiagp*PI*Diagp*λ

    CC = zeros(Float64,6)
    CC[1] = maximum(ψh - (μ*ones(Nbin,1)- Uvc  + diagm(eco.xsip .*Uccp)*(λ - R*λt)))
    CC[2] = maximum(P*ψh - β*R*P*PIt*ψh)
    CC[3] = maximum(Pc*λ)
    CC[4] = maximum(sum(Ly.*Sp.*ψh) - μ*Ltot*0*(FL-w)/w)
    CC[5] = sum( ψh.*vec(eco.abp)) - sum(Sp.*λt.*eco.xsip.*UcpV)
    CC[6] = sum(Sp.*ψh)

    Pl = Plan(vec(λ),vec(λt),vec(ψh),vec(ω),vec(Uvc),vec(CC),μ)
    return Pl
end
