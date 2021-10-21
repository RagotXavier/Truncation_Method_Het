#########################################################
 # Main file to run the results using the BKM Method
#########################################################

#############Packages
using LinearAlgebra
using Parameters
using IterativeSolvers
using FastGaussQuadrature
using ForwardDiff
using QuantEcon
using Plots
using Arpack
using BenchmarkTools
using JLD2
#using XLSX
using MATLAB
using Plots; pyplot();
using SparseArrays;
using MAT
using Statistics


#############Declaration of variables and basic functions. The structure is based on the EGM code of Alisdair McKay (errors ar ours)
struct AiyagariParametersEGM{T <: Real}
    β::T
    α::T
    δ::T
    γ::T
end


struct AiyagariModelEGM{T <: Real,I <: Integer}
    params::AiyagariParametersEGM{T} #parameters
    aGrid::Array{T,1} #policy grid
    aGridl::Array{T,1}
    na::I #number of grid points in policy function
    dGrid::Array{T,1} #distribution grid
    nd::I #number of grid points in distribution
    states::Array{T,1} #earning states
    ns::I #number of states
    Transv::Array{T,1} #transition for productivity
end


mutable struct AggVarsEGM{S <: Real,T <: Real}
    R::S #interes-rate
    w::T #wage
    TT::T #transfer
end


#############Gini Wealth
function gini2(wagedistarray)
    Swages = cumsum(wagedistarray[:,1].*wagedistarray[:,2])
    Gwages = Swages[1]*wagedistarray[1,2] +
                sum(wagedistarray[2:end,2] .*
                        (Swages[2:end]+Swages[1:end-1]))
    return 1 - Gwages/Swages[end]
end


function wealthD(ww::AbstractArray,NN::Integer)
    ww2 = sortslices(ww, dims=1)
    fracw = cumsum(ww2[:,1].*ww2[:,2])
    fracw =fracw ./fracw[end]
    fracS = cumsum(ww2[:,2])
    plot(fracS,fracw)
    vect = zeros(Float64,NN,2)
    is  = 0
    for i=1:NN
        is = maximum(findall(x->x<=i/NN,fracS))
        vect[i,1] = fracS[is] #- vect[i-1,1]
        vect[i,2] = fracw[is] #- vect[i-1,2]
    end
    vecf = zeros(Float64,NN,2)
    vecf[1,:] .= vect[1,:]
    for i=2:NN
        vecf[i,1] = vect[i,1] - vect[i-1,1]
        vecf[i,2] = vect[i,2] - vect[i-1,2]
    end
    return  vecf
end


##############################
#   MAIN
##############################
include("Aiyagari_solve.jl")
include("Parameters.jl")

K0= 44.98314527071115
AA0,pol0,Aggs0 = AiyagariEGM(K0,ns,Transv,states,amin,amax,curv,β,α,δ,γ,TT,na)
polA_ss,polC_ss,D_ss,A,Aggs_ss,Mat,upc,Ls,Ulevel = steady(K0,pol0,AA0,R,w,TT,1e-4)

Ltot = dot(vstan,states) #total labor
A = (1+tc)*A #total assets corrected for consumption taxes
B = A-A #total debt in the economy
Ctot = D_ss'*polC_ss #total consumption, it has to be equal to Ctot = Aggs_ss.w*Ltot + (Aggs_ss.R-1)*A + TT
Y = A^α*Ltot^(1-α) #total output
G = Y- δ*A - Ctot #total government expenses, it has to be close to zero
Iss = δ*A

Grid = repeat(AA0.aGrid,ns)
D2 = copy(D_ss)
ww = [Grid D2]
ww2 = sortslices(ww, dims=1)
Gini = gini2(ww2)
vW = wealthD(ww2,5)

d_initial_ss= D_ss


#############Computing the evolution of prices
function Pricesv(Ki::T,
        K::Array{T,1},
        Z::Array{T,1},
        params::AiyagariParametersEGM,
        TT::T) where{T <: Real}

    @unpack β,α,δ,γ = params

    Periods = length(K)

    R= zeros(Float64,Periods)
    w= zeros(Float64,Periods)
    L= zeros(Float64,Periods)

    L[1] =dot(vstan,states)
    w[1] = (1- α)*Z[1]*(Ki/L[1])^(α)
    R[1] = α*Z[1]*(Ki/L[1])^(α-1.0) + 1.0 - δ

    for i=2:Periods
        L[i] =dot(vstan,states)
        w[i] = (1- α)*Z[i]*(K[i-1]/L[1])^(α)
        R[i] = α*Z[i]*(K[i-1]/L[i])^(α-1.0) + 1.0 - δ
    end

    return R,w,L
end

Ki=A
Periods = 250 #number of periods for the simulation
K = Ki*ones(Float64,Periods)
Z = ones(Float64,Periods)
eps = 0.01
Z[1]=1+eps
K[1] = Ki
rhoz = 0.95
for i=2:Periods
    Z[i] = 1+rhoz*(Z[i-1]-1)
    K[i] = Ki #*Z[i]
end

R,w,L = Pricesv(Ki,K,Z,AA0.params,TT)

#plot(R)
#plot(w)


##############################
#   TRANSITION AND IMPULSE RESPONSE FUNCTIONS
##############################
include("MIT_functions.jl")

pol_ss = polA_ss #initial policy for assets
polc_ss =polC_ss #initial policy for consumption
dis_ss = D_ss #initial distribution
U_ss = Ulevel #initial policy for utility

K_ss = A
TimePeriods = 200 #time periods for the simulation
OffeqPathTime = 1

#steady state productivity
z_ss = 1.0

#############Create capital path guess and path for shocks
Kguess = K_ss*ones(TimePeriods+1) #add time 0
Zpath  = ones(Float64,TimePeriods+1) #add time 0

epsz = 0.00312 #standard deviation for the shock
rhoz = 0.95

Zpath[1]=1+epsz
σf = 0

for i=2:TimePeriods
    Zpath[i] = 1+rhoz*(Zpath[i-1]-1)
end


#############Finding the equilibrium
eq_apols,eq_cpols,eq_dpols,eq_Utpols,eq_aggKs = equilibrium(Kguess,Zpath,polc_ss,U_ss,pol_ss,dis_ss,AA0,TT,σf)


#############Constructing the IRF for the Aggregate variables
function Aggregates(K::AbstractArray,TT,Z::AbstractArray,params::AiyagariParametersEGM) #K is here K(-1)
#before tax
    @unpack β,α,δ,γ = params
    L =fill(dot(vstan,states),length(K),1)
    w = Z.*(1-α).*(K./L).^(α)
    R = α*Z.*(K).^(α-1.0) .+ 1.0 .- δ
    return R,w,L
end

Rv,wv,Lv = Aggregates(eq_aggKs,TT,Zpath,AA0.params);
Yv = Zpath.*Lv.^(1-α).*eq_aggKs.^α
Iv =  eq_aggKs[2:end] - (1 - AA0.params.δ)*eq_aggKs[1:end-1]
Cv =   Yv[1:end-1]- Iv


#############Determining the consumption in percentiles of the population
ww = [Grid D_ss collect(1:na*ns)]
ww2 = sortslices(ww, dims=1)

fracS2 = cumsum(ww2[:,2]) #determine the cumulative distribution of the population

b = findall(x -> x<=0.10, fracS2) #bottom 10 percent
b = Int.(ww2[:,3][b]) #find the points where the bottom 10 percent are located
Cbsa_ss = D_ss[b]'*polC_ss[b] #total consumption of the bottom 10 percent
Cbs_ss = Cbsa_ss/sum(D_ss[b])

ww3= ww2[sortperm(ww2[:, 1],rev=true), :]
fracS3 = cumsum(ww3[:,2])

m= findall(x -> x <=0.50, fracS3) #top 50 percent
m = Int.(ww3[:,3][m]) #find the points where the top 50 percent are located
Cmsa_ss = D_ss[m]'*polC_ss[m] #total consumption of the top 50 percent
Cms_ss = Cmsa_ss/sum(D_ss[m])

t= findall(x -> x <=0.10, fracS3) #top 10 percent
t = Int.(ww3[:,3][t]) #find the points where the top 10 percent are located
Ctsa_ss = D_ss[t]'*polC_ss[t] #total consumption of the top 10 percent
Cts_ss = Ctsa_ss/sum(D_ss[t])


#############Determining the consumption in percentiles of the population for all the periods
WW = [Grid eq_dpols collect(1:na*ns)]
WW2 = sortslices(WW, dims=1)

fracS2  = Array{Float64, 2}(undef, na*ns, TimePeriods)

Cbsa = Array{Float64, 2}(undef, TimePeriods, 1)
Cbs = Array{Float64, 2}(undef, TimePeriods, 1)
for i = 1: TimePeriods
    fracS2[:,i] = cumsum(WW2[:,i+2])
    global b = findall(x -> x<=0.10, fracS2[:,i])
    global b = Int.(WW2[:,end][b])
    Cbsa[i,1] = eq_dpols[b, i+1]'*eq_cpols[b,i]
    Cbs[i,1] = Cbsa[i,1]/sum(eq_dpols[b, i+1])
end
Cbs_eps = 100*(Cbs .-Cbs_ss)./(Cbs_ss)
plot(Cbs_eps)

WW3= WW2[sortperm(WW2[:, 1],rev=true), :]
fracS3  = Array{Float64, 2}(undef, na*ns, TimePeriods)

Cmsa = Array{Float64, 2}(undef, TimePeriods, 1)
Cms = Array{Float64, 2}(undef, TimePeriods, 1)
for i = 1: TimePeriods
    fracS3[:,i] = cumsum(WW3[:,i+2])
    global m = findall(x -> x<=0.50, fracS3[:,i])
    global m = Int.(WW3[:,end][m])
    Cmsa[i,1] = eq_dpols[m, i+1]'*eq_cpols[m,i]
    Cms[i,1] = Cmsa[i,1]/sum(eq_dpols[m, i+1])
end
Cms_eps = 100*(Cms .-Cms_ss)./(Cms_ss)
plot(Cms_eps)

Ctsa = Array{Float64, 2}(undef, TimePeriods, 1)
Cts = Array{Float64, 2}(undef, TimePeriods, 1)
for i = 1: TimePeriods
    fracS3[:,i] = cumsum(WW3[:,i+2])
    global t = findall(x -> x<=0.10, fracS3[:,i])
    global t = Int.(WW3[:,end][t])
    Ctsa[i,1] = eq_dpols[t, i+1]'*eq_cpols[t,i]
    Cts[i,1] = Ctsa[i,1]/sum(eq_dpols[t, i+1])
end
Cts_eps = 100*(Cts .-Cts_ss)./(Cts_ss)
plot(Cts_eps)


#############Computing total Welfare in the Steady state
CC=(Mat*D_ss).*polC_ss  #sum(CC) = Ctot
ii=findall(x -> x==0, D_ss)
ci = CC./D_ss
Ul =Uvv(ci)
Ul[ii] .=0
Ws =sum(D_ss.*Ul)

Wsv = zeros(Float64,TimePeriods)
for i=1:TimePeriods
    Wsv[i] = sum(eq_dpols[:,i+1].*eq_Utpols[:,i])
end


#############IMPULSE RESPONSE FUNCTION FOR TOTAL CONSUMPTION - TEST TO SEE WHETHER WHETHER IT WAS CORRECT
Csa = Array{Float64, 2}(undef, TimePeriods, 1)
Cs = Array{Float64, 2}(undef, TimePeriods, 1)
for i = 1: TimePeriods
    fracS2[:,i] = cumsum(WW2[:,i+2])
    global b = findall(x -> x<=1.00, fracS2[:,i])
    global b = Int.(WW2[:,end][b])
    Csa[i,1] = eq_dpols[b, i+1]'*eq_cpols[b,i]
    Cs[i,1] = Csa[i,1]/sum(eq_dpols[b, i+1])
end
Cs_eps = 100*(Cs .-Ctot)./(Ctot)


#############IMPULSE RESPONSE FUNCTIONS
Z_eps = 100*(Zpath.-1)
GDP_eps = 100*(Yv.-Y)./Y
C_eps = 100*(Cv.-(Y- δ*A))./(Y- δ*A)
Cbs_eps
Cms_eps
Cts_eps
K_eps = 100*(eq_aggKs.-A)./A
I_eps = 100*(Iv.-Iss)./Iss
#W_eps = 100*(Wsv .- Ws)./Ws
W_eps = 100*(Wsv .- Ws)
R_eps = 100*(Rv .- Aggs_ss.R)
w_eps = 100*(wv .- Aggs_ss.w)./Aggs_ss.w

plot(collect(0:1:TimePeriods),Z_eps,label="Shock")
plot(collect(0:1:TimePeriods),GDP_eps,label="Aggregate output path from Mit shock")
plot(collect(0:1:TimePeriods),K_eps,label="Aggregate capital path from Mit shock")
plot(collect(0:1:TimePeriods-1),C_eps,label="Aggregate consumption path from Mit shock")
plot(collect(0:1:TimePeriods-1),Cbs_eps,label="Consumption path for the bottom from Mit shock")
plot(collect(0:1:TimePeriods-1),Cms_eps,label="Consumption path for the median from Mit shock")
plot(collect(0:1:TimePeriods-1),Cts_eps,label="Consumption path for the top from Mit shock")
plot(collect(0:1:TimePeriods-1),I_eps,label="Aggregate investment path from Mit shock")
plot(collect(0:1:TimePeriods-1),W_eps,label="Aggregate Welfare path from Mit shock")
plot(collect(0:1:TimePeriods),R_eps,label="Interest rate path from Mit shock")
plot(collect(0:1:TimePeriods),w_eps,label="Wage path from Mit shock")


#############TIME SERIES FOR THE RELEVANT VARIABLES DEVIATION FROM THE STEADY STATE
#I construct the relevant multipliers, I have to divide by the size of the shock to normalize
coefZt = (Z_eps.-Z_eps[end])/epsz
coefCt = (C_eps.-C_eps[end])/epsz
coefKt = (K_eps.-K_eps[end])/epsz
coefYt = (GDP_eps.-GDP_eps[end])/epsz
coefIt = (I_eps.-I_eps[end])/epsz
coefWt = (W_eps.-W_eps[end])/epsz
coefRt = (R_eps.-R_eps[end])/epsz
coefwt = (w_eps.-w_eps[end])/epsz


#I read here the file of innovations from MATLAB use in other methods
file = matopen("MATepsv2_e.mat")
uv2=read(file,"uv2")
close(file)

#I compute the path of the variables for the aggregate shock, using the multiplier
SimulPeriods = length(uv2)
zv = zeros(Float64,SimulPeriods)
rhoz = 0.95
zv[1] = uv2[1]
Zt  = zeros(Float64,SimulPeriods-TimePeriods)
GDPt  = zeros(Float64,SimulPeriods-TimePeriods)
Consot    = zeros(Float64,SimulPeriods-TimePeriods)
Capitalt  = zeros(Float64,SimulPeriods-TimePeriods)
Investmentt   = zeros(Float64,SimulPeriods-TimePeriods)
Welfaret  = zeros(Float64,SimulPeriods-TimePeriods)
Rt  = zeros(Float64,SimulPeriods-TimePeriods)
wt  = zeros(Float64,SimulPeriods-TimePeriods)

for i=2:SimulPeriods
    zv[i] = rhoz*zv[i-1]+uv2[i]
end

for i=1:(SimulPeriods-TimePeriods)
    for j=1:TimePeriods
        Zt[i]     +=  coefZt[j]*uv2[TimePeriods-j+i]
        GDPt[i]     +=  coefYt[j]*uv2[TimePeriods-j+i]
        Consot[i]   +=  coefCt[j]*uv2[TimePeriods-j+i]
        Capitalt[i] +=  coefKt[j]*uv2[TimePeriods-j+i]
        Investmentt[i] +=  coefIt[j]*uv2[TimePeriods-j+i]
        Welfaret[i]  +=  coefWt[j]*uv2[TimePeriods-j+i]
        Rt[i]  +=  coefRt[j]*uv2[TimePeriods-j+i]
        wt[i]  +=  coefwt[j]*uv2[TimePeriods-j+i]

    end
end


#############TIME SERIES FOR THE RELEVANT VARIABLES
#I construct the relevant multipliers, I have to divide by the size of the shock to normalize
coefC = (Cv.-Cv[end])/epsz
coefK = (eq_aggKs.-A[end])/epsz
coefY = (Yv.-Yv[end])/epsz
coefL = (Lv.-Lv[end])/epsz
coefI = (Iv.-Iv[end])/epsz
coefW = (Wsv.-Wsv[end])/epsz
coefR = (Rv.-Rv[end])/epsz
coefw = (wv.-wv[end])/epsz


#I read here the file of innovations from MATLAB use in other methods
file = matopen("MATepsv2_e.mat")
uv2=read(file,"uv2")
close(file)

#I compute the path of the variables for the aggregate shock, using the multiplier
SimulPeriods = length(uv2)
zv = zeros(Float64,SimulPeriods)
rhoz = 0.95
zv[1] = uv2[1]
GDP  = zeros(Float64,SimulPeriods-TimePeriods)
Conso    = zeros(Float64,SimulPeriods-TimePeriods)
Capital  = zeros(Float64,SimulPeriods-TimePeriods)
Investment   = zeros(Float64,SimulPeriods-TimePeriods)
Labour  = zeros(Float64,SimulPeriods-TimePeriods)
Welfare  = zeros(Float64,SimulPeriods-TimePeriods)
Rsimul  = zeros(Float64,SimulPeriods-TimePeriods)
wsimul  = zeros(Float64,SimulPeriods-TimePeriods)

for i=2:SimulPeriods
    zv[i] = rhoz*zv[i-1]+uv2[i]
end

for i=1:(SimulPeriods-TimePeriods)
    for j=1:TimePeriods
        GDP[i]     +=  coefY[j]*uv2[TimePeriods-j+i]
        Conso[i]   +=  coefC[j]*uv2[TimePeriods-j+i]
        Capital[i] +=  coefK[j]*uv2[TimePeriods-j+i]
        Investment[i] +=  coefI[j]*uv2[TimePeriods-j+i]
        Labour[i]  +=  coefL[j]*uv2[TimePeriods-j+i]
        Welfare[i]  +=  coefW[j]*uv2[TimePeriods-j+i]
        Rsimul[i]  +=  coefR[j]*uv2[TimePeriods-j+i]
        wsimul[i]  +=  coefw[j]*uv2[TimePeriods-j+i]
    end
end

#Normalized Standard deviation and correlations
@show std(GDP)/mean(Yv)
@show std(Conso)/mean(Cv)
@show std(Capital)/mean(eq_aggKs)
@show std(Investment)/mean(Iv)
@show std(Labour)/mean(Lv)
@show std(Welfare)/mean(Wsv)
@show std(Welfare)

@show cor(GDP[1:end],Conso[1:end])
@show cor(GDP[1:end],Capital[1:end])
@show cor(GDP[1:end],Investment[1:end])

@show cor(GDP[1:end-1],GDP[2:end])
@show cor(Conso[1:end-1],Conso[2:end])
@show cor(Capital[1:end-1],Capital[2:end])
@show cor(Investment[1:end-1],Investment[2:end])



##############################
#   Saving to figure
##############################
file = matopen("tofigBKM.mat", "w")
write(file, "Z_eps", Z_eps)
write(file, "GDP_eps", GDP_eps)
write(file, "K_eps", K_eps)
write(file, "C_eps", C_eps)
write(file, "Cbs_eps", Cbs_eps)
write(file, "Cms_eps", Cms_eps)
write(file, "Cts_eps", Cts_eps)
write(file, "I_eps", I_eps)
write(file, "W_eps", W_eps)
write(file, "R_eps", R_eps)
write(file, "w_eps", w_eps)
write(file, "Zt", Zt)
write(file, "GDPt", GDPt)
write(file, "Consot", Consot)
write(file, "Capitalt", Capitalt)
write(file, "Investmentt", Investmentt)
write(file, "Welfaret", Welfaret)
write(file, "Rsimul", Rsimul)
write(file, "wsimul", wsimul)

close(file)


