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
using XLSX
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
    R::S #interest-rate
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


using FileIO
data = load("initial_dis_low.jld2")
d_initial = data["D_ss"]


using MAT
fileIn = matopen("tofigtruncation.mat")
dset = read(fileIn)
coef = dset["coef"]
TT = dset["TT_ss"]


##############################
#   MAIN
##############################
include("Aiyagari_solve.jl")
include("Parameters.jl")

tau = collect(0.050:0.0005:0.090)
θ = 0.23620745544082986

Y0 = copy(Y0)
TT = tau*Y0

polA_ss = zeros(Float64,na*ns,length(TT))
polC_ss = zeros(Float64,na*ns, length(TT))
D_ss = zeros(Float64,na*ns, length(TT))
A = zeros(Float64,length(TT))
A_t = zeros(Float64,length(TT))
Prices = zeros(Float64, 3, length(TT) )
Mat = zeros(Float64,na*ns, na*ns, length(TT))
upc = zeros(Float64,na*ns, length(TT))
Ls = zeros(Float64,length(TT))
Ulevel = zeros(Float64,na*ns, length(TT))

Ltot = zeros(Float64,length(TT))
B = zeros(Float64,length(TT))
Ctot = zeros(Float64,length(TT))
Ctot_t = zeros(Float64,length(TT))
Y = zeros(Float64,length(TT))
G = zeros(Float64,length(TT))
Iss = zeros(Float64,length(TT))

for i = 1:length(TT)

    global AA0,pol0,Aggs0 = AiyagariEGM(K0,ns,Transv,states,amin,amax,curv,β,α,δ,γ,TT[i],na)
    polA_ss[:,i],polC_ss[:,i],D_ss[:,i],A[i],Aggs_ss,Mat[:,:,i],upc[:,i],Ls[i],Ulevel[:,i] = steady(K0,pol0,AA0,R,w,TT[i],1e-4)

    Prices[1,i]  = Aggs_ss.w
    Prices[2,i]  = Aggs_ss.R
    Prices[3,i] = Aggs_ss.TT

    Ltot[i] = dot(vstan,states) #total labor
    B[i] = A[i]-A[i] #total debt
    Ctot[i] = D_ss[:,i]'*polC_ss[:,i] #total consumption
    Ctot_t[i] = Prices[1,i]*Ltot[i] + (Prices[2,i]-1)*A[i] - TT[i]
    Y[i] = A[i]^α*Ltot[i]^(1-α) #total output
    G[i] = Y[i]- δ*A[i] - Ctot[i] #total government expenses, it has to be close to TT
    Iss[i] = δ*A[i]
    A_t[i]= ((Prices[2,i]-1 + δ)/(α*Ltot[i]^(1-α)))^(1 /(α-1))
end

d_initial_ss = d_initial


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

Periods = 250 #number of periods for the simulation

Ki = zeros(Float64, length(TT))
Z = zeros(Float64,Periods,length(TT))
R = zeros(Float64,Periods,length(TT))
w = zeros(Float64,Periods,length(TT))
L = zeros(Float64,Periods,length(TT))


for i = 1:length(TT)
    Ki[i]=A[i]
    K = Ki[i]*ones(Float64,Periods, length(TT))
    Z = ones(Float64,Periods,length(TT))
    #eps = 0.01
    eps = 0.0
    Z[1,i]=1+eps
    K[1,i] = Ki[i]
    #rhoz = 0.95
    rhoz = 0
    for j=2:Periods
        Z[j,i] = 1+rhoz*(Z[j-1, i]-1)
        K[j,i] = Ki[i] #*Z[i]
    end
    R[:,i],w[:,i],L[:,i] = Pricesv(Ki[i],K[:,i],Z[:,i],AA0.params,TT[i])
end


##############################
#   TRANSITION AND IMPULSE RESPONSE FUNCTIONS
##############################
include("MIT_functions.jl")

K_ss = A
TimePeriods = 200 #time periods for the simulation
OffeqPathTime = 1

z_ss = 1.0 #steady state productivity

eq_apols = zeros(Float64,na*ns,TimePeriods+1,length(TT))
eq_cpols = zeros(Float64,na*ns,TimePeriods,length(TT))
eq_dpols = zeros(Float64,na*ns,TimePeriods+1,length(TT))
eq_Utpols = zeros(Float64,na*ns,TimePeriods,length(TT))
eq_aggKs = zeros(Float64,TimePeriods+1,length(TT))
Kguess = zeros(Float64,TimePeriods+1,length(TT))


#############Create capital path guess and path for shocks
for i = 1:length(TT)
    Kguess[:,i] = K_ss[i]*ones(TimePeriods+1) #add time 0
end

Zpath  = ones(Float64,TimePeriods+1) #add time 0
epsz = 0.00 #standard deviation for the shock
rhoz = 0.0

Zpath[1]=1+epsz
σf = 0


for i=2:TimePeriods
    Zpath[i] = 1+rhoz*(Zpath[i-1]-1)
end


for i = 1: length(TT) #creating one equilibrium for each value of TT
    Kguess_new = copy(Kguess[:,i])
    polC_ss_new = copy(polC_ss[:,i])
    Ulevel_new = copy(Ulevel[:,i])
    polA_ss_new = copy(polA_ss[:,i])
    D_ss_new = copy(d_initial_ss)
    TT_new = copy(TT[i])

    #############Finding the equilibrium
    eq_apols[:,:,i],eq_cpols[:,:,i],eq_dpols[:,:,i],eq_Utpols[:,:,i],eq_aggKs[:,i] = equilibrium(Kguess_new,Zpath,polC_ss_new,Ulevel_new,polA_ss_new,D_ss_new,AA0,TT_new,σf)
end


#############Constructing the IRF for the Aggregate variables
function Aggregates(K::AbstractArray,TT,Z::AbstractArray,params::AiyagariParametersEGM) #K is here K(-1)
#before tax
    @unpack β,α,δ,γ = params
    L =fill(dot(vstan,states),length(K),1)
    w = Z.*(1-α).*(K./L).^(α)
    R = α*Z.*(K).^(α-1.0) .+ 1.0 .- δ
    return R,w,L
end


Rv = zeros(Float64,TimePeriods,length(TT))
wv = zeros(Float64,TimePeriods,length(TT))
Lv = zeros(Float64,TimePeriods,length(TT))
Yv = zeros(Float64,TimePeriods,length(TT))
Iv = zeros(Float64,TimePeriods,length(TT))
Tv = zeros(Float64,TimePeriods,length(TT))
Cv = zeros(Float64,TimePeriods,length(TT))
Cv2 = zeros(Float64,TimePeriods,length(TT))
CC = zeros(Float64,na*ns,length(TT))
ci = zeros(Float64,na*ns,length(TT))
Ul = zeros(Float64,na*ns,length(TT))
Ws = zeros(Float64,length(TT))
Wsv = zeros(Float64,TimePeriods,length(TT))
Z_eps= zeros(Float64,TimePeriods,length(TT))
GDP_eps= zeros(Float64,TimePeriods,length(TT))
C_eps= zeros(Float64,TimePeriods,length(TT))
K_eps= zeros(Float64,TimePeriods,length(TT))
I_eps= zeros(Float64,TimePeriods,length(TT))
W_eps= zeros(Float64,TimePeriods,length(TT))
R_eps= zeros(Float64,TimePeriods,length(TT))
w_eps= zeros(Float64,TimePeriods,length(TT))
T_eps= zeros(Float64,TimePeriods,length(TT))


for i = 1:length(TT)
    Rv[:,i],wv[:,i],Lv[:,i] = Aggregates(eq_aggKs[2:end,i],TT[i],Zpath[1:end-1],AA0.params);
    Yv[:,i] = Zpath[1:end-1].*Lv[:,i].^(1-α).*eq_aggKs[1:end-1,i].^α
    Iv[:,i] =  eq_aggKs[2:end,i] - (1 - AA0.params.δ)*eq_aggKs[1:end-1,i]
    Tv[:,i] =  coef[1].+ coef[2]*Yv[:,i].+coef[3].*eq_aggKs[2:end,i]
    Cv[:,i] =  Yv[:,i]- Iv[:,i] - Tv[:,i]
    Cv2[:,i] = wv[:,i].*Lv[:,i] + (Rv[:,i]).*eq_aggKs[1:end-1,i] - eq_aggKs[2:end,i]- Tv[:,i]

    #############Computing total Welfare in the Steady state
    CC[:,i]=(Mat[:,:,i]*D_ss[:,i]).*polC_ss[:,i]  #sum(CC) = Ctot
    global ii=findall(x -> x==0, D_ss[:,i])
    ci[:,i] = CC[:,i]./D_ss[:,i]
    Ul[:,i] =Uvv(ci[:,i])
    Ul[ii,i] .=0
    Ws[i] =sum(D_ss[:,i].*(Ul[:,i])) .+ TT[i].^(θ)

    #############Computing the Welfare in the Transition
    for j=1:TimePeriods
        #Wsv[j, i] = sum(eq_dpols[:,j+1, i].*(eq_Utpols[:,j,i])).+ TT[i].^(θ)
        Wsv[j, i] = sum(eq_dpols[:,j+1, i].*(eq_Utpols[:,j,i])).+ TT[i].^(θ)
    end

    #############IMPULSE RESPONSE FUNCTIONS
    Z_eps[:,i] = 100*(Zpath[1:end-1].-1)
    GDP_eps[:,i] = 100*(Yv[:,i].-Y[i])./Y[i]
    C_eps[:,i] = 100*(Cv[:,i].- Ctot_t[i])/ Ctot_t[i]
    K_eps[:,i] = 100*(eq_aggKs[2:end,i].- A_t[i])./A_t[i]
    I_eps[:,i] = 100*(Iv[:,i].-Iss[i])./Iss[i]
    W_eps[:,i] = 100*(Wsv[:,i] .- Ws[i])./Ws[i]
    R_eps[:,i] = 100*(Rv[:,i] .- Prices[2,i])
    w_eps[:,i] = 100*(wv[:,i] .- Prices[1,i])./Prices[1,i]
    T_eps[:,i] = 100*(Tv[:,i].-(Tv[end,i]))./(Tv[end,i])
end


using FileIO

data = load("truncation.jld2")

dif_truncation = data["dif_truncation"]
Sψ = data["Sψ"]
VT = data["VT"]
VTp = data["VTp"]
Rt = data["Rt"]
wt = data["wt"]
Sλ = data["Sλ"]
tauopt_bew = data["tauopt_bew"]

##############################
#   COMPUTING THE WELFARE IN THE TRANSITION
##############################

β_vector = zeros(Float64,length(Wsv[:,1]))
for i = 1:length(Wsv[:,1])
    β_vector[i] = β.^(i-1)
end

Welfare = transpose(β_vector'*Wsv) .+ (β^(length(Wsv[:,1]))./(1 - β)).*Wsv[end,:]
Welfare_max, ind = findmax(Welfare)


##############################
#   COMPUTING THE WELFARE IN THE STEADY STATE
##############################

Welfare_steady = Ws./(1-β)
Welfare_steady_max, ind_steady = findmax(Welfare_steady)

#############Difference between marginal utility for consumption and goverment
UcpV = ci[:,ind].^-γ
UcpV[isnan.(UcpV)] .=0.0
Ucind= sum(D_ss[:,ind].*(UcpV))
Vgind= ((θ).*TT[ind]^(θ-1))
difind = Ucind - Vgind

#############Difference between marginal utility for consumption and goverment
UcpV = ci[:,ind_steady].^-γ
UcpV[isnan.(UcpV)] .=0.0
Ucind_steady=sum(D_ss[:,ind_steady].*(UcpV))
Vgind_steady= ((θ).*TT[ind_steady]^(θ-1))
difind_steady = Ucind_steady - Vgind_steady

UcpV = ci.^-γ
UcpV[isnan.(UcpV)] .=0.0
Uc = zeros(Float64, length(TT))
for i = 1:length(TT)
    Uc[i] = D_ss[:,i]'*UcpV[:,i]
end
Ucproj = Sψ.+Sλ
dif_transition = Uc.- VTp
findmin(abs.(dif_transition))

tauopt = TT./Y
tauopt[ind]
tauopt[ind_steady]


#############Plots
using LaTeXStrings
p1= plot(tauopt,Sψ ,label="Sψ", linewidth = 2, dpi=200, color = :blue)
p1= plot!(tauopt,VTp ,label="V'(T)", linewidth = 2, dpi=200,linestyle =:dash)
p1= plot!(tauopt,Uc ,label="U'(c)", linewidth = 2, dpi=200,linestyle =:dash)

using LaTeXStrings
p2= plot(tauopt,dif_truncation ,label="Δtruncation", linewidth = 2, dpi=200, color = :blue)
p2= plot!(tauopt,dif_transition ,label="Δtransition", linewidth = 2, dpi=200,linestyle =:dash)

Epsilon= ones(length(Welfare))
for i = 1: length(Welfare)
    Epsilon[i] = exp((1- β)*(Welfare[i] - Welfare[ind])) -1
end

Epsilon_steady= ones(length(Welfare))
for i = 1: length(Welfare)
    Epsilon_steady[i] = exp((1- β)*(Welfare_steady[i] - Welfare_steady[ind_steady])) -1
end

using LaTeXStrings
p3= plot(tauopt,100 .*Epsilon,label="Welfare_transition", linewidth = 2, dpi=200, color = :blue)
#p3= plot!(tau, Welfare_steady ,label="Welfare_ss", linewidth = 2, dpi=200,linestyle =:dash, legend=:bottomleft)
savefig(p3, "Welfare_optimal_low.png")

using LaTeXStrings
p4= plot(tauopt,100 .*Epsilon,label= "Consumption Equivalent", linewidth = 2, dpi=200, color = :blue)
#p4= plot!(tau, Welfare_steady ,label="Welfare_ss", linewidth = 2, dpi=200 ,legend=:bottomleft)
p4 = plot!([tauopt[ind],0.08], seriestype="vline", xticks = ([tauopt[ind],0.08],["\$ \\tau^{tr} = 0.064 \$", "\$ \\tau^{p} = 0.08 \$" ]), linewidth = 1, dpi=200,linestyle =:dash ,  label="")
savefig(p4, "Welfare_optimal_low_tau.png")

findmin(abs.(Welfare.-Welfare_steady))

using FileIO
save("Distribution_low.jld2", "Welfare", Welfare, "Welfare_steady", Welfare_steady, "ind", ind, "ind_steady", ind_steady)
