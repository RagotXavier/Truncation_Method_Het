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
data = load("initial_dis.jld2")


using MAT
fileIn = matopen("tofigtruncation.mat")
dset = read(fileIn)
coef = dset["coef"]
TT = dset["TT_ss"]


using MAT
fileIn = matopen("pathwelfare.mat")
dset = read(fileIn)
GDP_path = dset["GDP"]
TT_path = dset["TT"]
C_path = dset["C"]


##############################
#   MAIN
##############################
include("Aiyagari_solve.jl")
include("Parameters.jl")

θ = 0.23620745544082986

Y0 = copy(Y0)

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


Grid = repeat(AA0.aGrid,ns)
D2 = copy(D_ss)
ww = [Grid D2]
ww2 = sortslices(ww, dims=1)
Gini = gini2(ww2)
vW = wealthD(ww2,5)

d_initial_ss = D_ss


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
    eps = 0.0
    Z[1,i]=1+eps
    K[1,i] = Ki[i]
    rhoz = 0
    for j=2:Periods
        Z[j,i] = 1+rhoz*(Z[j-1, i]-1)
        K[j,i] = Ki[i]
    end
    R[:,i],w[:,i],L[:,i] = Pricesv(Ki[i],K[:,i],Z[:,i],AA0.params,TT[i])
end


##############################
#   TRANSITION AND IMPULSE RESPONSE FUNCTIONS
##############################
include("MIT_functions_path.jl")

K_ss = A
TimePeriods = 400 #time periods for the simulation
OffeqPathTime = 1

z_ss = 1.0 #steady state productivity


Zpath  = ones(Float64,TimePeriods+1) #add time 0
epsz = 0.00
rhoz = 0.00

Zpath[1]=1+epsz
σf = 0

for i=2:TimePeriods
    Zpath[i] = 1+rhoz*(Zpath[i-1]-1)
end


ap =20 #number of alternative policies
TTpath = TT_path[2:end]

TTpath_ap1 = ones(Float64,TimePeriods+1, ap) #trajectories for TT in the alternative policies
for i = 1: ap
    TTpath_ap1[:,i] = (i).*(TTpath .- TTpath[end]) .+ TTpath[end]
    #TTpath_ap1[:,i] = (1+(i/10)).*(TTpath .- TTpath[end]) .+ TTpath[end]
end

TTpath_ap1 = reverse(TTpath_ap1, dims=2)

TTpath_ap2 = ones(Float64,TimePeriods+1, ap) #trajectories for TT in the alternative policies
for i = 1: ap
    TTpath_ap2[:,i] = -(i-1).*(TTpath .- TTpath[end]) .+ TTpath[end]
    #TTpath_ap2[:,i] = (1-(i/10)).*(TTpath .- TTpath[end]) .+ TTpath[end]
end

TTpath_ap = [TTpath_ap1 TTpath_ap2]

eq_apols = zeros(Float64,na*ns,TimePeriods+1,length(TTpath_ap[1,:]))
eq_cpols = zeros(Float64,na*ns,TimePeriods,length(TTpath_ap[1,:]))
eq_dpols = zeros(Float64,na*ns,TimePeriods+1,length(TTpath_ap[1,:]))
eq_Utpols = zeros(Float64,na*ns,TimePeriods,length(TTpath_ap[1,:]))
eq_aggKs = zeros(Float64,TimePeriods+1,length(TTpath_ap[1,:]))
Kguess = zeros(Float64,TimePeriods+1,length(TTpath_ap[1,:]))


#############Create capital path guess and path for shocks
for i = 1:length(TTpath_ap[1,:])
    Kguess[:,i] = K_ss.*ones(TimePeriods+1) #add time 0
end


for i = 1: length(TTpath_ap[1,:]) #creating one equilibrium for each value of TT
    TTpath = TTpath_ap[:,i]
    Kguess_new = copy(Kguess[:,i])
    polC_ss_new = copy(polC_ss[:,1])
    Ulevel_new = copy(Ulevel[:,1])
    polA_ss_new = copy(polA_ss[:,1])
    D_ss_new = copy(d_initial_ss)
    TT_new = copy(TTpath)

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


Rv = zeros(Float64,TimePeriods,length(eq_aggKs[1,:]))
wv = zeros(Float64,TimePeriods,length(eq_aggKs[1,:]))
Lv = zeros(Float64,TimePeriods,length(eq_aggKs[1,:]))
Yv = zeros(Float64,TimePeriods,length(eq_aggKs[1,:]))
Iv = zeros(Float64,TimePeriods,length(eq_aggKs[1,:]))
Tv = zeros(Float64,TimePeriods,length(eq_aggKs[1,:]))
Cv = zeros(Float64,TimePeriods,length(eq_aggKs[1,:]))
Cv2 = zeros(Float64,TimePeriods,length(eq_aggKs[1,:]))
CC = zeros(Float64,na*ns,length(eq_aggKs[1,:]))
ci = zeros(Float64,na*ns,length(eq_aggKs[1,:]))
Ul = zeros(Float64,na*ns,length(eq_aggKs[1,:]))
Ws = zeros(Float64,length(eq_aggKs[1,:]))
Wsv = zeros(Float64,TimePeriods,length(eq_aggKs[1,:]))
Z_eps= zeros(Float64,TimePeriods,length(eq_aggKs[1,:]))
GDP_eps= zeros(Float64,TimePeriods,length(eq_aggKs[1,:]))
C_eps= zeros(Float64,TimePeriods,length(eq_aggKs[1,:]))
K_eps= zeros(Float64,TimePeriods,length(eq_aggKs[1,:]))
I_eps= zeros(Float64,TimePeriods,length(eq_aggKs[1,:]))
W_eps= zeros(Float64,TimePeriods,length(eq_aggKs[1,:]))
R_eps= zeros(Float64,TimePeriods,length(eq_aggKs[1,:]))
w_eps= zeros(Float64,TimePeriods,length(eq_aggKs[1,:]))
T_eps= zeros(Float64,TimePeriods,length(eq_aggKs[1,:]))


for i = 1:length(eq_aggKs[1,:])
    Rv[:,i],wv[:,i],Lv[:,i] = Aggregates(eq_aggKs[2:end,i],TTpath_ap[1:end-1, i],Zpath[1:end-1],AA0.params);
    Yv[:,i] = Zpath[1:end-1].*Lv[:,i].^(1-α).*eq_aggKs[1:end-1,i].^α
    Iv[:,i] =  eq_aggKs[2:end,i] - (1 - AA0.params.δ)*eq_aggKs[1:end-1,i]
    Tv[:,i] =  TTpath_ap[1:end-1, i]
    Cv[:,i] =  Yv[:,i]- Iv[:,i] - Tv[:,i]
    Cv2[:,i] = wv[:,i].*Lv[:,i] + (Rv[:,i]).*eq_aggKs[1:end-1,i] - eq_aggKs[2:end,i]- Tv[:,i]

    #############Computing total Welfare in the Steady state
    CC[:,1]=(Mat[:,:,1]*D_ss[:,1]).*polC_ss[:,1]  #sum(CC) = Ctot
    global ii=findall(x -> x==0, D_ss[:,1])
    ci[:,1] = CC[:,1]./D_ss[:,1]
    Ul[:,1] =Uvv(ci[:,1])
    Ul[ii,1] .=0
    Ws[i] =sum(D_ss[:,1].*(Ul[:,1])) .+ TTpath_ap[end, i].^(θ)

    #############Computing the Welfare in the Transition
    for j=1:TimePeriods
        #Wsv[j, i] = sum(eq_dpols[:,j+1, i].*(eq_Utpols[:,j,i])).+ TT[i].^(θ)
        Wsv[j, i] = sum(eq_dpols[:,j+1, i].*(eq_Utpols[:,j,i])).+ TTpath_ap[j, i].^(θ)
    end

    #############IMPULSE RESPONSE FUNCTIONS
    Z_eps[:,i] = 100*(Zpath[1:end-1].-1)
    GDP_eps[:,i] =100*(Yv[:,i].-Y[1])./Y[1]
    C_eps[:,i] = 100*(Cv[:,i].- Ctot_t[1])/ Ctot_t[1]
    K_eps[:,i] = 100*(eq_aggKs[2:end,i].- A_t[1])./A_t[1]
    I_eps[:,i] = 100*(Iv[:,i].-Iss[1])./Iss[1]
    W_eps[:,i] = 100*(Wsv[:,i] .- Ws[1])./Ws[1]
    R_eps[:,i] = 100*(Rv[:,i] .- Prices[2,1])
    w_eps[:,i] = 100*(wv[:,i] .- Prices[1,1])./Prices[1,1]
    T_eps[:,i] = 100*(Tv[:,i].-(Tv[end,i]))./(Tv[end,i])
end



##############################
#   COMPUTING THE WELFARE IN THE TRANSITION
##############################

β_vector = zeros(Float64,length(Wsv[:,1]))
for i = 1:length(Wsv[:,1])
    β_vector[i] = β.^(i-1)
end
Welfare = transpose(β_vector'*Wsv) .+ (β^(length(Wsv[:,1]))./(1 - β)).*Wsv[end,:]
Welfare_max, ind = findmax(Welfare)


Epsilon= ones(length(Welfare))
for i = 1: length(Welfare)
    Epsilon[i] = exp((1- β)*(Welfare[i] - Welfare[ap])) -1
end

plot(-(ap-1):1:ap, 100 .*Epsilon,
xlabel="Deviation from the optimal path", ylabel="Consumption equivalent", label = "")
savefig("CE.png")

100*(TTpath_ap[1,ind] - TTpath_ap[1,ap])./TTpath_ap[1,ap] #deviation in percentage terms in the taxes in period 1

tau_ind = TTpath_ap[1:end-1,ind]./Yv[:, ind] #optimal tax path in transition
tau_ap = TTpath_ap[1:end-1,ap]./Yv[:, ap] #optimal tax path in projection


index = ones(length(TTpath_ap[1,:]))

for i = 1:length(TTpath_ap[1,:])
    index[i] = 100*(TTpath_ap[1,i] - TTpath_ap[1,ap])./TTpath_ap[1,ap]
end

plot(index, 100 .*Epsilon,
xlabel="Percentage deviation from the optimal path", ylabel="Consumption equivalent", label = "")
savefig("CE_r.png")

plot(tau_ap, label="τ projection", lw = 2, legend=:bottomright)
plot!(tau_ind, label="τ transition", lw = 2, legend=:bottomright)
savefig("tau_r.png")

tau_real = TT_path[2:end-1]./Yv[:, ap]
plot(tau_ap, label="τ projection estimated", lw = 2, legend=:bottomright)
plot!(tau_real, label="τ projection real", lw = 2, legend=:bottomright)
savefig("tau_estimated_real.png")


Epsilon2= ones(length(Welfare))
for i = 1: length(Welfare)
   Epsilon2[i] = exp((Welfare[i] - Welfare[ap]))-1
end

plot(index, 100*Epsilon2,
xlabel="Deviation from the optimal path", ylabel="Consumption equivalent", label = "")
savefig("CE_r2.png")


using FileIO
save("CE.jld2", "Welfare", Welfare, "ind", ind, "TTpath_ap", TTpath_ap, "ap", ap, "tau_ind", tau_ind, "tau_ap", tau_ap, "Yv", Yv)
