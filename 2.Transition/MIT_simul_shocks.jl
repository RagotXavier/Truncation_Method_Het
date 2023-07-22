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

θ = 0.23620745544082986

Y0 = copy(Y0)
tau = TT/Y0

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
    eps = 0.01
    #eps = 0.0
    Z[1,i]=1+eps
    K[1,i] = Ki[i]
    rhoz = 0.95
    #rhoz = 0
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
epsz = 0.00312 #standard deviation for the shock
rhoz = 0.95

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
    W_eps[:,i] = 100*(exp.(Wsv[:,i] .- Ws[i]) .-1)
    R_eps[:,i] = 100*(Rv[:,i] .- Prices[2,i])
    w_eps[:,i] = 100*(wv[:,i] .- Prices[1,i])./Prices[1,i]
    T_eps[:,i] = 100*(Tv[:,i].-(Tv[end,i]))./(Tv[end,i])
end


#############READ MATLAB FILES FOR THE IRF AND PLOT THE IRF COMPARISON
using MAT
fileIn = matopen("tofigtruncation.mat")
dset = read(fileIn)
close(fileIn)
Zt_epsT = dset["Zt_eps"]
GDPt_epsT= dset["GDPt_eps"]
Ct_epsT = dset["Ct_eps"]
Kt_epsT= dset["Kt_eps"]
It_epsT= dset["It_eps"]
Wt_epsT= dset["Wt_eps"]
Rt_epsT= dset["Rt_eps"]
wt_epsT= dset["wt_eps"]
Tt_epsT= dset["Tt_eps"]

using MAT
fileIn = matopen("tofigreiter.mat")
dset = read(fileIn)
close(fileIn)
Zt_epsR = dset["Zt_eps"]
GDPt_epsR= dset["GDPt_eps"]
Ct_epsR = dset["Ct_eps"]
Kt_epsR= dset["Kt_eps"]
It_epsR= dset["It_eps"]
Wt_epsR= dset["Wt_eps"]
Rt_epsR= dset["Rt_eps"]
wt_epsR= dset["wt_eps"]
Tt_epsR= dset["Tt_eps"]

using MAT
fileIn = matopen("tofigRA.mat")
dset = read(fileIn)
close(fileIn)
Zt_epsRA = dset["Zt_eps"]
GDPt_epsRA= dset["GDPt_eps"]
Ct_epsRA = dset["Ct_eps"]
Kt_epsRA= dset["Kt_eps"]
It_epsRA= dset["It_eps"]
Wt_epsRA= dset["Wt_eps"]
Rt_epsRA= dset["Rt_eps"]
wt_epsRA= dset["wt_eps"]
Tt_epsRA= dset["Tt_eps"]


using LaTeXStrings
p1= plot(Zt_epsT ,label="Truncation", linewidth = 2, dpi=300, color = :blue)
p1= plot!(Zt_epsR ,label="Reiter", linewidth = 2, dpi=300,linestyle =:dash)
#p1= plot!(Zt_epsRA ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
#p1= plot!(Z_eps[1:TimePeriods], label="BKM", linewidth = 2, dpi=300, linestyle =:dash)
title!("\$Z\$")

p2= plot(GDPt_epsT ,label="Truncation", linewidth = 2, dpi=300, color = :blue)
p2= plot!(GDPt_epsR ,label="Reiter", linewidth = 2, dpi=300,linestyle =:dash)
#p2= plot!(GDPt_epsRA ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
#p2= plot!(GDP_eps[1:TimePeriods], label="BKM", linewidth = 2, dpi=300, linestyle =:dash)
title!("\$GDP\$")

p3= plot(Ct_epsT ,label="Truncation", linewidth = 2, dpi=300, color = :blue)
p3= plot!(Ct_epsR ,label="Reiter", linewidth = 2, dpi=300,linestyle =:dash)
#p3= plot!(Ct_epsRA ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
#p3= plot!(C_eps[1:TimePeriods], label="BKM", linewidth = 2, dpi=300, linestyle =:dash)
title!("\$C\$")

p4= plot(Kt_epsT ,label="Truncation", linewidth = 2, dpi=300, color = :blue)
p4= plot!(Kt_epsR ,label="Reiter", linewidth = 2, dpi=300,linestyle =:dash)
#p4= plot!(Kt_epsRA ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
#p4= plot!(K_eps[1:TimePeriods], label="BKM", linewidth = 2, dpi=300, linestyle =:dash)
title!("\$K\$")

p5= plot(It_epsT ,label="Truncation", linewidth = 2, dpi=300, color = :blue)
p5= plot!(It_epsR ,label="Reiter", linewidth = 2, dpi=300,linestyle =:dash)
#p5= plot!(It_epsRA ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
#p5= plot!(I_eps[1:TimePeriods], label="BKM", linewidth = 2, dpi=300, linestyle =:dash)
title!("\$I\$")

p6= plot(Wt_epsT ,label="Truncation", linewidth = 2, dpi=300, color = :blue)
p6= plot!(Wt_epsR ,label="Reiter", linewidth = 2, dpi=300,linestyle =:dash)
#p6= plot!(Wt_epsRA ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
#p6= plot!(W_eps[1:TimePeriods], label="BKM", linewidth = 2, dpi=300, linestyle =:dash)
title!("\$Welfare\$")

p7= plot(Rt_epsT ,label="Truncation", linewidth = 2, dpi=300, color = :blue)
p7= plot!(Rt_epsR ,label="Reiter", linewidth = 2, dpi=300,linestyle =:dash)
#p7= plot!(Rt_epsRA ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
#p7= plot!(R_eps[1:TimePeriods], label="BKM", linewidth = 2, dpi=300, linestyle =:dash)
title!("\$R\$")

p8= plot(wt_epsT ,label="Truncation", linewidth = 2, dpi=300, color = :blue)
p8= plot!(wt_epsR ,label="Reiter", linewidth = 2, dpi=300,linestyle =:dash)
#p8= plot!(wt_epsRA ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
#p8= plot!(w_eps[1:TimePeriods], label="BKM", linewidth = 2, dpi=300, linestyle =:dash)
title!("\$w\$")

p9= plot(Tt_epsT ,label="Truncation", linewidth = 2, dpi=300, color = :blue)
p9= plot!(Tt_epsR ,label="Reiter", linewidth = 2, dpi=300,linestyle =:dash)
#p9= plot!(Tt_epsRA ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
#p9= plot!(T_eps[1:TimePeriods], label="BKM", linewidth = 2, dpi=300, linestyle =:dash)
title!("\$T\$")

plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, layout=(3,3), size=(750, 750), dpi=300)
savefig("IRF")
