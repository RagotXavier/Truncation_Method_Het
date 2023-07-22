#########################################################
 # Main file to run the results using Reiter's Method
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
using Roots


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


struct Projection{T <: Real,I <: Integer} #be careful the following structure store information by bin (and not per capita). Divide by the size of the bin to have per capita terms
    Agg::AggVarsEGM{T,T}
    N::I #lenght of the truncation
    ns::I #number of states
    Nbin::I #number of bins
    R::T
    w::T
    A::T
    B::T
    states::Array{T,1}
    Transv::Array{T,1}
    Sp::Array{T,1} #size of each bin
    abp::Array{T,1}
    aep::Array{T,1}
    lb::Array{T,1}
    cp::Array{T,1}
    Ucv::Array{T,1} #vector of marginal utilities (before multiplication by xsip) by bin
    Matab::SparseMatrixCSC{T,I} #transition matrix
    Ucp::Array{T,1}
    Ul::Array{T,1}
    Res::Array{T,1} #check that res is small and Res =  Upvb - β*R*Matab*Upvb
    xsip::Array{T,1}
    Ntoi::Array{I,1}
    CC::Array{T,1}
    indnc::Array{I,1} #index of non-credit constrained histories
    indcc::Array{I,1} #index of credit constrained histories
    ytype::Array{I,1}
    xsyn::Array{Float64,1}
    xsyn1::Array{Float64,1}
    xsyn2::Array{Float64,1}
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
include("Projection_Reiter.jl")
include("Parameters.jl")

tau = 0.082
θ = 0.23620745544082986

TT = tau.*Y0

AA0,pol0,Aggs0 = AiyagariEGM(K0,ns,Transv,states,amin,amax,curv,β,α,δ,γ,TT,na)
polA_ss,polC_ss,D_ss,A,Aggs_ss,Mat,upc, Ls, Ulevel = steady(K0,pol0,AA0,R,w,TT,1e-4)

Ltot = dot(vstan,states) #total labor
B = A-A #total debt in the economy
Ctot = D_ss'*polC_ss #total consumption
Y = A^α*Ltot^(1-α) #total output
G = Y- δ*A - Ctot #total government expenses, it has to be close to TT

#Verifying
A_t= ((Aggs_ss.R-1 + δ)/(α*Ltot^(1-α)))^(1 /(α-1))
Ctot_t = Aggs_ss.w*Ltot + (Aggs_ss.R-1)*A - TT
tauopt = TT./Y

@unpack dGrid = AA0
dGridl_t = repeat(dGrid,ns)
D2 = copy(D_ss)
D2[1] = D_ss[1]
sum(D2)
ww = [dGridl_t D2]
ww2 = sortslices(ww, dims=1)
Gini = gini2(ww2)
vW = wealthD(ww2,5)

na = AA0.na
ns = AA0.ns
cpol = reshape(polC_ss,na,ns)
grida = reshape(dGridl_t,na,ns)
dist = reshape(D_ss,na,ns)
dc = cpol[2:na,:] .- cpol[1:na-1,:]
da = grida[2:na,:].- grida[1:na-1,:]
mpc = dc./da
ampc = sum(mpc.*(dist[1:na-1,:]))

println("K/Y: ",A/(Y*4), " Gini: ", Gini, " G/Y ",G/Y, " C/Y ",Ctot/Y, " I/Y ",δ*A/Y,  " MPC ", ampc)
vW = wealthD(ww2,10)


##############################
#   REITER'S METHOD
##############################
aGrid = AA0.aGrid
Vind,Wind,resE,apol = Projection_Reiter(AA0,polA_ss, polC_ss,Aggs_ss,Mat,states,Ls)

Lv = kron(Ls[1:ns],ones(na,1)) #vector of labor supplies
ytype = repeat(kron(collect(StepRange(1, Int8(1), ns)),ones(Int,na))) #vector of individual's types
Cpol =  - apol' + Aggs_ss.R*repeat(AA0.aGrid,ns) + Aggs_ss.w*states[ytype].*Lv- Aggs_ss.TT*ones(Float64,na*ns) #policy function for consumption
RR = repeat(aGrid,ns) + polC_ss - Aggs_ss.R*polA_ss - Aggs_ss.w*states[ytype].*Lv + Aggs_ss.TT*ones(Float64,na*ns)  #vector of residuals, it has to be zero
plot(resE)

#############Computing total Welfare in the Steady state
CC=(Mat*D_ss).*polC_ss  #sum(CC) = Ctot
ii=findall(x -> x==0, D_ss)
ci = CC./D_ss
Ul =Uvv(ci)
Ul[ii] .=0

Ws =sum(D_ss.*Ul) + TT^(θ)

#############Determining the consumption in percentiles of the population
ww = [dGridl_t D2 collect(1:na*ns)]
ww2 = sortslices(ww, dims=1)

fracS2 = cumsum(ww2[:,2]) #determine the cumulative distribution of the population

b = findall(x -> x<=0.10, fracS2) #bottom 10 percent
b = Int.(ww2[:,3][b]) #find the points where the bottom 10 percent are located
Cbsa = D_ss[b]'*polC_ss[b] #total consumption of the bottom 10 percent
Cbs = Cbsa/sum(D_ss[b]) #per capita consumption of the bottom 10 percent


##############################
#   Saving to Dynare
##############################
file = matopen("todynare_Reiter.mat", "w")
write(file, "alpha", AA0.params.α)
write(file, "beta", AA0.params.β)
write(file, "delta", AA0.params.δ)
write(file, "gamma", AA0.params.γ)
write(file, "theta", θ)
write(file, "tau", tau)
write(file, "Trans", reshape(AA0.Transv,ns,ns))
write(file, "abar", minimum(AA0.aGrid))
write(file, "states", AA0.states)
write(file, "TT", Aggs_ss.TT)
write(file, "Vind", Vind)
write(file, "Wind", Wind)
write(file, "resE", resE)
write(file, "na", AA0.na)
write(file, "ns", AA0.ns)
write(file, "ytype", ytype)
write(file, "apol", apol.parent) # policy rule ap(a)
write(file, "aGrid", repeat(AA0.aGrid,ns))
write(file, "D_ss", D_ss)
write(file, "Cpol",Cpol)
write(file, "polA_ss",polA_ss)
write(file, "w",Aggs_ss.w)
write(file, "R",Aggs_ss.R)
write(file, "K",A)
write(file, "Ltot",Ltot)
write(file, "Ctot",Ctot)
write(file, "b",b)
write(file, "Cbs",Cbs)
write(file, "Mat",Mat)
write(file, "Ws",Ws)
close(file)
