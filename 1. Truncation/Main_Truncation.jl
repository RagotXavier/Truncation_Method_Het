#########################################################
 # Main file to run the results using Truncation method
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
    #display(plot(fracS,fracw))
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


#############Variables for the Planner
struct Plan{T <: Real}
    lambda::Array{T,1}
    lambdat::Array{T,1}
    psih::Array{T,1}
    omega::Array{T,1} #weights
    Uv::Array{T,1}
    CC::Array{T,1} # vector of residuals, it must be 0 to check equations
    mu::T
end


##############################
#   MAIN
##############################
include("Aiyagari_solve.jl")
include("Projection_Truncation.jl")
include("Function_Optim.jl")
#include("Function_Optim2.jl") #such that the Pareto weights are all positive
include("Parameters.jl")

K0= 29.34474

AA0,pol0,Aggs0 = AiyagariEGM(K0,ns,Transv,states,amin,amax,curv,β,α,δ,γ,TT,na)
polA_ss,polC_ss,D_ss,A,Aggs_ss,Mat,upc = steady(K0,pol0,AA0,R,w,TT,1e-4)

Ltot = dot(vstan,states) #total labor
A = (1+tc)*A #total assets corrected for consumption taxes
B = A-A #total debt in the economy
Ctot = D_ss'*polC_ss #total consumption, it has to be equal to Ctot = Aggs_ss.w*Ltot + (Aggs_ss.R-1)*A + TT
Y = A^α*Ltot^(1-α) #total output
G = Y- δ*A - Ctot #total government expenses, it has to be close to zero

@unpack dGrid = AA0
dGridl_t = repeat(dGrid,ns)*(1+tc)
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

#
##############################
#   TRUNCATION METHOD
##############################

@time eco,upc,resE, stdΓ,resi,ytype = Projection_plan(N,AA0,polA_ss,polC_ss,upc,D_ss,Aggs_ss,Mat,A,B,vstan,states)
display(plot(eco.xsip))
display(plot(eco.xsyn))
display(plot(eco.xsyn2))
stdΓ

#############Computing total Welfare in the Steady state
Ws = sum(eco.xsyn.*eco.Sp.*eco.Ul)
Ws2 = sum(eco.xsyn2.*eco.Sp.*eco.Ul)


#############Computing the Pareto weights - call it Pl.omega
Pl = Optim(eco,AA0.params)


#############Computing the Pricing Kernel
M = sum(eco.Ucv.*eco.xsip.*Pl.omega.*eco.Sp)


#############Pareto weights versus Productivity
plot(eco.ytype,Pl.omega.*eco.Sp)
weight = zeros(Float64,eco.ns)
for i=1:eco.Nbin
    weight[eco.ytype[i]]=    weight[eco.ytype[i]]+eco.Sp[i]*Pl.omega[i]
end
weight=weight./stdist
plot(states,weight, ylims = (0,2),label="", linewidth = 2, dpi=300)
xlabel!("Productivity level")
ylabel!("Pareto weights")


#############Intertemporal Pareto weights
PI= Matrix(eco.Matab)
II = sparse(I,eco.Nbin,eco.Nbin)
II = Matrix(II)
Vp =inv(II - β*PI')*(Pl.omega.*eco.Ul)
Va= inv(II - β*PI')*eco.Ul
omega=Vp./Va
weight2 = zeros(Float64,eco.ns)
for i=1:eco.Nbin
    weight2[eco.ytype[i]]=    weight2[eco.ytype[i]]+eco.Sp[i]*omega[i]
end
weight2=weight2./stdist
plot(states,weight2, ylims = (0,2),label="", linewidth = 2, dpi=300)
xlabel!("Productivity level")
ylabel!("Pareto weights")


#############Determining the consumption in percentiles of the population
ww = [eco.abp./eco.Sp eco.Sp collect(1:eco.Nbin)]
ww2 = sortslices(ww, dims=1)

fracS2 = cumsum(ww2[:,2]) #determine the cumulative distribution of the population

b = findall(x -> x<0.10, fracS2) #bottom 10 percent
b = Int.(ww2[:,3][b]) #find the points where the bottom 10 percent are located
Cbsa = sum(eco.cp[b]) #total consumption of the bottom 10 percent
Cbs = Cbsa./sum(eco.Sp[b])

ww3= ww2[sortperm(ww2[:, 1],rev=true), :]
fracS3 = cumsum(ww3[:,2])

m= findall(x -> x < 0.50, fracS3) #top 50 percent
m = Int.(ww3[:,3][m]) #find the points where the top 50 percent are located
Cmsa = sum(eco.cp[m]) #total consumption of the top 50 percent
Cms = Cmsa./sum(eco.Sp[m])

t= findall(x -> x < 0.10, fracS3) #top 10 percent
t = Int.(ww3[:,3][t]) #find the points where the top 10 percent are located
Ctsa = sum(eco.cp[t]) #total consumption of the top 10 percent
Cts = Ctsa./sum(eco.Sp[t])


##############################
#   Saving to Dynare
##############################
file = matopen("todynare_Truncation.mat", "w")
write(file, "Pl", Pl)
write(file, "eco", eco)
write(file, "alpha", AA0.params.α)
write(file, "beta", AA0.params.β)
write(file, "delta", AA0.params.δ)
write(file, "gamma", AA0.params.γ)
write(file, "abar", minimum(AA0.aGrid))
write(file, "states", AA0.states)
write(file, "G", G)
write(file, "tk", tk)
write(file, "tl", tl)
write(file, "Ctot", Ctot)
write(file, "b",b)
write(file, "Cbs",Cbs)
write(file, "m",m)
write(file, "Cms",Cms)
write(file, "t",t)
write(file, "Cts",Cts)
write(file, "M", M)
write(file, "Ws", Ws)
write(file, "Ws2", Ws2)
close(file)
