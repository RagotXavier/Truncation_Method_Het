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

K0= 44.98314527071115
AA0,pol0,Aggs0 = AiyagariEGM(K0,ns,Transv,states,amin,amax,curv,β,α,δ,γ,TT,na)
polA_ss,polC_ss,D_ss,A,Aggs_ss,Mat,upc,Ls,Ulevel = steady(K0,pol0,AA0,R,w,TT,1e-4)

Ltot = dot(vstan,states) #total labor
A = (1+tc)*A #total assets corrected for consumption taxes
B = 0 #total debt
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


##############################
#   REITER'S METHOD
##############################
aGrid = AA0.aGrid
Vind,Wind,resE,apol = Projection_Reiter( AA0,polA_ss, polC_ss,Aggs_ss,Mat,states,Ls)

Lv = kron(Ls[1:ns],ones(na,1)) #vector of labor supplies
ytype = repeat(kron(collect(StepRange(1, Int8(1), ns)),ones(Int,na))) #vector of individual's types
Cpol =  - apol' + Aggs_ss.R*repeat(AA0.aGrid,ns) + Aggs_ss.w*states[ytype].*Lv #policy function for consumption
RR = repeat(aGrid,ns) + polC_ss - Aggs_ss.R*polA_ss - Aggs_ss.w*states[ytype].*Lv #vector of residuals, it has to be zero
plot(resE)

#############Computing total Welfare in the Steady state
CC=(Mat*D_ss).*Cpol  #sum(CC) = Ctot
ii=findall(x -> x==0, D_ss)
ci = CC./D_ss
Ul =Uvv(ci)
Ul[ii] .=0

Ws =sum(D_ss.*Ul)


#############Determining the consumption in percentiles of the population
ww = [dGridl_t D2 collect(1:na*ns)]
ww2 = sortslices(ww, dims=1)

fracS2 = cumsum(ww2[:,2]) #determine the cumulative distribution of the population

b = findall(x -> x<=0.10, fracS2) #bottom 10 percent
b = Int.(ww2[:,3][b]) #find the points where the bottom 10 percent are located
Cbsa = D_ss[b]'*polC_ss[b] #total consumption of the bottom 10 percent
Cbs = Cbsa/sum(D_ss[b]) #per capita consumption of the bottom 10 percent

ww3= ww2[sortperm(ww2[:, 1],rev=true), :]
fracS3 = cumsum(ww3[:,2]) #determine the inverse cumulative distribution of the population

m= findall(x -> x <=0.50, fracS3) #top 50 percent
m = Int.(ww3[:,3][m]) #find the points where the top 50 percent are located
Cmsa = D_ss[m]'*polC_ss[m] #total consumption of the top 50 percent
Cms = Cmsa/sum(D_ss[m]) #per capita consumption of the top 50 percent

t= findall(x -> x <=0.10, fracS3) #top 10 percent
t = Int.(ww3[:,3][t]) #find the points where the top 10 percent are located
Ctsa = D_ss[t]'*polC_ss[t] #total consumption of the top 10 percent
Cts = Ctsa/sum(D_ss[t]) #per capita consumption of the top 10 percent


##############################
#   Saving to Dynare
##############################
file = matopen("todynare_Reiter.mat", "w")
write(file, "alpha", AA0.params.α)
write(file, "beta", AA0.params.β)
write(file, "delta", AA0.params.δ)
write(file, "gamma", AA0.params.γ)
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
write(file, "m",m)
write(file, "Cms",Cms)
write(file, "t",t)
write(file, "Cts",Cts)
write(file, "Mat",Mat)
write(file, "Ws",Ws)
close(file)


##############################
#  Policy Functions
##############################

policyA = reshape(polA_ss,na,ns) #policy function for assets
policyVa = reshape(apol,na,ns) #policy function for assets
policyC=reshape(Cpol,na,ns) #policy function for consumption
policyL=reshape(Lv,na,ns) #policy function for labor
Dist = reshape(D_ss,na,ns) #policy function for distribution


#############Policy function for assets
p1 = plot(aGrid,policyA[:,1],title="Policy rule Assets",label="Pol 1",titlefont=font(7, "Courier"))
p1 = plot!(aGrid,policyA[:,2],title="Policy rule Assets",label="Pol 2",titlefont=font(7, "Courier"))
p1 = plot!(aGrid,policyA[:,3],title="Policy rule Assets",label="Pol 3",titlefont=font(7, "Courier"))
p1 = plot!(aGrid,policyA[:,4],title="Policy rule Assets",label="Pol 4",titlefont=font(7, "Courier"))
p1 = plot!(aGrid,policyA[:,5],title="Policy rule Assets",label="Pol 5",titlefont=font(7, "Courier"))

p1 = plot(aGrid,policyVa[:,1],title="Policy rule Assets",label="Pol 1",titlefont=font(7, "Courier"))
p1 = plot!(aGrid,policyVa[:,2],title="Policy rule Assets",label="Pol 2",titlefont=font(7, "Courier"))
p1 = plot!(aGrid,policyVa[:,3],title="Policy rule Assets",label="Pol 3",titlefont=font(7, "Courier"))
p1 = plot!(aGrid,policyVa[:,4],title="Policy rule Assets",label="Pol 4",titlefont=font(7, "Courier"))
p1 = plot!(aGrid,policyVa[:,5],title="Policy rule Assets",label="Pol 5",titlefont=font(7, "Courier"))
#savefig(p1,"Policy_rule_assets.pdf")


#############Policy function for consumption
p2 = plot(aGrid,policyC[:,1],title="Policy rule Consumption",label="Pol 1",titlefont=font(7, "Courier"))
p2 = plot!(aGrid,policyC[:,2],title="Policy rule Consumption",label="Pol 2",titlefont=font(7, "Courier"))
p2 = plot!(aGrid,policyC[:,3],title="Policy rule Consumption",label="Pol 3",titlefont=font(7, "Courier"))
p2 = plot!(aGrid,policyC[:,4],title="Policy rule Consumption",label="Pol 4",titlefont=font(7, "Courier"))
p2 = plot!(aGrid,policyC[:,5],title="Policy rule Consumption",label="Pol 5",titlefont=font(7, "Courier"))
#savefig(p2,"Policy_rule_consumption.pdf")


#############Policy function for labor
p3 = plot(aGrid,policyL[:,1],title="Policy rule Labor",label="Pol 1",titlefont=font(7, "Courier"))
p3 = plot!(aGrid,policyL[:,2],title="Policy rule Labor",label="Pol 2",titlefont=font(7, "Courier"))
p3 = plot!(aGrid,policyL[:,3],title="Policy rule Labor",label="Pol 3",titlefont=font(7, "Courier"))
p3 = plot!(aGrid,policyL[:,4],title="Policy rule Labor",label="Pol 4",titlefont=font(7, "Courier"))
p3 = plot!(aGrid,policyL[:,5],title="Policy rule Labor",label="Pol 5",titlefont=font(7, "Courier"))
#savefig(p3,"Policy_rule_labor.pdf")


#############Policy function for distribution
p4 = plot(aGrid,Dist[:,1],title="Distribution",label="Pol 1",titlefont=font(7, "Courier"))
p4 = plot!(aGrid,Dist[:,2],title="Distribution",label="Pol 2",titlefont=font(7, "Courier"))
p4 = plot!(aGrid,Dist[:,3],title="Distribution",label="Pol 3",titlefont=font(7, "Courier"))
p4 = plot!(aGrid,Dist[:,4],title="Distribution",label="Pol 4",titlefont=font(7, "Courier"))
p4 = plot!(aGrid,Dist[:,5],title="Distribution",label="Pol 5",titlefont=font(7, "Courier"))
#savefig(p4,"Distribution.pdf")
