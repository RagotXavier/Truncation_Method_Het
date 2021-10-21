#########################################################
 # Defining Parameters
#########################################################

#############Parameters

α = 0.36 #capital-share
δ = 0.025 #depreciation in quarterly terms
γ = 1.0001 #concavity of the utility function
ns = 5 #number of grid points for the idiosyncratic process
na = 100 #number of grids points for the assets
amin = 1e-9 #minimum value for the assets
amax = 1000.0 #maximum value for the assets
curv = 4 #curvature for the asset grid
pers = 0.9617 #persistence component of income in yearly terms
sigma2 = sqrt(0.0396)
vary = (sigma2^2/(1 - pers^2)) #standard deviation component of income in yearly terms
rho    = pers^0.25; #persistence component of income in quarterly terms
sige = ((sigma2^2/(1 + rho^2 + rho^4 + rho^6 ))^(0.5)) #standard deviation component of income in quarterly terms

#############Variables in the steady-state

tk,tl,tc = 0.00,0.00,0.00 #tax on capital, labor, and consumption
KsY = 10.26 #capital to output ratio
#β = (1+KsY^(-1)*α-δ)^-1 #discount factor
β = 0.98 #discount factor
R = 1/β #before-tax gross interest rate
R = 1 + (1 - tk)*(R-1) #after-tax interest rate
w = (1-α)*( (R - (1-δ))/α )^(α/(α-1)) #before-tax gross wage
w =  w*(1-tl)/(1+tc) #after-tax wage corrected for consumption taxes
Ltot =1 #total labor

K0 = Ltot*( (R - (1-δ))/α )^(1/(α-1)) #capital
Y = K0^α*Ltot^(1-α) #output

TT = 0.0*Y #transfers
TT = TT/(1+tc) #transfers corrected for consumption taxes


#########################################################
 # Functions
#########################################################

#############Utility function
if γ==1.0
    Uvv(X) = log.(X).+10
    Ulinv(X) = exp.(X .- 10.)
else
    Uvv(X) =  (X.^(1-γ).-1)./(1-γ) .+10.
    Ulinv(X) = (1. .+(1-γ).* ( X .- 10.) ).^(1/(1-γ))
end


#############Derivative of the utility function
uPrime(c,γ) = c.^(-γ)
uPrimeInv(up,γ) = up.^(-1.0/γ)


#############Idiosyncratic process
mc     = rouwenhorst(ns,rho,sige)
endow  = exp.(mc.state_values)
Trans = collect(mc.p') #redefining the Transition Matrix by putting additional weight on the diagonal
Trans[findall(x->x<=5*10^-5,Trans)].=0
for i =1:ns
    Trans[i,i] = Trans[i,i]+ (1 - sum(Trans,dims=1)[i])
end

stdist = (Trans^100000)[:,1] #stationary distribution for the idosyncratic process
Transv = vec(Trans)


#############Renormalization
lbar   = dot(stdist,endow .^(1.0))
states = endow./ (lbar)
vstan = stdist; #vector of stationary distribution #dot(stdist,states) has to be equal to one
