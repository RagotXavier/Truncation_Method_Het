#########################################################
 # Projection for Reiter's Method
#########################################################

function Projection_Reiter(
        AA0::AiyagariModelEGM,
        polA_ss::Array{T,1},
        polC_ss::Array{T,1},
        Aggs_ss::AggVarsEGM{T,T},
        Mat::Array{T,2},
        states::Array{T,1},
        Ls::Array{T,1}) where {T <: Real}

        Vmat = Mat[1:na,:]
        aGrid = AA0.aGrid
        aGrid[1] = 0
        Va = repeat(aGrid,ns)'*Mat # policy rule in a a'(a)
        pol = reshape(Va,na,ns)
        Vind = zeros(Int,na*ns)
        Wind = zeros(Float64,na*ns)

        for i = 1:length(Va)
                indcc=findall(x->x<=Va[i],aGrid)
                Vind[i] = maximum(indcc)
                Wind[i] = 1-(Va[i] - aGrid[Vind[i]])/(aGrid[Vind[i]+1]-aGrid[Vind[i]])
        end

        γ =  AA0.params.γ
        β =  AA0.params.β
        R = Aggs_ss.R
        Trans   = reshape(AA0.Transv,ns,ns)

        Ls  = fill(dot(vstan,states),ns) #vector of labor supplies
        Lv = kron(Ls[1:ns],ones(na,1))

        ytype = repeat(kron(collect(StepRange(1, Int8(1), ns)),ones(Int,na))) #vector of individual's types
        Cpol =  - Va' + Aggs_ss.R*repeat(aGrid,ns) + Aggs_ss.w*states[ytype].*Lv

        Xe = (Aggs_ss.w .*states)
        Xc = 0*Xe
        Xc = kron(Xc,ones(na,1))

        resE = zeros(Float64,na*ns)
        cimp = zeros(Float64,na*ns)
        for i=1:na
                for j=1:ns
                        ind = i+(j-1)*na
                        EUc = 0
                        for jp=1:ns
                                EUc = EUc + Trans[jp,j]*(
                                (Wind[ind]*Cpol[Vind[ind]+(jp-1)*na] +
                                (1-Wind[ind])*Cpol[Vind[ind]+1+(jp-1)*na] -
                                Xc[Vind[ind]+(jp-1)*na])^-γ)
                        end
                        resE[ind] = (Cpol[ind] -Xc[ind] )^-γ - β*R*EUc
                        cimp[ind] = (β*R*EUc)^(-1/γ)+Xc[ind]
                end
        end
        return Vind,Wind,resE,Va
end
