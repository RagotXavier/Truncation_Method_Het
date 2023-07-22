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
p3= plot!(tauopt, 100 .*Epsilon_steady ,label="Welfare_ss", linewidth = 2, dpi=200, legend=:bottomleft)
savefig(p3, "Welfare_optimal_steady.png")

data = load("Distribution_low.jld2")
ind_low = data["ind"]

data = load("Distribution_large.jld2")
ind_large = data["ind"]


using LaTeXStrings
p4= plot(tauopt,100 .*Epsilon,label="Consumption Equivalent", linewidth = 2, dpi=200, legend=:topleft,color = :blue, foreground_color_legend = nothing,background_color_legend=nothing)
#p4= plot!(tau, Welfare_steady ,label="Welfare_ss", linewidth = 2, dpi=200, legend=:bottomleft)
p4 = plot!([tauopt[ind_low], tauopt[ind],0.0802, tauopt[ind_large]],xticks = ([tauopt[ind_low], tauopt[ind], 0.0802, tauopt[ind_large]],["\$ \\tau^{low} = 6.4 \$","\$ \\tau^{c} = 7.8 \\quad \\quad \$","\$\\quad \\quad \\quad \\tau^{p}_{\\infty} = 8\$", "\$ \\quad \\quad \\quad\\tau^{high} = 8.45 \$"]), seriestype="vline", linewidth = 1, dpi=200,linestyle =:dash ,  label="")
#legend(p4,frameon=false)
savefig(p4, "Welfare_optimal.png")

findmin(abs.(Welfare.-Welfare_steady))


using FileIO
save("Distribution.jld2", "Welfare", Welfare, "Welfare_steady", Welfare_steady, "ind", ind, "ind_steady", ind_steady)
