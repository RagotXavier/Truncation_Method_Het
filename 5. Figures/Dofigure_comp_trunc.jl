 ###########READ MATLAB FILES FOR THE IRF AND PLOT THE IRF COMPARISON
using MAT
using Plots; pyplot();



TimePeriods = 200 #time periods for the simulation

#
fileIn = matopen("../1. Truncation/MATLAB_DYNAMICS/tofigtruncation.mat")
#fileIn = matopen("tofigtruncation.mat")
dset = read(fileIn)
close(fileIn)
mut_epsT = dset["ut_eps"]
Zt_epsT = dset["Zt_eps"]
GDPt_epsT= dset["GDPt_eps"]
Ct_epsT = dset["Ct_eps"]
Cbt_epsT = dset["Cbt_eps"]
Cmt_epsT = dset["Cmt_eps"]
Ctt_epsT = dset["Ctt_eps"]
Kt_epsT= dset["Kt_eps"]
It_epsT= dset["It_eps"]
Wt_epsT= dset["Wt_eps"]
Rt_epsT= dset["rt_eps"]
wt_epsT= dset["wt_eps"]


fileIn = matopen("../1. Truncation/MATLAB_DYNAMICS/tofigtruncation2.mat")
#fileIn = matopen("tofigtruncation.mat")
dset = read(fileIn)
close(fileIn)
mut_epsT2 = dset["ut_eps"]
Zt_epsT2 = dset["Zt_eps"]
GDPt_epsT2= dset["GDPt_eps"]
Ct_epsT2 = dset["Ct_eps"]
Cbt_epsT2 = dset["Cbt_eps"]
Cmt_epsT2 = dset["Cmt_eps"]
Ctt_epsT2 = dset["Ctt_eps"]
Kt_epsT2= dset["Kt_eps"]
It_epsT2= dset["It_eps"]
Wt_epsT2= dset["Wt_eps"]
Rt_epsT2= dset["rt_eps"]
wt_epsT2= dset["wt_eps"]


fileIn = matopen("../4. RA/tofigRa.mat")
dset = read(fileIn)
close(fileIn)
mut_epsRa = dset["ut_eps"]
Zt_epsRa = dset["Zt_eps"]
GDPt_epsRa = dset["GDPt_eps"]
Ct_epsRa = dset["Ct_eps"]
Kt_epsRa = dset["Kt_eps"]
It_epsRa = dset["It_eps"]
Wt_epsRa = dset["Wt_eps"]
Rt_epsRa = dset["rt_eps"]
wt_epsRa = dset["wt_eps"]



using LaTeXStrings
p1= plot(Zt_epsT ,label="Truncation \$N=5\$", linewidth = 2, dpi=300, color = :blue)
p1= plot!(Zt_epsT2 ,label="Truncation \$N=2\$", linewidth = 2, dpi=300, color = :red)
p1= plot!(Zt_epsRa ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
title!("\$Z\$")

p2= plot(GDPt_epsT ,label="Truncation \$N=5\$", linewidth = 2, dpi=300, color = :blue)
p2= plot!(GDPt_epsT2 ,label="Truncation \$N=2\$", linewidth = 2, dpi=300, color = :red)
p2= plot!(GDPt_epsRa ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
title!("\$GDP\$")

p3= plot(Ct_epsT ,label="Truncation \$N=5\$", linewidth = 2, dpi=300, color = :blue)
p3= plot!(Ct_epsT2 ,label="Truncation \$N=2\$", linewidth = 2, dpi=300, color = :red)
p3= plot!(Ct_epsRa ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
title!("\$C\$")

p4= plot(Cbt_epsT ,label="Truncation \$N=5\$", linewidth = 2, dpi=300, color = :blue)
p4= plot!(Cbt_epsT2 ,label="Truncation \$N=2\$", linewidth = 2, dpi=300, color = :red)
p4= plot!(Ct_epsRa ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
title!("\$C_{bot}\$")


p7= plot(Kt_epsT ,label="Truncation \$N=5\$", linewidth = 2, dpi=300, color = :blue)
p7= plot!(Kt_epsT2 ,label="Truncation \$N=2\$", linewidth = 2, dpi=300, color = :red)
p7= plot!(Kt_epsRa ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
title!("\$K\$")

p8= plot(It_epsT ,label="Truncation \$N=5\$", linewidth = 2, dpi=300, color = :blue)
p8= plot!(It_epsT2 ,label="Truncation \$N=2\$", linewidth = 2, dpi=300, color = :red)
p8= plot!(It_epsRa ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
title!("\$I\$")

p9= plot(Wt_epsT ,label="Truncation \$N=5\$", linewidth = 2, dpi=300, color = :blue)
p9= plot!(Wt_epsT2 ,label="Truncation \$N=2\$", linewidth = 2, dpi=300, color = :red)
p9= plot!(Wt_epsRa ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
title!("\$Welfare\$")

p10= plot(Rt_epsT ,label="Truncation \$N=5\$", linewidth = 2, dpi=300, color = :blue)
p10= plot!(Rt_epsT2 ,label="Truncation \$N=2\$", linewidth = 2, dpi=300, color = :red)
p10= plot!(Rt_epsRa ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
title!("\$R\$")

p11= plot(wt_epsT ,label="Truncation \$N=5\$", linewidth = 2, dpi=300, color = :blue)
p11= plot!(wt_epsT2 ,label="Truncation \$N=2\$", linewidth = 2, dpi=300, color = :red)
p11= plot!(wt_epsRa ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
title!("\$w\$")

plot(p2, p3, p4,p7, p8, p9, p10, p11, layout=(4,2), size=(750, 750), dpi=300)
savefig("IRFcomp")
