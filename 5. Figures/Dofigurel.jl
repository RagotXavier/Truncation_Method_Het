 ###########READ MATLAB FILES FOR THE IRF AND PLOT THE IRF COMPARISON
using MAT
using Plots; pyplot();



TimePeriods = 200 #time periods for the simulation



fileIn = matopen("../2. BKM/tofigBKM.mat")
dset = read(fileIn)
Z_eps = dset["Z_eps"]
GDP_eps= dset["GDP_eps"]
C_eps = dset["C_eps"]
Cbs_eps = dset["Cbs_eps"]
K_eps= dset["K_eps"]
I_eps= dset["I_eps"]
W_eps= dset["W_eps"] # be careful : welfare is in in level variation
R_eps= dset["R_eps"]
w_eps= dset["w_eps"]
Zt= dset["Zt"]
GDPt= dset["GDPt"]
Consot= dset["Consot"]
Capitalt= dset["Capitalt"]
Investmentt= dset["Investmentt"]
Welfaret= dset["Welfaret"]
Rsimul= dset["Rsimul"]
wsimul= dset["wsimul"]
close(fileIn)



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
Kt_epsT= dset["Kt_eps"]
It_epsT= dset["It_eps"]
Wt_epsT= dset["Wt_eps"]
Rt_epsT= dset["rt_eps"]
wt_epsT= dset["wt_eps"]


fileIn = matopen("../3. Reiter/MATLAB_DYNAMICS/tofigreiter.mat")
dset = read(fileIn)
close(fileIn)
mut_epsR = dset["ut_eps"]
Zt_epsR = dset["Zt_eps"]
GDPt_epsR= dset["GDPt_eps"]
Ct_epsR = dset["Ct_eps"]
Cbt_epsR = dset["Cbt_eps"]
Kt_epsR= dset["Kt_eps"]
It_epsR= dset["It_eps"]
Wt_epsR= dset["Wt_eps"]
Rt_epsR= dset["rt_eps"]
wt_epsR= dset["wt_eps"]

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
p1= plot!(Zt_epsR ,label="Reiter", linewidth = 2, dpi=300,linestyle =:dash)
p1= plot!(Z_eps[1:TimePeriods], label="BKM", linewidth = 2, dpi=300, color = :green, linestyle =:dash)
p1= plot!(Zt_epsRa ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
title!("\$Z\$")

p2= plot(GDPt_epsT ,label="Truncation \$N=5\$", linewidth = 2, dpi=300, color = :blue)
p2= plot!(GDPt_epsR ,label="Reiter", linewidth = 2, dpi=300,linestyle =:dash)
p2= plot!(GDP_eps[1:TimePeriods], label="BKM", linewidth = 2, dpi=300, color = :green, linestyle =:dash)
p2= plot!(GDPt_epsRa ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
title!("\$GDP\$")

p3= plot(Ct_epsT ,label="Truncation \$N=5\$", linewidth = 2, dpi=300, color = :blue)
p3= plot!(Ct_epsR ,label="Reiter", linewidth = 2, dpi=300,linestyle =:dash)
p3= plot!(C_eps[1:TimePeriods], label="BKM", linewidth = 2, dpi=300, color = :green, linestyle =:dash)
p3= plot!(Ct_epsRa ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
title!("\$C\$")

p4= plot(Cbt_epsT ,label="Truncation \$N=5\$", linewidth = 2, dpi=300, color = :blue)
p4= plot!(Cbt_epsR ,label="Reiter", linewidth = 2, dpi=300,linestyle =:dash)
p4= plot!(Cbs_eps[1:TimePeriods], label="BKM", linewidth = 2, dpi=300, color = :green, linestyle =:dash)
p4= plot!(Ct_epsRa ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
title!("\$C_{bot}\$")

p7= plot(Kt_epsT ,label="Truncation \$N=5\$", linewidth = 2, dpi=300, color = :blue)
p7= plot!(Kt_epsR ,label="Reiter", linewidth = 2, dpi=300,linestyle =:dash)
p7= plot!(K_eps[1:TimePeriods], label="BKM", linewidth = 2, dpi=300, color = :green, linestyle =:dash)
p7= plot!(Kt_epsRa ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
title!("\$K\$")

p8= plot(It_epsT ,label="Truncation \$N=5\$", linewidth = 2, dpi=300, color = :blue)
p8= plot!(It_epsR ,label="Reiter", linewidth = 2, dpi=300,linestyle =:dash)
p8= plot!(I_eps[1:TimePeriods], label="BKM", linewidth = 2, dpi=300, color = :green, linestyle =:dash)
p8= plot!(It_epsRa ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
title!("\$I\$")

p9= plot(Wt_epsT ,label="Truncation \$N=5\$", linewidth = 2, dpi=300, color = :blue)
p9= plot!(Wt_epsR ,label="Reiter", linewidth = 2, dpi=300,linestyle =:dash)
p9= plot!(W_eps[1:TimePeriods], label="BKM", linewidth = 2, dpi=300, color = :green, linestyle =:dash)
p9= plot!(Wt_epsRa ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
title!("\$Welfare\$")

p10= plot(Rt_epsT ,label="Truncation \$N=5\$", linewidth = 2, dpi=300, color = :blue)
p10= plot!(Rt_epsR ,label="Reiter", linewidth = 2, dpi=300,linestyle =:dash)
p10= plot!(R_eps[1:TimePeriods], label="BKM", linewidth = 2, dpi=300, color = :green, linestyle =:dash)
p10= plot!(Rt_epsRa ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
title!("\$R\$")

p11= plot(wt_epsT ,label="Truncation \$N=5\$", linewidth = 2, dpi=300, color = :blue)
p11= plot!(wt_epsR ,label="Reiter", linewidth = 2, dpi=300,linestyle =:dash)
p11= plot!(w_eps[1:TimePeriods], label="BKM", linewidth = 2, dpi=300, color = :green, linestyle =:dash)
p11= plot!(wt_epsRa ,label="RA", linewidth = 2, dpi=300,linestyle =:dash)
title!("\$w\$")

plot(p2, p3, p4,p7, p8, p9, p10, p11, layout=(4,2), size=(750, 750), dpi=300)
savefig("IRF")


#############READ MATLAB FILES FOR THE TIME SERIES AND PLOT THE TIME SERIES COMPARISON
using MAT


fileIn = matopen("../1. Truncation/MATLAB_DYNAMICS/tofigtruncationsimt.mat")
#fileIn = matopen("tofigtruncationsimt.mat")
dset = read(fileIn)
close(fileIn)
mutT = dset["ut"]
ZtT = dset["Zt"]
GDPtT= dset["GDPt"]
CtT = dset["Ct"]
CbtT = dset["Cbt"]
KtT= dset["Kt"]
ItT= dset["It"]
WtT= dset["Wt"]
RtT= dset["rt"]
wtT= dset["wt"]

using MAT
fileIn = matopen("../3. Reiter/MATLAB_DYNAMICS/tofigreitersimt.mat")
dset = read(fileIn)
close(fileIn)
mutR = dset["ut"]
ZtR = dset["Zt"]
GDPtR= dset["GDPt"]
CtR = dset["Ct"]
CbtR = dset["Cbt"]
KtR= dset["Kt"]
ItR= dset["It"]
WtR= dset["Wt"]
RtR= dset["rt"]
wtR= dset["wt"]

using LaTeXStrings
p1= plot(ZtT[TimePeriods:end] ,label="Truncation \$N=5\$", linewidth = 0.2, dpi=300, color = :blue)
p1= plot!(ZtR[TimePeriods:end] ,label="Reiter", linewidth = 0.2, dpi=300,linestyle =:dash)
p1= plot!(Zt, label="BKM", linewidth = 0.2, dpi=300, color = :green, linestyle =:dash)
title!("\$Z\$")

p2= plot(GDPtT[TimePeriods:end] ,label="Truncation \$N=5\$", linewidth = 0.2, dpi=300, color = :blue)
p2= plot!(GDPtR[TimePeriods:end] ,label="Reiter", linewidth = 0.2, dpi=300,linestyle =:dash)
p2= plot!(GDPt, label="BKM", linewidth = 0.2, dpi=300, color = :green, linestyle =:dash)
title!("\$GDP\$")

p3= plot(CtT[TimePeriods:end] ,label="Truncation \$N=5\$", linewidth = 0.2, dpi=300, color = :blue)
p3= plot!(CtR[TimePeriods:end] ,label="Reiter", linewidth = 0.2, dpi=300,linestyle =:dash)
p3= plot!(Consot, label="BKM", linewidth = 0.2, dpi=300, color = :green, linestyle =:dash)
title!("\$C\$")

p4= plot(KtT[TimePeriods:end] ,label="Truncation \$N=5\$", linewidth = 0.2, dpi=300, color = :blue)
p4= plot!(KtR[TimePeriods:end] ,label="Reiter", linewidth = 0.2, dpi=300,linestyle =:dash)
p4= plot!(Capitalt, label="BKM", linewidth = 0.2, dpi=300, color = :green, linestyle =:dash)
title!("\$K\$")

p5= plot(ItT[TimePeriods:end] ,label="Truncation", linewidth = 0.2, dpi=300, color = :blue)
p5= plot!(ItR[TimePeriods:end] ,label="Reiter", linewidth = 0.2, dpi=300,linestyle =:dash)
p5= plot!(Investmentt, label="BKM", linewidth = 0.2, dpi=300, color = :green, linestyle =:dash)
title!("\$I\$")

p6= plot(WtT[TimePeriods:end] ,label="Truncation", linewidth = 0.2, dpi=300, color = :blue)
p6= plot!(WtR[TimePeriods:end] ,label="Reiter", linewidth = 0.2, dpi=300,linestyle =:dash)
p6= plot!(Welfaret, label="BKM", linewidth = 0.2, dpi=300, color = :green, linestyle =:dash)
title!("\$Welfare\$")

p7= plot(RtT[TimePeriods:end] ,label="Truncation", linewidth = 0.2, dpi=300, color = :blue)
p7= plot!(RtR[TimePeriods:end] ,label="Reiter", linewidth = 0.2, dpi=300,linestyle =:dash)
p7= plot!(100*Rsimul, label="BKM", linewidth = 0.2, dpi=300, color = :green, linestyle =:dash)
title!("\$R\$")

p8= plot(wtT[TimePeriods:end] ,label="Truncation", linewidth = 0.2, dpi=300, color = :blue)
p8= plot!(wtR[TimePeriods:end] ,label="Reiter", linewidth = 0.2, dpi=300,linestyle =:dash)
p8= plot!(100*wsimul, label="BKM", linewidth = 0.2, dpi=300, color = :green, linestyle =:dash)
title!("\$w\$")



display(plot(p1, p2, p3, p4, p5, p6,p7,p8,layout=(4,2), size=(750, 750), dpi=300))
savefig("Time_series")

