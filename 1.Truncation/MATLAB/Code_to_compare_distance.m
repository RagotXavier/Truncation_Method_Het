%must be in 1.Reiter/MATLAB

load result_dynare_N6.mat
G6 = oo_.endo_simul(14,:);
Y6 = oo_.endo_simul(11,:);

save case6 G6 Y6

load result_dynare_N5.mat
G5 = oo_.endo_simul(14,:);
Y5 = oo_.endo_simul(11,:);

save case5 G5 Y5


load result_dynare_N4.mat
G4 = oo_.endo_simul(14,:);
Y4 = oo_.endo_simul(11,:);

save case4 G4 Y4

load result_dynare_N3.mat
G3 = oo_.endo_simul(14,:);
Y3 = oo_.endo_simul(11,:);

save case3 G3 Y3

load result_dynare_N2.mat
G2 = oo_.endo_simul(14,:);
Y2 = oo_.endo_simul(11,:);

save case2 G2 Y2

clear 

load case2
load case3
load case4
load case5
load case6

dist5=max(abs( (G5-G6)./G6));
dist4=max(abs( (G4-G6)./G6));
dist3=max(abs( (G3-G6)./G6));
dist2=max(abs( (G2-G6)./G6));



% load ../../3.Reiter/MATLAB/irfs_Reiter.mat
% 
% %
% figure; 
% plot(GDP_eps)
% hold on
% plot(tmp.GDP_eps)
% %hold off;
% 
% distanceG = (GDP_eps'-tmp.GDP_eps);
% plot(distanceG)
% dG = max(abs(distanceG))
% 
% figure;
% plot(C_eps)
% hold on
% plot(tmp.C_eps)
% %hold off;
% 
% distanceC = (C_eps'-tmp.C_eps);
% plot(distanceC)
% dC = max(abs(distanceC))
% 
% figure;
% plot(TT_eps)
% hold on
% plot(tmp.TT_eps)
% %hold off;
% 
% distanceTT = (TT_eps'-tmp.TT_eps);
% plot(distanceTT)
% dTT = max(abs(distanceTT))