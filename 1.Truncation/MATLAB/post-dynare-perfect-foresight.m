fprintf('\n')


var_list = {"GDP","TT"};
indices = containers.Map;

for k = 1:length(M_.endo_names)
    var = M_.endo_names{k};
    if ismember(var, var_list);
        indices(var) = k; 
    end
end
taus = oo_.endo_simul(indices('TT'),:)./oo_.endo_simul(indices('GDP'),:);

save('taus.mat','taus')
