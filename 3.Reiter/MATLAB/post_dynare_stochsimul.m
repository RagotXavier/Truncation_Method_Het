tmp = oo_.irfs;
save('irfs_Reiter.mat','tmp');
                                                                                                            
fprintf('\n')

var_list = {"GDP","C","K", "W"};
means = containers.Map;
for k = 1:length(M_.endo_names)
    var = M_.endo_names{k};
    if ismember(var, var_list);
        means(var) = oo_.steady_state(k);
    end
end


fprintf("%15s %15s %15s\n", "Variable", "Mean", "Std-dev/mean");
for k = 1:length(oo_.var_list)
    var = oo_.var_list{k};
    if ismember(var, var_list);
        fprintf("%15s %15.3f %15.3f %15.3f\n", var, means(var), 100*oo_.var(k,k)^0.5/oo_.mean(k), 100*oo_.var(k,k)^0.5/means(var));
    end
end

fprintf('\n\n')

%auto correlations
var_list = {"GDP","C"};
indices = containers.Map;

for k = 1:length(oo_.var_list)
    var = oo_.var_list{k};
    if ismember(var, var_list);
        indices(var) = k; 
    end
end

fprintf("%15s %15s\n", "Variable","Autocorrel");
for i=1:length(var_list)
    k = indices(var_list{i});
    fprintf("%15s %15.3f\n", oo_.var_list{k}, 100*oo_.autocorr{1}(k,k));
end

fprintf('\n\n')
fprintf("%15s %15s %15s\n", "Variable 1","Variable 2","Correlation");
for i=1:length(var_list)    
    for j=1:length(var_list)
        if j<=i
            continue
        end
        k1 = indices(var_list{i});
        k2 = indices(var_list{j});
        fprintf("%15s %15s %15.3f\n", oo_.var_list{k1}, oo_.var_list{k2}, 100*oo_.var(k1,k2)/(oo_.var(k1,k1)*oo_.var(k2,k2))^0.5);
    end
end

fprintf('\n')
