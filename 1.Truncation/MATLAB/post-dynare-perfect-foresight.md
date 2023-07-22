This file has to be run _after_ Dynare simulations have been computed. It is automatically launched in [`Dynare_Truncation.m`](./Dynare_Truncation.html). The file is not meant to be launched otherwise.

The file saves the tax-to-GDP data (called $\tau$ in the paper) in the mat-file `taus.mat` for later plotting.

```matlab
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
```
