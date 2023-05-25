function rev= createRev(old_model)
%%
n_rxns = length(old_model.rxns); % number of reactions in the model
rev = zeros(n_rxns, 1); % initialize the 'rev' vector
for i = 1:n_rxns
    if old_model.lb(i) == 0 && old_model.ub(i) > 0
        rev(i) = 1; % reaction is reversible
    elseif old_model.lb(i) < 0 && old_model.ub(i) == 0
        rev(i) = 1; % reaction is reversible
    elseif old_model.lb(i) < 0 && old_model.ub(i) > 0
        rev(i) = 0; % reaction is partially reversible
    end
end