function [updated_rxns, indices, exchange]= findUptakeRxns(models)
% Initialize the cell arrays for exchange reactions and uptakes
exc_rxns = cell(1, length(models));
up_takes = cell(1, length(models));
updated_rxns = cell(1, length(models));

% Loop through the models and find exchange reactions and uptakes
for i = 1:length(models)
    [exc_rxns{i}, up_takes{i}] = findExcRxns(models{i});
    nonzero_indices = find(up_takes{i});
    nonzero_indices_exchange = find(exc_rxns{i});
    
    % Initialize a temporary cell array to store uptake reactions
    temp_updated_rxns = cell(length(nonzero_indices), 1);
    temp_exchange_rxns = cell(length(nonzero_indices_exchange), 1);
    
    % Store the uptake reactions for the current model
    for j = 1:length(nonzero_indices)
        temp_updated_rxns{j} = models{i}.rxns(nonzero_indices(j));
    end

    for k=1:length(nonzero_indices_exchange)
        temp_exchange_rxns{k}=models{i}.rxns(nonzero_indices_exchange(k)); 
    end 
    % Store the temporary cell array in the updated_rxns cell array
    updated_rxns{i} = temp_updated_rxns;
    indices{i} = nonzero_indices;
    exchange{i} = temp_exchange_rxns;
end
end 