%% Retrieve the indices and names of metabolites that have a name similar to the metabolite of interest

function[metabolites, metabolite_names, strict_metabolites, strict_metabolite_names]= findMetsFromList(models, metabolites_of_interest)

%models need to be a cell array with all the models stored in each cell 
metabolites = cell(length(models), 1);
metabolite_names = cell(length(models), 1);
strict_metabolites= cell(length(models), 1);
strict_metabolite_names= cell(length(models), 1);
metabolites= cell(length(models), 1);

for i = 1:length(models)
    read_mets = models{i}.metNames; % finds the names of all metabolites from metNames for each model
    mets = models{i}.mets;% find of all metabolites from mets from each model
    mets_idx = cell(length(metabolites_of_interest), 1);
   for l=1:length(read_mets)
    read_mets{l} = regexprep(read_mets{l}, '_c0$', '');
    read_mets{l} = regexprep(read_mets{l}, '_e0$', '');
   end 
   for j=1:length(metabolites_of_interest)
    mets_idxs{j, i}= ismember(lower(read_mets), lower(metabolites_of_interest{j}));
    mets_idx{j, i}= find(contains(lower(read_mets), lower(metabolites_of_interest{j}), "IgnoreCase",true)); % find the index of the metabolite of interest
    if ~isempty(mets_idx{j, i})
        metabolites{j, i} = mets(mets_idx{j, i});
        metabolite_names{j, i} = read_mets(mets_idx{j, i});
    end
    if ~isempty(mets_idxs{j, i})
        strict_metabolites{j, i}=mets(mets_idxs{j, i});
        strict_metabolite_names{j, i} = read_mets(mets_idxs{j, i});
    end 
   end 
end