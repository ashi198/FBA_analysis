function[model_with_formulas, formulas] = createMetaboliteFormulas(old_model, compounds)

%% loop to check for compounds from ModelSeed and Kbase files and extract formulas
% mets have the suffix "[c0], [e0], [b]" at the end of every metabolite ID. The
% following code is to remove this 

% Convert the structure to a cell array
dataCell = struct2cell(old_model);

% Get the index of the "mets" field
metsIndex = find(strcmp(fieldnames(old_model), 'mets'));
dataStruct= [];

% Remove the suffix "[c0], [e0], [b]" from all entries in the "mets" field
dataCell{metsIndex} = cellfun(@(x) strrep(x, '[c0]', ''), dataCell{metsIndex}, 'UniformOutput', false);
dataCell{metsIndex} = cellfun(@(x) strrep(x, '[e0]', ''), dataCell{metsIndex}, 'UniformOutput', false);
dataCell{metsIndex} = cellfun(@(x) strrep(x, '[b]', ''), dataCell{metsIndex}, 'UniformOutput', false);

% Convert the cell array back to a structure
dataStruct = cell2struct(dataCell, fieldnames(old_model), 1);

%% Define the metabolite IDs to search from compound database
%extract query metabolite IDs
met_ME= dataStruct.mets;

%extract metabolite IDs from the database 
compound_id= compounds.id;

%extract fromulas from the database 
formula_com= compounds.formula;

%extract indexs of memebers
[~, idx] = ismember(met_ME, compound_id);

% Extract the formulas for the IDs that were found
formulas = cell(size(idx));
for i = 1:numel(idx)
    if idx(i) > 0 && idx(i) <= numel(formula_com)
        formulas{i} = formula_com{idx(i)};
    end
end
% 
% add the extracted formulas in Model Seed
old_model.metFormulas=formulas;
model_with_formulas=old_model;