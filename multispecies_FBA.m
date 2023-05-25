% initialize cobra toolbox 
initCobraToolbox('false')
solverOK=changeCobraSolver('gurobi','LP');
%% Upload Kbase and CarveMe models 
model_names = cell(2, 1);
models_names={'1A01_Kbase.sbml'; '3B05_Kbase.sbml'; 'B2R09_Kbase.sbml'; '5D01_Kbase.sbml'; 'AS40_Kbase.sbml'; ...
    'I3M17_Kbase.sbml'; '6B07_Kbase.sbml'; 'C1R06_Kbase.sbml'};

for i=1:length(models_names)
    models{i}= readCbModel(models_names{i})
end 
models_original=models;

%% Find the index of the metabolite of interest
%Dulcitol=Galactitol
metabolites_of_interest={'Acetate'; 'glucose';'fructose';'glycerol'; 'galactose'; 'glucosamine'; 'alanine'; 'mannose'; 'Glycine'; 'Histidine'; 'Citrate'; 'Arginine';'Glutamine'; 'Lactate'; 'proline'; 'Lactose'; 'maltose'; 'mannitol'; 'succinate'; 'taurine'; 'cellobiose'; 'ribose'; 'xylose';....
    'Arabinose'; 'Sucrose'; 'Urea'; 'Valine'; 'Acetaldehyde'; 'Ethanol'; 'Adenine'; 'Sorbitol'; 'melibiose'; 'threonine'; 'Oxaloacetate'; 'Propionate'; 'Pyruvate'; 'Beta-alanine'; 'cysteine'; 'Galactitol'; 'lysine'; 'Leucine'; 'Methionine';....
    'tyrosine'; 'methanol'; 'O2'};

[metabolites, metabolite_names, strict_metabolites, strict_metabolite_names]= findMetsFromList(models, metabolites_of_interest);
metabolites_for_uptake=strict_metabolites;

%% findUptakeRxns will produce a list of the uptake reactions present in all models  
[updated_rxns]= findUptakeRxns(models); 

%% find out all uptake reactions that are associated with only the metabolites of interest 
matching_uptake_reactions=cell(length(updated_rxns), 1);
meta_temp=cell(length(metabolites_for_uptake), 1); 

for i=1:length(models)
  % have all metabolites for one models sequentially
  meta_temp=cell(length(metabolites_for_uptake), 1); 
  for s=1:length(metabolites_for_uptake)
  meta_temp{s, i}=metabolites_for_uptake{s, i}; 
  end  

react_temp=updated_rxns{i};
% remove extra characters from the metabolite name
   for l=1:length(meta_temp)
    meta_temp = meta_temp(~cellfun('isempty', meta_temp));
    meta_temp{l} = regexprep(meta_temp{l}, '\[c0\]$', '');
    meta_temp{l} = regexprep(meta_temp{l}, '\[e0\]$', '');
    meta_temp{l} = regexprep(meta_temp{l}, '\[C_c\]$', '');
    meta_temp{l} = regexprep(meta_temp{l}, '\[C_e\]$', ''); 
    meta_temp{l} = regexprep(meta_temp{l}, '\[c\]$', '');
    meta_temp{l} = regexprep(meta_temp{l}, '\[e\]$', '');
    meta_temp{l} = regexprep(meta_temp{l}, '\[p\]$', '');
    meta_temp{l} = regexprep(meta_temp{l}, '\[C_p\]$', '');
   end
% pick all the uptake reactions associated with desired metabolites
   for j=1:length(meta_temp)
       idx=meta_temp{j};
       matching_reactions= {}; 
     for k=1:length(react_temp)
       if contains(lower(react_temp{k}), lower(idx))
            matching_reactions{end+1}= react_temp{k};
       end
     end
    matching_uptake_reactions{j, i} = matching_reactions;
   end
end

%% Perform FBA
models=models_original;
biomass_yield=cell(length(matching_uptake_reactions), length(models)); 
biomass_reactions={'bio1';'bio1';'bio1';'bio1';'bio1';'bio1'; 'bio1'; 'bio1'};

   for i=1:length(models)
    for j = 1:length(matching_uptake_reactions)
        %models{i}=models_original{i};
        if isempty(matching_uptake_reactions{j, i}) || isempty(biomass_reactions{i})
        biomass_yield{j, i}= 'NA';
        else 
        models{i} = changeRxnBounds(models{i}, matching_uptake_reactions{j, i}, -5, 'l'); % Set lower bound to -10
        models{i}=changeRxnBounds(models{i}, matching_uptake_reactions{45, i}, -1000, 'l'); %set to index of oxygen 
        models{i} = changeObjective(models{i}, biomass_reactions{i}); % Set objective to biomass
        sol = optimizeCbModel(models{i}, 'max'); % Perform FBA
        fprintf('Biomass yield on %s: %.6f\n', matching_uptake_reactions{j, i}, sol.f);
        biomass_yield{j, i}=sol.f;
        end 
    end
   end

   %% graphing heatplot
modelNames = {'1A01_Kbase.sbml'; '3B05_Kbase.sbml'; 'B2R09_Kbase.sbml'; '5D01_Kbase.sbml'; 'AS40_Kbase.sbml'; ...
    'I3M17_Kbase.sbml'; '6B07_Kbase.sbml'; 'C1R06_Kbase.sbml'};

temp=zeros(45, 8);
for i = 1:length(models)
    for j = 1:length(biomass_yield)
        if isempty(biomass_yield{j, i}) || strcmp(biomass_yield{j, i}, 'NA')  % Check for empty or 'NA' values
            temp(j, i) = 0;
        else
            temp(j, i) = biomass_yield{j, i};
        end
    end
end

h = heatmap(modelNames,metabolites_of_interest,temp);
h.Title = 'Biomass yields on substrates';
h.XLabel = 'Organisms';
h.YLabel = 'Substrates';


