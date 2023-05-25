%% Flux balance analysis workflow 
% in order to perform FBA, uptake reactions corresponding to the metabolite
% of interests are required. The following workflow finds that. 
%% 
% initialize cobra toolbox 
initCobraToolbox(false)
solverOK=changeCobraSolver('gurobi','LP');

%% upload models [Model retrieved from Raven should be kept on the 5th position]
model_names= {'1A01_Kbase.sbml'; 'ModelSeed_1A01.sbml';....
    '1A01_carveMe.xml'; 'Demeter_1A01_x.xml'; 'Raven_edited_main.xml'};
for i=1:length(model_names)
    models{i, 1}= readCbModel(model_names{i});
end
models_original= models;
models_original{5, 1}= importModel(model_names{5});

%% Find the index of the metabolite of interest
%Dulcitol=Galactitol
metabolites_of_interest={'Acetate'; 'glucose';'fructose';'glycerol'; 'galactose'; 'glucosamine'; 'alanine'; 'mannose'; 'Glycine'; 'Histidine'; 'Citrate'; 'Arginine';'Glutamine'; 'Lactate'; 'proline'; 'Lactose'; 'maltose'; 'mannitol'; 'succinate'; 'taurine'; 'cellobiose'; 'ribose'; 'xylose';....
    'Arabinose'; 'Sucrose'; 'Urea'; 'Valine'; 'Acetaldehyde'; 'Ethanol'; 'Adenine'; 'Sorbitol'; 'melibiose'; 'threonine'; 'Oxaloacetate'; 'Propionate'; 'Pyruvate'; 'Beta-alanine'; 'cysteine'; 'Galactitol'; 'lysine'; 'Leucine'; 'Methionine';....
    'tyrosine'; 'methanol'; 'O2'};

% findMetsFromList retrieve the indices and names of metabolites that have
% a name similar to the metabolite of interest. Metabolites and
% metabolite_names will produce a list of all metabolite that contains the
% name of the metabolite z.b the metabolites file for "acetate" will give
% results such as "". Strict_metabolites and strict_metabolite_names
% produce a list of metabolites that matches exactly with the query
% metabolite. 

%metabolite_for_uptake should be updated such that only one entry per
%metabolite is present in every cell. 

[metabolites, metabolite_names, strict_metabolites, strict_metabolite_names]= findMetsFromList(models, metabolites_of_interest);
metabolites_for_uptake=strict_metabolites; %% select all the appropriate metabolites before proceeding further

%% findUptakeRxns will produce a list of the uptake reactions present in all models  
[updated_rxns]= findUptakeRxns(models); 

%% find out all uptake reactions that are associated with only the metabolites of interest 
matching_uptake_reactions=cell(length(updated_rxns), 1);
meta_temp=cell(length(metabolites_for_uptake), 1); 

%% warning--> 
%Raven models have a different format of mets. Please compare formulas of
%uptake reactions instead of names to find a match

mets_raven=cell(length(metabolites_for_uptake), 1);
mets_raven= metabolites_for_uptake{:, 5}; %assuming number 5 is the column for Raven 
for_raven_reactions=updated_rxns{5}; % to get all uptake reactions relevant to raven 

temp={};
for i=1:length(for_raven_reactions)
  temp=for_raven_reactions{i};
  only_for_raven_formulas{i, 1}=printRxnFormulaOri(models{5}, temp);
  only_for_raven_formulas{i, 2}= temp; 
end 

%% 
for i=1:length(models)
  
  meta_temp=cell(length(metabolites_for_uptake), 5); 
  % have all metabolites for one models sequentially
   for s=1:length(metabolites_for_uptake)
    meta_temp{s, i}=metabolites_for_uptake{s, i}; 
   end  
   
  % Condition for inputing raven with formulas instead of reaction names 
  if i ~=5
    react_temp=updated_rxns{i}; % all uptake reactions associated with one model 
  else
     react_temp = cell(length(only_for_raven_formulas), 1);
     react_temp= only_for_raven_formulas; 
  end 
 
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

%% Flux balance analysis
models=models_original;
biomass_yield=cell(length(matching_uptake_reactions), length(models)); 
%set biomass reactions 
biomass_reactions={'bio1';'bio1';'Growth';'bio1'; 'growth'};
%% adding biomass reaction for Raven 
 models_original{5} = addReaction(models_original{5},'growth',{'G3P[c]','ACETYL-COA[c]','ATP[c]','ERYTHROSE-4P[c]','C00085[c]','GAP[c]','C00092[c]','GLN[c]', ...
    'GLT[c]','WATER[c]','NAD[c]','NADPH[c]','OXALACETIC_ACID[c]','PHOSPHO-ENOL-PYRUVATE[c]','PYRUVATE[c]','C00117[c]','ADP[c]','2-KETOGLUTARATE[c]','Coenzyme A','PROTON[c]','NADH[c]',' NADP[c]','Pi[c]'}, [-1.496 -3.7478 -59.81 -0.361 -1 -0.129 -0.205 -0.2557 -4.9414 -59.81 -3.547 -13.0279 -1.7867 -0.5191 -2.8328 -0.8977 59.81 ....
    4.1182 3.7478 59.81 3.547 13.0279 59.81],false);

%% Perform FBA 
   for i=1:length(models)
    for j = 1:(length(matching_uptake_reactions)-1)
        %models{i}=models_original{i};
        if matching_uptake_reactions{j, i} == '0'
        biomass_yield{j, i}= 'NA';
        else 
        models{i} = changeRxnBounds(models{i}, matching_uptake_reactions{j, i}, -10, 'l'); % Set lower bound to -10
        models{i}=changeRxnBounds(models{i}, matching_uptake_reactions{45, i}, -1000, 'l'); %set to index of oxygen 
        models{i} = changeObjective(models{i}, biomass_reactions{i}); % Set objective to biomass
        sol = optimizeCbModel(models{i}, 'max'); % Perform FBA
        fprintf('Biomass yield on %s: %.6f\n', matching_uptake_reactions{j, i}, sol.f);
        biomass_yield{j, i}=sol.f;
        end
    end
  end 
 %% map a heat plot 
names=model_names;
temp_1 = zeros(size(biomass_yield));  % Initialize temp_1 with zeros

for i = 1:length(models)
    for j = 1:length(biomass_yield)
        if isempty(biomass_yield{j, i}) || strcmp(biomass_yield{j, i}, 'NA')  % Check for empty or 'NA' values
            temp_1(j, i) = 0;
        else
            temp_1(j, i) = biomass_yield{j, i};
        end
    end
end
h = heatmap(names,metabolites_of_interest,temp_1);
h.Title = 'Biomass yields on substrates';
h.XLabel = 'Organisms';
h.YLabel = 'Substrates';

%% for acetate secretation to check for certain substrates (glucosamine, glucose, galactose, pyruvate)
sub_check={'[]' '[]' 'R_EX_gam_e' 'EX_gam(e)' '[]'; 'EX_cpd00027_e0' 'EX_cpd00027_e0' 'R_EX_glc__D_e' 'EX_glc_D(e)' '[]';....
    '[]' 'R_EX_gal_e' 'EX_galt(e)' '[]' '[]'; '[]' '[]' 'R_EX_pyr_e' 'EX_pyr(e)' '[]'};
%% find indexes associated with acetate

acetate={'cpd00029[e0]'	'cpd00029[c0]'	'ac[C_e]'	'ac[e]'	'ACET[c]'};
for j=1:length(acetate)
 temp_2{j}=acetate{j};
 acetateIndex(j) = find(strcmp(temp_2{j}, models{j}.mets));
end  

% check acetate secretation for the following carbon sources 

for i=1:length(models)
   for j=1:(length(sub_check)-1)
    if isempty(sub_check{j, i})
    acetate_flux{j, i}= 'NA';
    else 
    models{i}=models_original{i};
    models{i} = changeRxnBounds(models{i}, sub_check{j, i}, -10, 'l'); % Set lower bound to -10
    models{i}=changeRxnBounds(models{i}, matching_uptake_reactions{45, i}, -1000, 'l'); 
    models{i} = changeObjective(models{i}, biomass_reactions{i}); % Set objective to biomass
    sol = optimizeCbModel(models{i}, 'max');
    acetate_flux{j, i} = sol.v(acetateIndex(i));
    end
   end 
end
