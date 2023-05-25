% For fillGap: only use importModel for this 
initCobraToolbox(false)
solverOK=changeCobraSolver('gurobi','LP');
%% Define names of all models 
model_names= {'1A01_carveMe.xml'; 'Demeter_1A01_x.xml'; 'Raven_edited_main.xml';'PK_iJN1463.xml'; 'SE_iJB785.xml'};

for i=1:length(model_names)
    model{i}= importModel(model_names{i});
end

%% extra stuff for updating fields 
rev=createRev(model{1});
model{1}.rev=rev;
model{1}= readCbModel(model_names{1});

%add metabolite formula before running this 
% Add metabolite formulas 
compounds = readtable('output.csv');
for i=1:2 % only for modelSeed/kbase 
% use createMetaboliteFormulas to add metabolite formulas 
[model{i}, ~] = createMetaboliteFormulas(model{i}, compounds);
end 

%% find blocked reactions and deadend metabolites
for i=1:3
deadEnds{i}=detectDeadEnds(model{i});
blocked{i}=findBlockedReaction(model{i});
end 
%% 
%Try to fill gaps using the full KEGG/ targeted models to see if that gives a
%significantly higher number

for k=4
keggModels{k}=importModel(model_names{k});
end 
%use this for Kegg reference model 
%keggModel=getModelFromKEGG([],false,false,false,false);

%Remove genes as they will not be used for the gapfilling, so they are removed to make this a little
%faster
for w=1:length(keggModels)
keggModels{w}=rmfield(keggModels{w},'genes');
keggModels{w}=rmfield(keggModels{w},'rxnGeneMat');
balanceStructure=getElementalBalance(keggModels{w});
keggModels{w}=removeReactions(keggModels{w},balanceStructure.balanceStatus~=1,true,true);
end 

%The function fillGaps with these settings will try to include reactions in
%order to have flux through all reactions in the model. There are other
%settings as well. The first flag says that production of all metabolites
%should be allowed.
%check if KBase/ModelSeed contain genes, ids, and grRules

params.relGap=0.6; %Lower number for a more exhaustive search
params.printReport=true;

for i=1:3
[newConnected_new{i}, cannotConnect_new{i}, addedRxns_new{i}, newModel]=fillGaps(model{i},keggModels{1},true,false,false,[],params);
newModels{i}=newModel;
end  

%Continue to improve the connectivity of the model by identifying
%metabolites that should be connected. A convenient way to get an overview
%of how connected the model is, and at the same time getting a lot of
%useful data, is to use gapReport. Note that it can take several to many
%hours to run, depending on the number of gaps in the model.
[noFluxRxns, noFluxRxnsRelaxed, subGraphs, notProducedMets, minToConnect,...
    neededForProductionMat]=gapReport(newModels{4});

% ensure that the new model grows 
model{1} = changeObjective(model{1}, 'Growth'); % Set objective to biomass
sol = optimizeCbModel(model{1}, 'max');

newModels{1} = changeObjective(newModels{1}, 'Growth'); % Set objective to biomass
sol = optimizeCbModel(newModels{1}, 'max');