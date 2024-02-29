%% Information
% Author: Diego Rodriguez Esperante
% Contact: d.rodriguezesperante@student.maastrichtuniversity.nl
% Last update: 28/01/2024

%clear; clc; close all;

%% Initialize
initCobraToolbox(false);

solverName = 'glpk';
solverType = 'LP';
changeCobraSolver(solverName, solverType);

clear solverName solverType
%% Read files
% Import Merged model
modelFileName = ['Model files' filesep 'H1R3MergedModel.mat'];
modelFileName= [pwd filesep modelFileName];
mergedModel = readCbModel(modelFileName);

% Import our mtbs of interest
mtbs_interest = readtable(['CSV' filesep 'Reactions - Metabolites.csv']);
rxns_interest = readtable(['CSV' filesep 'Reactions - Rxn-Sub Pairs.csv']);

clear modelFileName

%% Find missing, unnamed, extra metabolites
merged_mets = extract(mergedModel.mets, lettersPattern(3) + digitsPattern(5));

% Missing (found in H1 and R3 database but not extracted correctly in COBRA)
[mtbs.missing.names, mtbs.missing.ids] = setdiff(mtbs_interest{:, 'Human1'}, merged_mets);
mtbs.missing.ids(strcmp(mtbs.missing.names, '')) = [];
mtbs.missing.names(strcmp(mtbs.missing.names, '')) = [];

% Unnamed (not in either H1 or R3)
mtbs.unnamed.ids = find((strcmp(mtbs_interest{:, 'Human1'}, '').*strcmp(mtbs_interest{:, 'Recon3D'}, ''))); % Find metabolites with no names in H1 or R3
mtbs.unnamed.names = mtbs_interest{mtbs.unnamed.ids, 'NameInDiagram'};

% Extracted from H1 and R3 but not of interest
[mtbs.extra.names, mtbs.extra.ids] = setdiff(merged_mets, mtbs_interest{:, 'Human1'});

%% Find missing, extra rxns
% Missing
[rxns.missing.names, rxns.missing.ids] = setdiff(rxns_interest{:,"Database"}, mergedModel.rxns);
rxns.missing.ids(strcmp(rxns.missing.names, '')) = [];
rxns.missing.names(strcmp(rxns.missing.names, '')) = [];

% Extra
[rxns.extra.names, rxns.extra.ids] = setdiff(mergedModel.rxns, rxns_interest{:,"Database"});

%% Add metabolite information
completeModel = mergedModel;
% Missing metabolites
completeModel = addMultipleMetabolites(completeModel, strcat(mtbs.missing.names, 'c'), 'metNames', mtbs_interest{mtbs.missing.ids, "NameInDiagram"});

% Unnamed metabolites
completeModel = addMultipleMetabolites(completeModel, mtbs.unnamed.names, 'metNames', mtbs_interest{mtbs.unnamed.ids, "NameInDiagram"});

%% Add reactions
rxn_list = rxns_interest{strcmp(rxns_interest{:,'Redundancy'}, 'No redundancy'), ["Name", "Enzyme", "Substrate", "Product", "Direction"]};
rxn_list = cell2table(rxn_list, "VariableNames", ["Name", "Enzyme", "Substrate", "Product", "Direction"]);

% Design stoichiometric matrix
mtb_list = unique([rxn_list{:, 'Substrate'}; rxn_list{:, 'Product'}]);
mtb_list = [mtb_list, strings(length(mtb_list), 1)];
for i = 1:length(mtb_list)
    id = find(strcmp(mtb_list(i, 1), completeModel.metNames), 1);
    mtb_list(i, 2) = completeModel.mets(id);
end

S_mat = zeros(length(mtb_list), height(rxn_list));
for i = 1:length(mtb_list)
    S_mat(i,:) = S_mat(i,:)' - strcmp(mtb_list{i, 1}, rxn_list{:,'Substrate'});
    S_mat(i,:) = S_mat(i,:)' + strcmp(mtb_list{i, 1}, rxn_list{:,'Product'});
end

% Add reactions to model, add grRules and genes
completeModel = addMultipleReactions(completeModel, rxn_list{:,"Name"}, mtb_list(:, 2), S_mat);

% Check for duplicates rxns, remove unused genes
completeModel = checkDuplicateRxn(completeModel, 'S', 1);
completeModel = removeUnusedGenes(completeModel);

%% Save model to file
% save(['Model files' filesep 'completeModel.mat'], 'completeModel');