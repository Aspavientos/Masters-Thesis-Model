%% Information
% Author: Diego Rodriguez Esperante
% Contact: d.rodriguezesperante@student.maastrichtuniversity.nl
% Last update: 21/01/2024

%clear; clc; close all;

%% Initialize
initCobraToolbox(false);

solverName = 'glpk';
solverType = 'LP';
changeCobraSolver(solverName, solverType);

clear solverName solverType
%% Read files
% Import Merged model
modelFileName = 'Model files/H1R3MergedModel.mat';
modelFileName= [pwd filesep modelFileName];
mergedModel = readCbModel(modelFileName);

% Import our mtbs of interest
mtbs_interest = readtable("CSV/Reactions - Metabolites.csv");
rxns_interest = readtable("CSV/Reactions - Rxn-Sub Pairs.csv");

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
rxn_list = rxns_interest{strcmp(rxns_interest{:,'Redundancy'}, 'No redundancy'), ["Enzyme", "Substrate", "Product", "Direction"]};
rxn_list = cell2table(rxn_list, "VariableNames", ["Enzyme", "Substrate", "Product", "Direction"]);
rxn_list = addvars(rxn_list, zeros(height(rxn_list), 1), zeros(height(rxn_list), 1), zeros(height(rxn_list), 1), 'NewVariableNames', {'SubID', 'ProdID', 'Dir'});
% Design stoichiometric matrix
% Find indices
for i = 1:height(mtbs_interest)
    mtb = mtbs_interest{i,'NameInDiagram'};
    id = find(strcmp(mtb, completeModel.metNames), 1);
    if ~isempty(id)
        rxn_list{strcmp(mtb, rxn_list{:,"Substrate"}), 'SubID'} = id;
        rxn_list{strcmp(mtb, rxn_list{:,"Product"}), 'ProdID'} = id;
    end
    clear id mtb
end

% Establish directions
rxn_list{strcmp('One-way', rxn_list{:,'Direction'}), 'Dir'} = 1;
rxn_list{strcmp('Two-way', rxn_list{:,'Direction'}), 'Dir'} = 0;

% Add to stoichiometric matrix
last_rxn = length(completeModel.rxns);
completeModel.rxns = cat(1, completeModel.rxns, rxn_list{:,"Enzyme"});
for i = 1:height(rxn_list)
    completeModel.S(rxn_list{i, "SubID"}, i + last_rxn) = -1;
    completeModel.S(rxn_list{i, "ProdID"}, i + last_rxn) = 1;
    if (rxn_list{i, "Dir"} == 0)
        % Cannot implement reversible reactions just yet
    end
end

%% Save model to file
save("Model files\completeModel.mat", 'completeModel');