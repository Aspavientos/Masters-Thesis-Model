%% Information
% Author: Diego Rodriguez Esperante
% Contact: d.rodriguezesperante@student.maastrichtuniversity.nl
% Last update: 03/03/2024

% clear; clc; close all;

% Following the FBA tutorial here: https://opencobra.github.io/cobratoolbox/stable/tutorials/tutorialFBA.html

%% Initialize
initCobraToolbox(false);

solverName = 'glpk';
solverType = 'LP';
changeCobraSolver(solverName, solverType);

clear solverName solverType
%% Read files
% Import polished model
modelFileName = ['Model files' filesep 'polishedModel.mat'];
modelFileName = [pwd filesep modelFileName];
polishedModel = readCbModel(modelFileName);

% Experimental data
folder = ['CSV' filesep 'Expression data'];
expA_data = readtable([folder filesep 'exprA_ENSEMBL.csv']);
expA_meta = readtable(['CSV' filesep 'exprA_meta.csv'], 'TextType','string');

% Model construction data
rxns_interest = readtable(['CSV' filesep 'Reactions - Rxn-Sub Pairs.csv']);

clear modelFileName folder
%% Common information
% List of genes
genelist = findGenesFromRxns(polishedModel, polishedModel.rxns);

% Meaningful data
subgroup = strcmp(expA_meta{:,"Ectopic"}, 'TRUE');
exprData.gene = expA_data(:,2:end).Properties.VariableNames';
exprData.med = median(expA_data{subgroup,2:end})';
exprData.max = max(expA_data{subgroup,2:end})';
exprData.min = min(expA_data{subgroup,2:end})';

% Source metabolites and objective reactions
sourcemet = 'MAM01450';

% Sink metabolites and objective reactions
sinkmet = 'MAM01787';
objctv = ['sink_' sinkmet];

%% Pseudo eFlux
% Set all default bounds
boundaryModel = changeRxnBounds(polishedModel, polishedModel.rxns, 10000, 'u');
boundaryModel = changeRxnBounds(boundaryModel, boundaryModel.rxns, -10000, 'l');

% Calculate and change boundary values
[expr_med, ~] = selectGeneFromGPR(boundaryModel, exprData.gene, exprData.med, GPRparser(boundaryModel));
[expr_max, ~] = selectGeneFromGPR(boundaryModel, exprData.gene, exprData.max, GPRparser(boundaryModel));
[expr_min, ~] = selectGeneFromGPR(boundaryModel, exprData.gene, exprData.min, GPRparser(boundaryModel));
boundaryModel = changeRxnBounds(boundaryModel, boundaryModel.rxns, expr_med, 'u');
boundaryModel = changeRxnBounds(boundaryModel, boundaryModel.rxns, 0, 'l');

% Add demand/sink reactions
boundaryModel = addSinkReactions(boundaryModel, sinkmet, -10000, 10000);
boundaryModel = addSinkReactions(boundaryModel, sourcemet, -10000, 10000);

% Change objective function
boundaryModel = changeObjective(boundaryModel, objctv);

% Flux-balance analysis
optCB_sol = optimizeCbModel(boundaryModel, 'max');

% Print results
bound_v_sol = [boundaryModel.ub, optCB_sol.v];
printFluxVector(boundaryModel, optCB_sol.v, 1)

% Save results to csv
folder = ['CSV' filesep 'Flux Balance Results'];
testname = 'PlusMedZero';
csv_table = table(boundaryModel.rxns, boundaryModel.ub, optCB_sol.v, 'VariableNames', {'Reactions', 'Upper bounds', 'Flux Balance'});

writetable(csv_table, [folder filesep testname '.csv']);

% Save model for display later
folder = ['Model files' filesep 'Flux Balance Results'];;
save(['Model files' filesep testname '.mat'], 'boundaryModel');