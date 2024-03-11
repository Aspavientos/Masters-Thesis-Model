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
expA_data = readtable(['CSV' filesep 'exprA_ENSEMBL.csv']);
expA_meta = readtable(['CSV' filesep 'exprA_meta.csv'], 'TextType','string');

% Model construction data
rxns_interest = readtable(['CSV' filesep 'Reactions - Rxn-Sub Pairs.csv']);

clear modelFileName
%% Common information
% List of genes
genelist = findGenesFromRxns(polishedModel, polishedModel.rxns);

% Meaningful data
subgroup = strcmp(expA_meta{:,"Ectopic"}, 'TRUE');
exprData.gene = expA_data(:,2:end).Properties.VariableNames';
exprData.med = median(expA_data{subgroup,2:end})';
exprData.max = max(expA_data{subgroup,2:end})';
exprData.min = min(expA_data{subgroup,2:end})';

% Sink metabolites and objective reactions
sinkmet = 'MAM01787';
objctv = ['sink_' sinkmet];

%% Expression as proxy for boundaries
% Set all default bounds
boundaryModel = changeRxnBounds(polishedModel, polishedModel.rxns, 10000, 'u');
boundaryModel = changeRxnBounds(boundaryModel, boundaryModel.rxns, 0, 'l');

% Calculate and change boundary values
[expr_max, ~] = selectGeneFromGPR(boundaryModel, exprData.gene, exprData.max, GPRparser(boundaryModel));
[expr_min, ~] = selectGeneFromGPR(boundaryModel, exprData.gene, exprData.min, GPRparser(boundaryModel));
boundaryModel = changeRxnBounds(boundaryModel, boundaryModel.rxns, expr_max, 'u');
boundaryModel = changeRxnBounds(boundaryModel, boundaryModel.rxns, 0, 'l');

% Add demand/sink reactions
boundaryModel = addSinkReactions(boundaryModel, sinkmet, 0, 10000);

% Change objective function
boundaryModel = changeObjective(boundaryModel, objctv);

% Flux-balance analysis
optCB_sol = optimizeCbModel(boundaryModel, 'max');
%enOpt_sol = enumerateOptimalSolutions(boundaryModel);
rndObj_sol = randomObjFBASol(boundaryModel, {'max', 0});

% Print results
printFluxVector(boundaryModel, optCB_sol.v, 1)