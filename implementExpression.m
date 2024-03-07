%% Information
% Author: Diego Rodriguez Esperante
% Contact: d.rodriguezesperante@student.maastrichtuniversity.nl
% Last update: 03/03/2024

%clear; clc; close all;

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

expA_data = readtable(['CSV' filesep 'exprA_ENSEMBL.csv']);

expA_meta = readtable(['CSV' filesep 'exprA_meta.csv'], 'TextType','string');

clear modelFileName
%% Set boundaries
% Set all bounds to default
polishedModel = changeRxnBounds(polishedModel, polishedModel.rxns, 0, 'b');
genelist = findGenesFromRxns(polishedModel, polishedModel.rxns);

subgroup = strcmp(expA_meta{:,"Ectopic"}, 'TRUE');
exprData.gene = expA_data(:,2:end).Properties.VariableNames';
exprData.value = median(expA_data{subgroup,2:end})';
% Calculate boundary values
[expr, ~] = selectGeneFromGPR(polishedModel, exprData.gene, exprData.value, GPRparser(polishedModel));
expr_bound = expr./median(expr)*1000;
polishedModel = changeRxnBounds(polishedModel, polishedModel.rxns, expr_bound, 'u');


%% Change objective function

%% Add demand/sink reactions