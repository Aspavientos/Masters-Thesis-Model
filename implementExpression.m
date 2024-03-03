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

clear modelFileName
%% Set boundaries
% Set all bounds to default
polishedModel = changeRxnBounds(polishedModel, polishedModel.rxns, 0, 'b');
genelist = findGenesFromRxns(polishedModel, polishedModel.rxns);

% Calculate boundary values
% exprbound = ;
polishedModel = changeRxnBounds(polishedModel, polishedModel.rxns, exprbound, 'b');

%% Change objective function

%% Add demand/sink reactions