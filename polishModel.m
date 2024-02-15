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
%% Read filesk
% Import complete model
modelFileName = ['Model files' filesep 'completeModel.mat'];
modelFileName = [pwd filesep modelFileName];
completeModel = readCbModel(modelFileName);

clear modelFileName
%% Find same metabolites in different reactions and sequester them
db_mets.pat = contains(completeModel.mets, lettersPattern(3) + digitsPattern(5));
db_mets.extr = extract(completeModel.mets(db_mets.pat), lettersPattern(3) + digitsPattern(5));
[db_mets.red, db_mets.red_ia, ~] = unique(db_mets.extr);

% Work with no dupes model
modelNoDupes = removeMetabolites(completeModel, completeModel.mets(setdiff(1:length(db_mets.extr), db_mets.red_ia)));
dbndps_mets.pat = contains(modelNoDupes.mets, lettersPattern(3) + digitsPattern(5));
dbndps_mets.extr = extract(modelNoDupes.mets(dbndps_mets.pat), lettersPattern(3) + digitsPattern(5));
modelNoDupes.mets(dbndps_mets.pat) = dbndps_mets.extr;

% Work with all dupes model
modelAllDupes = removeMetabolites(completeModel, completeModel.mets([db_mets.red_ia; find(~db_mets.pat)]));
dbadps_mets.pat = contains(modelAllDupes.mets, lettersPattern(3) + digitsPattern(5));
dbadps_mets.extr = extract(modelAllDupes.mets(dbadps_mets.pat), lettersPattern(3) + digitsPattern(5));
modelAllDupes.mets(dbadps_mets.pat) = dbadps_mets.extr;

% Merge models
trueNoDupes = mergeTwoModels(modelNoDupes, modelAllDupes, 1);