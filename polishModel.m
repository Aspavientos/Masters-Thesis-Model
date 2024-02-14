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
modelFileName = ['Model files' filesep 'completeModel.mat'];
modelFileName= [pwd filesep modelFileName];
mergedModel = readCbModel(modelFileName);