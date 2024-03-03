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

%% Output network in Cytoscape
folder = 'Network files';
outputNetworkCytoscape(polishedModel, [folder filesep 'networkPolished']);