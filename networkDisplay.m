%% Information
% Author: Diego Rodriguez Esperante
% Contact: d.rodriguezesperante@student.maastrichtuniversity.nl
% Last update: 03/03/2024

% clear; clc; close all;

%% Initialize
initCobraToolbox(false);

solverName = 'glpk';
solverType = 'LP';
changeCobraSolver(solverName, solverType);

clear solverName solverType

%% Read files
% Import polished model
foldername = 'expTurn_ChoDHEA-sinkAN';
filename = 'expTurn_Turner_ChoDHEA-sinkAN';
folder = ['Model files' filesep 'Flux balance results' filesep foldername];
modelFileName = [folder filesep filename '.mat'];
modelFileName = [pwd filesep modelFileName];
displayModel = readCbModel(modelFileName);

if(isfield(displayModel, 'subSystems'))
    displayModel = rmfield(displayModel, 'subSystems');
end

clear modelFileName folder

%% Output network in Cytoscape
filename = 'networkDisplay';
folder = ['Network files' filesep filename filesep foldername];

if ~exist(folder, 'dir')
    mkdir(folder);
end

outputNetworkCytoscape(displayModel, [folder filesep filename]);

% outputNetworkCytoscape(displayModel, [folder filesep filename], ...
%     displayModel.rxns, displayModel.rxnNames, ...
%     displayModel.mets, displayModel.metNames);

% Cytoscape node data
varnames = {'Node', 'Type', 'ProperName'};
tablemets = table(displayModel.mets, repmat({'met'}, length(displayModel.mets), 1), displayModel.metNames, 'VariableNames',varnames);
tablerxns = table(displayModel.rxns, repmat({'rxn'}, length(displayModel.rxns), 1), displayModel.rxnNames, 'VariableNames',varnames);
tablegenes = table(displayModel.genes, repmat({'gene'}, length(displayModel.genes), 1), displayModel.geneNames, 'VariableNames',varnames);

totaltable = [tablemets; tablerxns; tablegenes];

writetable(totaltable, [folder filesep filename '_NodeData.csv']);