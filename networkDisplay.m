%% Information
% Author: Diego Rodriguez Esperante
% Contact: d.rodriguezesperante@student.maastrichtuniversity.nl
% Last update: 03/03/2024

% clear; clc; close all;

%% Initialize
try prepareTest()
catch
    initCobraToolbox(false);

    solverName = 'glpk';
    solverType = 'LP';
    changeCobraSolver(solverName, solverType);

    clear solverName solverType

    disp('COBRA Toolbox initialized!')
end
%% Read files
% Import polished model
[fileName, pathname] = uigetfile();
modelFileName = [pathname fileName];
displayModel = readCbModel(modelFileName);

if(isfield(displayModel, 'subSystems'))
    displayModel = rmfield(displayModel, 'subSystems');
end

clear modelFileName
%% Output network in Cytoscape
[~, file] = fileparts(fileName);
filesections = split(file, '_');
folder_sub = regexp(pathname,'\','split');
folder = ['Network files' filesep 'Network Display' filesep folder_sub{end-2} filesep folder_sub{end-1} filesep filesections{2}];

if ~exist(folder, 'dir')
    mkdir(folder);
end

outputNetworkCytoscape(displayModel, [folder filesep file]);

% outputNetworkCytoscape(displayModel, [folder filesep filename], ...
%     displayModel.rxns, displayModel.rxnNames, ...
%     displayModel.mets, displayModel.metNames);

% Cytoscape node data
varnames = {'Node', 'Type', 'ProperName'};
tablemets = table(displayModel.mets, repmat({'met'}, length(displayModel.mets), 1), displayModel.metNames, 'VariableNames',varnames);
tablerxns = table(displayModel.rxns, repmat({'rxn'}, length(displayModel.rxns), 1), displayModel.rxnNames, 'VariableNames',varnames);
tablegenes = table(displayModel.genes, repmat({'gene'}, length(displayModel.genes), 1), displayModel.geneNames, 'VariableNames',varnames);

totaltable = [tablemets; tablerxns; tablegenes];

writetable(totaltable, [folder filesep file '_NodeData.csv']);

clear fileName pathname folder folder_sub file varnames

disp('#### Finished! ####')