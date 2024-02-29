%% Information
% Author: Diego Rodriguez Esperante
% Contact: d.rodriguezesperante@student.maastrichtuniversity.nl
% Last update: 21/01/2024

%clear; clc; close all;

%% Initialize
% Following tutorial from: https://opencobra.github.io/cobratoolbox/stable/tutorials/tutorialPFBA.html
initCobraToolbox(false);

solverName = 'glpk';
solverType = 'LP';
changeCobraSolver(solverName, solverType);

%% Read files
addpath(pwd);

% Import Human1-GEM
modelFileName = ['Model files' filesep 'Human-GEM.mat'];
modelFileName= [pwd filesep modelFileName];
human1Model = readCbModel(modelFileName);

% Import Recon3D
modelFileName = ['Model files' filesep 'Recon3DModel_301.mat'];
modelFileName= [pwd filesep modelFileName];
recon3dModel = readCbModel(modelFileName);

% Import our rxns of interest
rxns_interest = readtable(['CSV' filesep 'Reactions - Rxn-Sub Pairs.csv']);
mtbs_interest = readtable(['CSV' filesep 'Reactions - Metabolites.csv']);

clear modelFileName solverName solverType

%% Manipulate Human1 model to fit our interest
% Find reactions
h1rxnID = unique(findRxnIDs(human1Model, rxns_interest{:,"Database"}));
h1rxnID(h1rxnID==0) = [];

% Find metabolites
h1mtbID = findMetsFromRxns(human1Model, human1Model.rxns(h1rxnID));
h1mtbID = findMetIDs(human1Model, h1mtbID);

% Find extra rxns and metabolites to exclude later
extraRxnh1 = setdiff(human1Model.rxns, human1Model.rxns(h1rxnID));
sh1Model = removeRxns(human1Model, extraRxnh1);

%% Manipulate Human1 model to fit our interest
% Find reactions
r3rxnID = unique(findRxnIDs(recon3dModel, rxns_interest{:,"Database"}));
r3rxnID(r3rxnID==0) = [];

% Find metabolites
r3mtbID = findMetsFromRxns(recon3dModel, recon3dModel.rxns(r3rxnID));
r3mtbID = findMetIDs(recon3dModel, r3mtbID);

% Find extra rxns and metabolites to exclude later
extraRxnr3 = setdiff(recon3dModel.rxns, recon3dModel.rxns(r3rxnID));
sr3Model = removeRxns(recon3dModel, extraRxnr3);

%% Merge the two models together
% Using COBRA's native mergeTwoModels

% Find duplicate metabolites between models, removed from R3 since it's the
% supporting database
mtbs_shared = ~(strcmp(mtbs_interest{:, 'Human1'}, '')+strcmp(mtbs_interest{:, 'Recon3D'}, '')); % Find metabolites with names in H1 and R3
mtbs_rosetta = mtbs_interest(mtbs_shared,:);

sr3Model.mets = erase(sr3Model.mets, "[" + lettersPattern(1) + "]");

for i = 1:length(mtbs_rosetta{:, 'Recon3D'})
    mtb_temp = strcmp(sr3Model.mets, mtbs_rosetta{i, 'Recon3D'});
    if ~isempty(mtb_temp)
        [sr3Model.mets{mtb_temp}] = deal(strcat(char(mtbs_rosetta{i, 'Human1'}), 'c'));
    end
    clear mtb_temp
end

clear mtbs_shared i
% Merge together
mergedModel = mergeTwoModels(sh1Model, sr3Model, 1);

%% Standardize names in model to be the same as CSV
merge_mets = extract(mergedModel.mets, lettersPattern(3) + digitsPattern(5));
for i = 1:length(merge_mets)
    std_id = find(strcmp(mtbs_interest{:, 'Human1'}, merge_mets(i)));
    if isempty(std_id)
        std_id = find(strcmp(mtbs_interest{:, 'Recon3D'}, merge_mets(i)));
    end
    if ~isempty(std_id)
        mergedModel.metNames(i) = mtbs_interest{std_id, "NameInDiagram"};
    end
end

%% Save model to file
% save(['Model files' filesep 'H1R3MergedModel.mat'], 'mergedModel');