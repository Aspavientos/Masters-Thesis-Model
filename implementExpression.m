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
file_data = [folder filesep 'GDS5211_expr.csv'];
file_meta = [folder filesep 'GDS5211_meta.csv'];

opts_data = detectImportOptions(file_data);
[opts_data.VariableTypes{2:end}] = deal('double');

exp_data = readtable(file_data, opts_data);
exp_meta = readtable(file_meta, 'TextType','string');

% Model construction data
rxns_interest = readtable(['CSV' filesep 'Reactions - Rxn-Sub Pairs.csv']);

clear modelFileName folder opts_data file_meta file_data
%% Common information
% List of genes
genelist = findGenesFromRxns(polishedModel, polishedModel.rxns);

% Meaningful data
subgroup = strcmp(exp_meta{:,"TurnerSyndrome"}, 'normal euploid');
exprData.gene = exp_data(:,2:end).Properties.VariableNames';
exprData.med = median(exp_data{subgroup,2:end}, "omitnan")';
exprData.max = max(exp_data{subgroup,2:end}, [], "omitnan")';
exprData.min = min(exp_data{subgroup,2:end}, [], "omitnan")';

% Source metabolites and objective reactions
sourcemet = {'MAM01450', 'MAM01660'};

% Sink metabolites and objective reactions
sinkmet = {'MAM01338'};
objctv = strcat('sink_', sinkmet{1});

clear subgroup
%% Pseudo eFlux
% Set all default bounds
boundaryModel = changeRxnBounds(polishedModel, polishedModel.rxns, 10000, 'u');
boundaryModel = changeRxnBounds(boundaryModel, boundaryModel.rxns, -10000, 'l');

% Calculate and change boundary values
[expr_med, ~] = selectGeneFromGPR(boundaryModel, exprData.gene, exprData.med, GPRparser(boundaryModel));
[expr_max, ~] = selectGeneFromGPR(boundaryModel, exprData.gene, exprData.max, GPRparser(boundaryModel));
[expr_min, ~] = selectGeneFromGPR(boundaryModel, exprData.gene, exprData.min, GPRparser(boundaryModel));
boundaryModel = changeRxnBounds(boundaryModel, boundaryModel.rxns, expr_max, 'u');
boundaryModel = changeRxnBounds(boundaryModel, boundaryModel.rxns, 0, 'l');

% Add demand/sink reactions
for i = 1:length(sourcemet)
    boundaryModel = addSinkReactions(boundaryModel, sourcemet{i}, -10000, 0);
end
for i = 1:length(sinkmet)
    boundaryModel = addSinkReactions(boundaryModel, sinkmet{i}, 0, 10000);
end
% Change objective function
boundaryModel = changeObjective(boundaryModel, objctv);

% Flux-balance analysis
optCB_sol = optimizeCbModel(boundaryModel, 'max');

% Print results
bound_v_sol = [boundaryModel.ub, optCB_sol.v];
printFluxVector(boundaryModel, optCB_sol.v, 1)

clear sinkmet sourcemet objctv 
%% Save data
% Save results to csv
experiment = 'expTurn';
cohort = 'NotTurner';
test = 'ChoDHEA-sinkAN';
foldername = [experiment '_' test];
filename = [experiment '_' cohort '_' test];

folder = ['CSV' filesep 'Flux Balance Results' filesep foldername];
if ~exist(folder, 'dir')
    mkdir(folder);
end

fract = optCB_sol.v./boundaryModel.ub;
fract(isnan(fract)) = 0;
csv_table = table(boundaryModel.rxns, boundaryModel.ub, optCB_sol.v, fract, ...
    'VariableNames', {'Reactions', ['UB ' cohort], ['FBA ' cohort], ['Frc ' cohort]});

writetable(csv_table, [folder filesep filename '_FBA.csv']);

% Save model for display later
folder = ['Model files' filesep 'Flux Balance Results' filesep foldername];
if ~exist(folder, 'dir')
    mkdir(folder);
end
save([folder filesep filename '.mat'], 'boundaryModel');

clear experiment cohort test filename folder fract i