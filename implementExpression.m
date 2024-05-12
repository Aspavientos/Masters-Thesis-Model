%% Information
% Author: Diego Rodriguez Esperante
% Contact: d.rodriguezesperante@student.maastrichtuniversity.nl
% Last update: 03/03/2024

% clear; clc; close all;

% Following the FBA tutorial here: https://opencobra.github.io/cobratoolbox/stable/tutorials/tutorialFBA.html

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
modelFileName = ['Model files' filesep 'polishedModel.mat'];
modelFileName = [pwd filesep modelFileName];
polishedModel = readCbModel(modelFileName);

% Experimental data
folder = ['CSV' filesep 'Expression data'];
file_data = [folder filesep 'GTEx_expr.csv'];
file_meta = [folder filesep 'GTEx_meta.csv'];

opts_data = detectImportOptions(file_data);
[opts_data.VariableTypes{2:end}] = deal('double');

exp_data = readtable(file_data, opts_data);
exp_meta = readtable(file_meta, 'TextType','string', 'NumHeaderLines', 0, 'Delimiter', 'comma');

% Model construction data
rxns_interest = readtable(['CSV' filesep 'Reactions - Rxn-Sub Pairs.csv']);

% Choose experiment options
experiment = 'expGTEx';
cohort = {'Adrenal', 'Testis', 'Ovary'};
% {'Lung', 'Heart', 'Brain', 'Skin', 'Visceral', 'Prostate', 'Mammary', 'Uterus', 'Muscle'}
% {'Adrenal', 'Testis', 'Ovary'}
sinkMode = 'overr'; % 'simple', 'DRAIN', 'overr'

% Source metabolites and objective reactions
sourcemet = {'MAM01450'}; % DHEA: MAM01660. Cho: MAM01450

% Sink metabolites and objective reactions
sinkmet = {'MAM01338', 'MAM02969', 'MAM01787', 'MAM01069', 'MAM01615', '3αDIOL', '3βDIOL', 'allopregnandiol', 'MAM01614', 'MAM01660'};

clear modelFileName folder opts_data file_meta file_data
%% FBA calculations
firstTissue = 1;
for k = 1:length(cohort)
    %% Common information
    % List of genes
    genelist = findGenesFromRxns(polishedModel, polishedModel.rxns);

    % Aggregate experimental data according to cohorts
    subgroup = strcmp(exp_meta{:,"Tissue"}, cohort{k});
    [~, sub_ids] = intersect(exp_data{:, 1}, exp_meta{subgroup, 1});
    exprData.gene = exp_data(:,2:end).Properties.VariableNames';
    exprData.med = median(exp_data{sub_ids,2:end}, 1, "omitnan")';
    exprData.max = max(exp_data{sub_ids,2:end}, [], 1, "omitnan")';
    exprData.min = min(exp_data{sub_ids,2:end}, [], 1, "omitnan")';

    clear subgroup sub_ids
    %% Pseudo eFlux
    % Set all default bounds
    boundaryModel = changeRxnBounds(polishedModel, polishedModel.rxns, 10000, 'u');
    boundaryModel = changeRxnBounds(boundaryModel, boundaryModel.rxns, 0, 'l');

    % Calculate and change boundary values
    [expr_med, ~] = selectGeneFromGPR(boundaryModel, exprData.gene, exprData.med, GPRparser(boundaryModel));
    [expr_max, ~] = selectGeneFromGPR(boundaryModel, exprData.gene, exprData.max, GPRparser(boundaryModel));
    [expr_min, ~] = selectGeneFromGPR(boundaryModel, exprData.gene, exprData.min, GPRparser(boundaryModel));
    boundaryModel = changeRxnBounds(boundaryModel, boundaryModel.rxns, expr_med, 'u');
    boundaryModel = changeRxnBounds(boundaryModel, boundaryModel.rxns, 0, 'l');

    % Add demand reactions
    for i = 1:length(sourcemet)
        boundaryModel = addSinkReactions(boundaryModel, sourcemet{i}, -10000, 0);
    end

    switch sinkMode
        case {'simple', 'overr'}
            % Without faucet reactions
            for i = 1:length(sinkmet)
                boundaryModel = addSinkReactions(boundaryModel, sinkmet{i}, 0, 10000);
            end

            objctv = strcat('sink_', sinkmet{1});
        case 'DRAIN'
            % Add sink reactions all connected to a joined faucet reaction
            boundaryModel = addMetabolite(boundaryModel, 'DRAIN', 'DRAIN');
            boundaryModel = addMultipleReactions(boundaryModel, strcat({'drain_'}, sinkmet), [sinkmet 'DRAIN'], [-eye(length(sinkmet)); ones(1, length(sinkmet))]);
            boundaryModel = changeRxnBounds(boundaryModel, strcat({'drain_'}, sinkmet), 10000, 'u');
            boundaryModel = changeRxnBounds(boundaryModel, strcat({'drain_'}, sinkmet), 0, 'l');
            boundaryModel = addSinkReactions(boundaryModel, 'DRAIN', 0, 10000);

            objctv = 'sink_DRAIN';
    end

    % Flux-balance analysis
    switch sinkMode
        case {'simple', 'DRAIN'}
            % Change objective function
            boundaryModel = changeObjective(boundaryModel, objctv);
            % Solve model
            optCB_sol = optimizeCbModel(boundaryModel, 'max');
        case 'overr'
            optCB_sol = struct('v', zeros(length(boundaryModel.rxns), 1), ...
                'overr', zeros(length(boundaryModel.rxns), 1));
            for i = 1:length(sinkmet)
                objctv = strcat('sink_', sinkmet{i});
                % Change objective function
                boundaryModel = changeObjective(boundaryModel, objctv);
                % Solve model
                optCB_soltemp = optimizeCbModel(boundaryModel, 'max');
                optCB_sol.(genvarname(sinkmet{i})) = optCB_soltemp.v;
                optCB_sol.v = optCB_sol.v + optCB_soltemp.v;
                optCB_sol.overr = optCB_sol.overr + (optCB_soltemp.v ~= 0);
            end
            clear i optCB_soltemp
    end


    % Print results
    bound_v_sol = [boundaryModel.ub, optCB_sol.v];
    printFluxVector(boundaryModel, optCB_sol.v, 1)

    clear objctv
    %% Save data
    % Save results to csv
    foldername = [sinkMode filesep experiment];
    filename = [experiment '_' cohort{k} '_' sinkMode];

    folder = ['CSV' filesep 'Flux Balance Results' filesep foldername];
    if ~exist(folder, 'dir')
        mkdir(folder);
    end

    csv_table = table(boundaryModel.rxns, boundaryModel.ub, optCB_sol.v, ...
        'VariableNames', {'Reactions', ['UB ' cohort{k}], ['FBA ' cohort{k}]});

    switch sinkMode
        case 'overr'
            overr_table = table(optCB_sol.overr, 'VariableNames', {['Overr ' cohort{k}]});
            csv_table = [csv_table, overr_table];
    end
    writetable(csv_table, [folder filesep filename '_FBA.csv']);

    % If this is the first tissue examined, save model for display later
    if firstTissue
        folder = ['Model files' filesep 'Flux Balance Results' filesep foldername];
        modelfilename = [experiment '_' strjoin(extractBetween(cohort, 1, 1), '') '_' sinkMode];
        if ~exist(folder, 'dir')
            mkdir(folder);
        end
        save([folder filesep modelfilename '.mat'], 'boundaryModel');
        firstTissue = 0;
    end
    clear filename folder fract i modelfilename
end
clear experiment cohort sinkMode sinkmet sourcemet firstTissue k
disp('#### Finished! ####')