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
% Import complete model
modelFileName = ['Model files' filesep 'completeModel.mat'];
modelFileName = [pwd filesep modelFileName];
completeModel = readCbModel(modelFileName);

% Import metabolite data
mtbs_interest = readtable(['CSV' filesep 'Reactions - Metabolites.csv']);

% Import gene association data
enz_interest = readtable(['CSV' filesep 'Reactions - Enzymes.csv'], 'Delimiter', 'comma');

% Import reaction data
rxns_interest = readtable(['CSV' filesep 'Reactions - Rxn-Sub Pairs.csv']);

clear modelFileName
%% Remove unused metabolites
mtbs_keep = mtbs_interest{strcmp(mtbs_interest{:,"OfInterest"},'Yes'),'NameInDiagram'};
[mtbs_rem, ia] = setdiff(completeModel.metNames, mtbs_keep);

completeModel = removeFieldEntriesForType(completeModel, ia, 'mets', numel(completeModel.mets));

%% Fix gene-rxn associations
% Remove previous genes to start from blank slate
generxnAss = removeGenesFromModel(completeModel, completeModel.genes);

% Add genes of interest
gen = enz_interest.EnsemblID;
gen(strcmp(gen, '')) = enz_interest.NameInDiagram(strcmp(gen, ''));
generxnAss = addGenes(generxnAss, gen, 'geneNames', enz_interest.NameInDiagram);

% Create grRules manually
grules = strings(length(generxnAss.rxns),1);
for i = 1:length(generxnAss.rxns)
    orig = find(strcmp(rxns_interest.Name, generxnAss.rxns(i)), 1);
    sub = rxns_interest.Substrate{orig};
    prod = rxns_interest.Product{orig};

    redun = find(strcmp(rxns_interest.Substrate, sub) & strcmp(rxns_interest.Product, prod));
    geneid = find(strcmp(rxns_interest.Enzyme(redun(1)), enz_interest.NameInDiagram));
    if ~strcmp(enz_interest.EnsemblID(geneid), '')
        grules(i) = enz_interest.EnsemblID(geneid);
    else
        grules(i) = enz_interest.NameInDiagram(geneid);
    end

    for j = 2:length(redun)
        geneid = find(strcmp(rxns_interest.Enzyme(redun(j)), enz_interest.NameInDiagram));
        if ~strcmp(enz_interest.EnsemblID(geneid), '')
            grules(i) = append(grules(i), ' or ', enz_interest.EnsemblID(geneid));
        else
            grules(i) = append(grules(i), ' or ', enz_interest.NameInDiagram(geneid));
        end
    end
    generxnAss = changeGeneAssociation(generxnAss, generxnAss.rxns(i), grules{i});
    generxnAss.rxnNames{i} = [sub ' to ' prod];
end

clear gen i j orig sub prod redun geneid
%% Edit model metadata
% Remove non-important fields
metaModel = generxnAss;
kept_fields = {'id', ...
    'mets', 'metNames', 'b', 'csense',  'c', 'lb', 'ub', ...
    'rxns', 'rxnNames', ...
    'genes', 'geneNames', ...
    'rxnGeneMat', 'grRules', 'rules', 'S'};
exclude_fields = setdiff(fieldnames(metaModel), kept_fields);
missing_fields = setdiff(kept_fields, fieldnames(metaModel));

for i = 1:length(missing_fields) % Add missing
    metaModel.(missing_fields{i}) = [];
end
metaModel = rmfield(metaModel, exclude_fields); % Remove extra
metaModel = orderfields(metaModel, kept_fields); % Order to keep nice

% Add relevant metadata
metaModel.id = char('Human1-Recon3D Endometrium Subset', 'Developed by: Diego Rodriguez', 'Mail: diegoeldelccm@gmail.com');

clear i
%% Perform sanity checks
verifyModel(metaModel);
polishedModel = metaModel;

%% Save to file
% save(['Model files' filesep 'polishedModel.mat'], 'polishedModel');