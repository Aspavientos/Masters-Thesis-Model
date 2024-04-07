%% Information
% Author: Diego Rodriguez Esperante
% Contact: d.rodriguezesperante@student.maastrichtuniversity.nl
% Last update: 03/03/2024

% clear; clc; close all;

%% Load model, data, metadata
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

clear modelFileName folder opts_data file_meta file_data



%% Select metadata options
fig = uifigure('Name', 'My GUI', 'NumberTitle', 'off', 'Position', [100 100 300 200]);

% Create a panel to hold the buttons and labels
gridLayout = uigridlayout(fig, [7 2]);
gridLayout.ColumnWidth = {'1x', 'fit'};
gridLayout.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
gridLayout.Layout.Row = [1 2 2 3 4 4 5 6 6 7];
gridLayout.Layout.Column = [1 1 2 1 1 2 1 1 2 2];

% Model elements
labelModelTitle = uilabel(gridLayout, 'Position', [1 1], 'Text', 'Label 1');
labelModel = uilabel(gridLayout, 'Position', [1 1], 'Text', 'Label 1');
buttonModel = uibutton(gridLayout, 'Text', 'File...');

% Data elements
labelDataTitle = uilabel(gridLayout, 'Position', [120 50 100 30], 'Text', 'Label 2');
labelData = uilabel(gridLayout, 'Position', [120 50 100 30], 'Text', 'Label 2');
buttonData = uibutton(gridLayout, 'Position', [10 50 100 30], 'Text', 'File...');

% Meta elemenets
labelMetaTitle = uilabel(gridLayout, 'Position', [120 90 100 30], 'Text', 'Label 3');
labelMeta = uilabel(gridLayout, 'Position', [120 90 100 30], 'Text', 'Label 3');
buttonMeta = uibutton(gridLayout, 'Position', [10 90 100 30], 'Text', 'File...');

% Add callbacks to the buttons
set(buttonModel, 'ButtonPushedFcn', @(src, event) updateLabel(src, labelModel));
set(buttonData, 'ButtonPushedFcn', @(src, event) updateLabel(src, labelData));
set(buttonMeta, 'ButtonPushedFcn', @(src, event) updateLabel(src, labelMeta));

% Function to update the label text
function updateLabel(src, label)
    disp('Click!')
end


%% Select network & generate network