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



%% Select metadata options
clear; clc;
fig = uifigure('Name', 'My GUI', 'NumberTitle', 'off', 'Position', [100 100 400 250]);
fig.UserData = struct('Model', [], 'Data', [], 'Meta', []);

% Create a panel to hold the buttons and labels
gridLayout = uigridlayout(fig, [7 2]);
gridLayout.ColumnWidth = {'1x', 'fit'};
gridLayout.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};

% Model elements
labelModelTitle = uilabel(gridLayout, 'Text', 'MODEL TITLE', ...
    'FontSize', 16, 'FontWeight','bold');
labelModelTitle.Layout.Column = [1,2];
labelModel = uilabel(gridLayout, 'Text', 'Select Model');
buttonModel = uibutton(gridLayout, 'Text', 'File...', 'Tag', 'Model');

% Data elements
labelDataTitle = uilabel(gridLayout, 'Text', 'DATA TITLE', ...
    'FontSize', 16, 'FontWeight','bold');
labelDataTitle.Layout.Column = [1,2];
labelData = uilabel(gridLayout, 'Text', 'Select Data');
buttonData = uibutton(gridLayout, 'Text', 'File...', 'Tag', 'Data');

% Meta elemenets
labelMetaTitle = uilabel(gridLayout,'Text', 'META TITLE', ...
    'FontSize', 16, 'FontWeight','bold');
labelMetaTitle.Layout.Column = [1,2];
labelMeta = uilabel(gridLayout, 'Text', 'Select Meta');
buttonMeta = uibutton(gridLayout, 'Text', 'File...', 'Tag', 'Meta');

% Next button
buttonNext = uibutton(gridLayout, 'Text', 'Next');
buttonNext.Layout.Column = [1,2];

% Add callbacks to the buttons
set(buttonModel, 'ButtonPushedFcn', @(src, event) selectFile(src, labelModel, fig));
set(buttonData, 'ButtonPushedFcn', @(src, event) selectFile(src, labelData, fig));
set(buttonMeta, 'ButtonPushedFcn', @(src, event) selectFile(src, labelMeta, fig));
set(buttonNext, 'ButtonPushedFcn', @(src, event) continueNext(src, fig));

% Functions
function selectFile(src, label, figure)
    [file, path] = uigetfile('*.*');
    if file == 0
        return
    end
    set(label, 'Text', file, 'FontAngle', 'italic', 'FontColor', 'b');
    new_struct = figure.UserData;
    new_struct.(src.Tag) = {file, path};
    set(figure, 'UserData' , new_struct);
    disp([src.Tag ' selected!'])
end

function continueNext(src, figure)
    if any(structfun(@isempty, figure.UserData))
        set(src, 'Text', 'Please select all files', 'FontColor', 'r');
        %return
    end
    set(src, 'Text', 'Loading...', 'Enable', 'off');
    try
        prepareTest();
    catch
        initCobraToolbox(false);
        solverName = 'glpk';
        solverType = 'LP';
        changeCobraSolver(solverName, solverType);
    end
    modelFileName = [figure.UserData.Model{2} figure.UserData.Model{1}];
    polishedModel = readCbModel(modelFileName);
end

%% Select network & generate network