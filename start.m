% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CoSOILcal
% Cosmogenic soil production rate calculator
% for samples in a variable-density profile
% Angel Rodes, SUERC, 2019
% angelrodes@gmail.com
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clean
clear
close all hidden
% clc

%% choose profile data
xls_files = dir(fullfile('*.csv'));
XLS_files = dir(fullfile('*.CSV'));
xls_str = [{xls_files.name} {XLS_files.name}];
[s,v] = listdlg('PromptString','Select a .csv file:',...
    'SelectionMode','multiple',...
    'ListString',xls_str);
profilefile=xls_str(s);

for file=profilefile
    file2=file{1};
    soil_solver( file2 )
end