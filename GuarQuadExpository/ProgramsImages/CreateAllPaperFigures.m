%Run all functions to produce figures for the paper submitted to the
%American Mathematical Monthly

%% Garbage collection
format compact %remove blank lines from output
format long %lots of digits
clearvars %delete all variables
close all %close all figures
set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
set(0,'defaultLineLineWidth',3) %thick lines
set(0,'defaultTextInterpreter','latex') %latex axis labels
set(0,'defaultLegendInterpreter','latex') %latex axis labels
set(0,'defaultLineMarkerSize',40) %latex axis labels

%% Create color figures
bwcolor='color';
GaussProbTrapExample(bwcolor)
SpikesFlukes(bwcolor)
FoolIntegral(bwcolor)
TrianglePeak(bwcolor)

%% Create color figures
bwcolor='bw';
GaussProbTrapExample(bwcolor)
SpikesFlukes(bwcolor)
FoolIntegral(bwcolor)
TrianglePeak(bwcolor)