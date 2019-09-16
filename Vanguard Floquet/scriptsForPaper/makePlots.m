clearvars
close all

addpath('../DS/lib/')
addpath('spectralFEconverge\')

global fonttype
global fontsize
global markerSZ

fonttype = 'Times New Roman';
fontsize = 14;
markerSZ = 80;

rawResults
diffNVDP
goal1StiffPlot
FSRK4ConsistencyPlot
PSConsistencyPlot
accuracyPlot
plotError