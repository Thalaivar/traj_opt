clearvars
close all

addpath('../DS/lib/')
addpath('spectralFEconverge\')

global fonttype
global fontsize
global markerSZ

fonttype = 'Computer Modern';
fontsize = 14;
markerSZ = 80;

rawResults
diffNVDP
goal1StiffPlot
FSRK4ConsistencyPlot
plotError