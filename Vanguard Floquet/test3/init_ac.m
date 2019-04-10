clear all

addpath('./AllData')
load data
load prm
% data_traj = data_lin;
data_traj = data_expo;

m=data_traj.params.m;
b=data_traj.params.b;
S=data_traj.params.S;
CD0=data_traj.params.Cd0;
CD1=data_traj.params.Cd1;
CD2=data_traj.params.Cd2;
VR=data_traj.VR;
p_exp=data_traj.p;
NCT=1;
ac = aircraft(m,b,S,CD0,CD1,CD2,NCT,p_exp,VR);

prm.m=data_traj.params.m;
prm.b=data_traj.params.b;
prm.S=data_traj.params.S;
prm.CD0=data_traj.params.Cd0;
prm.CD1=data_traj.params.Cd1;
prm.CD2=data_traj.params.Cd2;
prm.p_exp=data_traj.p;
prm.VR = data_traj.VR;

save('./AllData/AllData','prm','data_traj','ac');

rmpath('./AllData')