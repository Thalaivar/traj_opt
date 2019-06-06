clear all

% ref: Bower et al.
% m = 4;
% b = 2.492;
% S = 0.485;
% CD0 = 0.017;
% CD1 = 0.001854;
% CD2 = 0.04388;

% ref: Deittert et al.
m = 4.5;
b = 3;
S = 0.473;
CD0 = 0.0173;
CD1 = -0.0337;
CD2 = 0.0517;
p_exp = NaN;

ac = aircraft(m,b,S,CD0,CD1,CD2,p_exp);

save('acmod','ac');