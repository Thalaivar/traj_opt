clearvars
load('trajectories\eight_fft_001.mat')
n = (length(X)-2)/3;
T = X(3*n+1); VR = X(3*n+2);
p = 0.25;
[state, control] = get_traj([X(1:n), X(n+1:2*n), X(2*n+1:3*n)], T, VR, p);
X = [state(:,1); state(:,2); state(:,3); state(:,4); state(:,5); state(:,6); control(:,1); control(:,2); control(:,3); T; VR];
plotting(X, p);