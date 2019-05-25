% to generate initial guess of iniclined circles
% generates the save files

addpath('../lib/')
N = 500;

% circle guess
VR_0 = 0.1; tf_0 = 10;
[~,x] = fourierdiff(N);
t = tf_0*x/(2*pi);
x = 20*(-1 + cos(2*pi*t/tf_0));
y = -20*(2^0.5)*sin(2*pi*t/tf_0);
z = -20*(1 - cos(2*pi*t/tf_0)) - 0.15;   

%% full state discretization
p = 0.25;
[X,U] = getTrajFFT([x', y',z'], tf_0, VR_0, p);
sol = [X(:,1); unwrap(X(:,2)); X(:,3); X(:,4); X(:,5); X(:,6); U(:,1); U(:,2); tf_0; VR_0];
sol(8*N+3) = -2*pi/tf_0;
saveFile = ['../full_state_discrete/solutions/IG_circle_expo_', num2str(N)];
save(saveFile, 'sol', 'N', 'p');

p = 1;
[X,U] = getTrajFFT([x', y',z'], tf_0, VR_0, 1);
sol = [X(:,1); unwrap(X(:,2)); X(:,3); X(:,4); X(:,5); X(:,6); U(:,1); U(:,2); tf_0; VR_0];
sol(8*N+3) = -2*pi/tf_0;
saveFile = ['../full_state_discrete/solutions/IG_circle_linear_', num2str(N)];

%% differential flatness discretization
sol = [x'; y'; z'; tf_0; VR_0];

p = 0.25;
saveFile = ['../diff_flat_discrete/solutions/IG_circle_expo_', num2str(N)];
save(saveFile, 'sol', 'N', 'p');

p = 1;
saveFile = ['../diff_flat_discrete/solutions/IG_circle_linear_', num2str(N)];
save(saveFile, 'sol', 'N', 'p');

%% eight guess
N = 60; 

VR_0 = 0.1; tf_0 = 20;
[~, x] = fourierdiff(N);
t = tf_0*x/(2*pi);
x = -20*sin(4*pi*t/tf_0);
y = -40*sqrt(2)*sin(2*pi*t/tf_0);
z = -20*(sin(4*pi*t/tf_0) + 1) - 0.11; 

%% full state discretization
p = 0.25;
[X,U] = getTrajFFT([x', y',z'], tf_0, VR_0, 0.25);
sol = [X(:,1); unwrap(X(:,2)); X(:,3); X(:,4); X(:,5); X(:,6); U(:,1); U(:,2); tf_0; VR_0];
saveFile = ['../full_state_discrete/solutions/IG_eight_expo_', num2str(N)];
save(saveFile, 'sol', 'N', 'p');

p = 1;
[X,U] = getTrajFFT([x', y',z'], tf_0, VR_0, 1);
sol = [X(:,1); unwrap(X(:,2)); X(:,3); X(:,4); X(:,5); X(:,6); U(:,1); U(:,2); tf_0; VR_0];
saveFile = ['../full_state_discrete/solutions/IG_eight_linear_', num2str(N)];
save(saveFile, 'sol', 'N', 'p');

%% differential flatness discretization
sol = [x'; y'; z'; tf_0; VR_0];

p = 0.25;
saveFile = ['../diff_flat_fourier/solutions/IG_eight_expo_', num2str(N)];
save(saveFile, 'sol', 'N', 'p');

p = 1;
saveFile = ['../diff_flat_fourier/solutions/IG_eight_linear_', num2str(N)];
save(saveFile, 'sol', 'N', 'p');