N = 24;
[D, x] = fourierDiff(N);

% smooth function in [0, 2*pi]
f = exp(sin(x)); fprime = cos(x).*f;
% plot(x, D*f' - fprime', '-m', 'Marker', 'o');

% smooth function in [0, T]
T = 3;
t = 3*x/(2*pi);
f = exp(sin(t*2*pi/3)); fprime = (2*pi/3)*cos(t*2*pi/3).*f;
plot(t, (2*pi*D*f')/3 - fprime', '-m', 'Marker', 'o');