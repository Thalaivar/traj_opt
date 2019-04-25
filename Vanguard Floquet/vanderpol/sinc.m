function y = sinc(x,h)
    y = sin(pi*x/h)/((2*pi/h)*tan(x/2));
end