function vv = interp_sinc(t,v,tt)
    N = length(t);
    if mod(N,2)
        error('Length of t must be even.');
    end
    dt = diff(t);
    if abs(dt(1)-mean(dt))>1e-7&&abs(dt(end)-mean(dt))>1e-7
        error('t should be equi-spaced.');
    end
    if abs(v(1)-v(end))<1e-7
        error('v(1) must not be the same as v(N). Remove the leftmost data-point');
    end
    
    tmin = t(1)-(t(2)-t(1));
    tmax = t(end);
    x = 2*pi*(t-tmin)/(tmax-tmin);
    
    vv = zeros(size(tt));
    for j = 1:length(tt)
        vv(j) = 0;
        xx = 2*pi*(tt(j)-tmin)/(tmax-tmin);
        for i = 1:N            
            vv(j) = vv(j) + v(i)*cardFun(xx-x(i),N);
%             vv(j)
%             v(i)*cardFun(xx-x(i),N)
        end
    end
end
function c = cardFun(y,N)
    if abs(y)<1e-7||abs(y+2*pi)<1e-7
        c=1;
    else
        c=sin(0.5*N*y)/(N*tan(0.5*y));
    end
end