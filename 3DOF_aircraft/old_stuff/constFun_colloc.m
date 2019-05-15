function [c, ceq] = constFun_colloc(X, ac)
    d = 4;
    N = (length(X)-2)/8;
    
    th_eigvec_re = zeros(N,d);
    th_eigvec_im = zeros(N,d);
    for i = 1:d
        j = (i-1)*N + 1;
        th_eigvec_re(:,i) = X(j:j+(N-1),1);
        th_eigvec_im(:,i) = X(4*N+j:4*N+j+(N-1),1);
    end
    eigval = complex(X(8*N+1,1), X(8*N+2,1));
    
    [D, cheb_x] = cheb_diff(N-1);
    t = 0.5*ac.tf*(1-cheb_x);
    
    % get derivatives from dynamics
    dot_eigvec_re = zeros(N,d);
    dot_eigvec_im = zeros(N,d);
    for i = 1:N
        A = ac.get_jac(t(i), 'FD');
        A = [A(1:3,1:3),A(1:3,6);A(6,1:3),A(6,6)];
        dot_eigvec_re(i,:) = ((A - real(eigval)*eye(4))*th_eigvec_re(i,:)' + imag(eigval)*th_eigvec_im(i,:)')';
        dot_eigvec_im(i,:) = ((A - real(eigval)*eye(4))*th_eigvec_im(i,:)' - imag(eigval)*th_eigvec_re(i,:)')';
    end
    Xdot = zeros(N,1);
    for i = 1:d
        j = (i-1)*N + 1;
        Xdot(j:j+(N-1),1) = dot_eigvec_re(:,i);
        Xdot(4*N+j:4*N+j+(N-1),1) = dot_eigvec_im(:,i);
    end
    
    % get spectral derivative
    Xdot_spec = zeros(N,1);
    for i = 1:d
        j = (i-1)*N + 1;
        Xdot_spec(j:j+(N-1),1) = (-2/ac.tf)*D*th_eigvec_re(:,i);
        Xdot_spec(4*N+j:4*N+j+(N-1),1) = (-2/ac.tf)*D*th_eigvec_im(:,i);
    end
    
    ceq = Xdot_spec - Xdot;
    
    ceq(end+1:end+4,1) = th_eigvec_re(1,:)' - th_eigvec_re(end,:)';
    ceq(end+1:end+4,1) = th_eigvec_im(1,:)' - th_eigvec_im(end,:)';
    
    ceq(end+1,1) = norm(complex(th_eigvec_re(1,:)', th_eigvec_im(1,:)')) - 1;
    c = [];
end