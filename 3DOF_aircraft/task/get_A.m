function A = get_A(ac, t, X)
    I = eye(6); A = zeros(6);
    for i = 1:6
        h = 1e-6;
        f1 = ac.non_flat_model(t, X - h*I(:,i)); f2 = ac.non_flat_model(t, X + h*I(:,i));
        A(:,i) = (f2-f1)/(2*h);
    end
end