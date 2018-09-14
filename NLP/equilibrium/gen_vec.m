function v = gen_vec(comp, N)
    v = zeros(N,N);
    j = 1;
    for i = 1:N:length(comp)
        v(:,j) = comp(i:i+N-1,1);
        j = j+1;
    end
end