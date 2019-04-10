clear all

load eigData

[U,S,V] = svd(Dmat-M,'econ');

Z = U(:,1:d)*S(1:d,1:d)*(V(:,1:d)');

ee = eig(Z);