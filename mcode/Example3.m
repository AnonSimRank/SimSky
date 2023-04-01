function Example3()

% test the actual error and upper bound, iterative error in Eq.16

c=0.8;
theta=1;
A = [0 0 1/2 0 0 0;1 0 0 1 0 0;0 1/2 0 0 0 0;0 0 1/2 0 0 1;0 1/2 0 0 0 0;0 0 0 0 1 0];
n=size(A,1);
sigma = max(abs(eigs(sqrt(c)*A)))+eps;

j=3;
k=5;

v=sparse(j,1,1,n,1);
theta^2 * sigma^(2*k+1) / (1-sigma^2) *normest(A * v) * c^0.5 + theta^2 * c^(k+2) * c^0.5 * (1-sigma^(2*k+2)) / (sigma - sigma^3) * normest(A*v)
soar_ss_kk = single_source_simrank(A, v, c, 30);
soar_ss_k = single_source_simrank(A, v, c, 5);
norm(soar_ss_kk - soar_ss_k)
