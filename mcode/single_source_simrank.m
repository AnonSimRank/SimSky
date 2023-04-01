function s_k_j = single_source_simrank(A, v, c, kmax)
% single_source_simrank compute Optimized single-source SimRank 
% High Quality SimRank-Based Similarity Search

% The input arguments are:
%
%  - W,      n-by-n matrix;
%  - v,      n-by-1 vector;
%  - c,      scalar; 
%  - kmax,   summation item number;


x=cell(kmax+1,1);
x{1}=v; %n*1 vector where j-th element is 1
for t=1:kmax
   x{t+1}=A*x{t};
end

D=DiagEst_update_3(A,c,kmax);

y=cell(kmax+1,1);
y{1}=D(:,1).*x{kmax+1};
for t=1:kmax
    y{t+1}=D(:,t+1).*x{kmax+1-t}+c*(y{t}'*A)';
end

s_k_j=y{kmax+1};


