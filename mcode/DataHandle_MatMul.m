function [A_norm] = DataHandle_MatMul(fn)
%DataHandle function convert .txt file to matrix and column-normalised A
%If index from 0,then index add 1
%
%Input: fn -- .txt file
%Output: A -- column normalised matrix 
A=readmatrix(fn);
n1=max(max(A(:,1)),max(A(:,2)));
try
   A=sparse(A(:,1),A(:,2),1,n1,n1);
catch MException
    disp(MException);
    A=sparse(A(:,1)+1,A(:,2)+1,1,n1+1,n1+1);
end

n=size(A,1);
d=sum(A,1)';
inv_d = 1./d;
inv_d(~isfinite(inv_d)) = 0;
A_norm = A*spdiags(inv_d,0,n,n);

% A_norm = A;
% col_sum = sum(A);
% idx_nonzero_sum = find(col_sum);
% for i = idx_nonzero_sum
%     A_norm(:, i) = A(:, i) / col_sum(1, i);
% end


