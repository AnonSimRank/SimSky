function [final_vec_seq,beta]=second_order_vec_seq(W,b,k,m,reorthog)
%
% second_order_vec_seq computes W^k*b, b, W*b, W^2*b, W^3*b, ..., W^{k-1}*b via
% standard krylov subspace method 
%
% Input:  W -- n-by-n sparse matrix
%         b -- n-by-1 initial vector
%         k -- iterations
%         m -- low-dimension krylov subspace, m << n
%         reorthog -- 1 if full reorthogonalization and 0 if not
%
%Output:  final_vec_seq -- (k+1)*n-by-1 vector [W^k*b/beta b/beta W*b/beta W^2*b/beta ... W^(k-1)*b/beta]
%         beta -- scalar norm(W^k*b)

n=size(W,1);
final_vec_seq=zeros((k+1)*n,1);
[Q,H,m_] = arnoldi(W,b,m,reorthog);
H_=H(1:m_,1:m_);
e_m_=sparse(1,1,1,m_,1);
beta=norm(H*H_^(k-1)*e_m_);
if beta < eps
    error('norm of the first element of vec_seq is %d', beta);
end
final_vec_seq(1:n)=Q*H*H_^(k-1)*e_m_./beta;
final_vec_seq(n+1:2*n)=b./beta;
for i=2:k
    final_vec_seq(i*n+1:(i+1)*n)=Q*H*H_^(i-2)*e_m_./beta;
end