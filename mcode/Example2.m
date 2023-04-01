function Example2()

% test SimSky_.m simplified version, the default diagonal correction matrix is the identity matrix.

A=[0 0 1/2 0 0 0;1 0 0 1 0 0;0 1/2 0 0 0 0;0 0 1/2 0 0 1;0 1/2 0 0 0 0;0 0 0 0 1 0];
A
n=size(A,1);
v=[1 0 0 0 0 0]';
e1=[1 0 0]';
c=0.8;
m=3;
k=2;
B=c*A';

[U, T, m] = arnoldi(A,v,m,1);
norm(A*U(:,1:m)-U*T)
r0=A*v;
norm(r0 - U*T*e1)
r1=v;
norm(r1 - U(:,1:m)*e1)
r2=A^2*v;
norm(r2 - U*T*T(1:m,1:m)*e1)
r3=B*r2 + r0;
r4=B*r3 + r1;
[final_vec_seq,beta]=second_order_vec_seq(A,v,k,m,1);
[Q,P,H,m]=SimSky_(B,m,k,final_vec_seq,1);
norm(B*Q(:,1:m)+P(:,1:m)-Q*H)
norm(r3-beta*Q*H*e1)
norm(r4-beta*Q*H*H(1:m,1:m)*e1)



