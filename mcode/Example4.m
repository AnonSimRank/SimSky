function Example4()
% test SimSky.m 

A=[0 0 1/2 0 0 0;1 0 0 1 0 0;0 1/2 0 0 0 0;0 0 1/2 0 0 1;0 1/2 0 0 0 0;0 0 0 0 1 0];

n=size(A,1);
n
s=n;
c=0.8;
k=2;
m1=3;
m2=3;
j=1;
v=sparse(j,1,1,n,1);

W=zeros(n,s);
W(1:s,1:s)=speye(s,s);
x=ones(n-s,1);
[row,col,val]=find(x);
x(mod(row,2)==1)=-1;
W(s+1:n,s)=x;

baseline_ss = single_source_simrank(A, v, c, k);
[soar_ss,Q,P,H,m2] = SimSky(A,W,m1,m2,k,c,v,1);
norm(c*A'*Q(:,1:m2) + P(:,1:m2) - Q*H)
norm(baseline_ss - soar_ss)

