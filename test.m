function test()

A=[0 0 1/2 0 0 0;1 0 0 1 0 0;0 1/2 0 0 0 0;0 0 1/2 0 0 1;0 1/2 0 0 0 0;0 0 0 0 1 0];
n=size(A,1);
s=n;
v=[1 0 0 0 0 0]';
e1=[1 0 0]';
c=0.8;
m1=3;
m2=3;
k=3;
B=c*A';
W=eye(n,n);

W=zeros(n,s);
W(1:s,1:s)=speye(s,s);
x=ones(n-s,1);
[row,col,val]=find(x);
x(mod(row,2)==1)=-1;
W(s+1:n,s)=x;

baseline = single_source_simrank(A, v, c, k);
[soar_ss,Q,P,H,m2] = SimSky(A,W,m1,m2,k,c,v,1);
norm(c*A'*Q(:,1:m2) + P(:,1:m2) - Q*H)
norm(soar_ss - baseline)





