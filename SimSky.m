function [soar_ss,Q,P,H,m2] = SimSky(A,W,m1,m2,k,c,query_vec,reorthog)

D=ApproDiag(A,c,k,W);

n=size(A,1);
final_vec_seq=zeros((k+1)*n,1);
[U,T,m1] = arnoldi(A,query_vec,m1,reorthog);
T_=T(1:m1,1:m1);
e_m_=sparse(1,1,1,m1,1);
beta=norm(T*T_^(k-1)*e_m_);
if beta < eps
    error('norm of the first element of vec_seq is %d', beta);
end
final_vec_seq(1:n)=D(:,1).* (U*T*T_^(k-1)*e_m_./beta);
final_vec_seq(n+1:2*n)=D(:,k+1).*(query_vec./beta);
for i=2:k
    final_vec_seq(i*n+1:(i+1)*n)=D(:,k+2-i).*(U*T*T_^(i-2)*e_m_./beta);
end

Q=zeros(n,m2+1);
P=zeros(n,m2+1);
H=zeros(m2+1,m2);
B=c*A';

V=zeros((k+1)*n,m2+1);
V(:,1)=final_vec_seq;
tol=n*eps;
for i=1:m2
    r=B*V(1:n,i)+V(k*n+1:(k+1)*n,i);
    t=V(1:k*n,i);
    s=[r;t];
    ow=norm(r);
    for j=1:i
        temp=r'*V(1:n,j);
        s=s-temp*V(:,j);
        H(j,i)=temp;
    end
    
    if reorthog==1
        for j=1:i
            temp=s(1:n)'*V(1:n,j);
            H(j,i)=H(j,i)+temp;
            s=s-temp*V(:,j);
        end
    end
    
    H(i+1,i)=norm(s(1:n));
    if H(i+1,i) < tol*ow
        m=i;
        H=H(1:1+m,1:m);
        V=V(:,1:1+m);
        Q=V(1:n,1:1+m);
        P=V(k*n+1:(k+1)*n,1:1+m);
    end
    V(:,i+1)=s/H(i+1,i);
end

Q=V(1:n,:);
P=V(k*n+1:(k+1)*n,:);

H1=H(1:m2,1:m2);
e_m=sparse(1,1,1,m2,1);
soar_ss=beta*Q*H*H1^(k-1)*e_m;











