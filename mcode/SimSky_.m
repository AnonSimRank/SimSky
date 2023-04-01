function [Q,P,H,m]=superior_SOAR(B,m,k,final_vec_seq,reorthog)
%Input:   
%         B -- n-by-n matrix
%         final_vec_seq -- (k+1)*n-by-1 vector [W^k*v0/norm(W^k*v0) v0/norm(W^k*v0) W*v0/norm(W^k*v0) ... W^(k-1)*v0/norm(W^k*v0)]
%         m -- low-dimension krylov subspace, m << n
%         k -- iterations
%         reorthog -- 1 if full reorthogonalization and 0 if not
%Output: 
%         Q -- n-by-(m+1) orthogonal matrix
%         P -- n-by-(m+1) matrix
%         H -- (m+1)-by-m upper Hessenberg matrix
%         B*Q(:,1:m) + P(:,1:m) = Q*H

n=size(B,1);
Q=zeros(n,m+1);
P=zeros(n,m+1);
H=zeros(m+1,m);

V=zeros((k+1)*n,m+1);
V(:,1)=final_vec_seq;
tol=n*eps;
for i=1:m
    r=B*V(1:n,i)+V(k*n+1:(k+1)*n,i);
    t=V(1:k*n,i);
    s=[r;t];
    ow=norm(r);
    
    for j=1:i
        temp=r'*V(1:n,j);
%         r=r-temp*V(1:n,j);
%         t=t-temp*V(n+1:(m+1)*n,j);
        s=s-temp*V(:,j);
        H(j,i)=temp;
    end
    
    if reorthog==1
        for j=1:i
%             temp=r'*V(1:n,j);
            temp=s(1:n)'*V(1:n,j);
            H(j,i)=H(j,i)+temp;
%             r=r-temp*V(1:n,j);
%             t=t-temp*V(n+1:(m+1)*n,j);
            s=s-temp*V(:,j);
        end
    end
    
%     H(i+1,i)=norm(r);
    H(i+1,i)=norm(s(1:n));
    if H(i+1,i) < tol*ow
        m=i;
        H=H(1:1+m,1:m);
        V=V(:,1:1+m);
        Q=V(1:n,1:1+m);
        P=V(k*n+1:(k+1)*n,1:1+m);
        return;
    end
%     V(:,i+1)=[r/H(i+1,i);t/H(i+1,i)];
    V(:,i+1)=s/H(i+1,i);
end

Q=V(1:n,:);
P=V(k*n+1:(k+1)*n,:);























