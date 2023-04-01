function D=Varied_D_all(A,c,k)

%Varied_D_all compute diagonal vector all
%Input: A -- column normalised adjacency matrix
%       c -- decay factor
%       k -- number of iteration
%Output: diagonal vector D_k

n=size(A,1);
D=zeros(n,k+1);
D(:,1)=ones(n,1);

for i=1:k
    temp=zeros(n,1);
    P=A;
    for ll=1:i
        temp=temp + c^(ll)*(P.*P)'*D(:,i+1-ll);
        P=A*P;
    end
    total=ones(n,1)-temp;
    D(:,i+1)=total;
end

% D_k=D(:,k+1);
