function D=DiagEst_update_3(A,c,k)

%DiagEst compute main diagonal estimator
%Input: A -- column normalised adjacency matrix
%       c -- decay factor
%       k -- numbers of iteration
%       V -- row orthogonal matrix
%Output:diagonal matrix sequences

n=size(A,2);

D=zeros(n,k+1);
one_vec_n=ones(n,1);
D(:,1)=one_vec_n;

for j=1:k
    nume=zeros(n,1);
    for i=1:n
        v_i=sparse(i,1,1,n,1);
        u_vec_cell=cell(j+1,1);
        u_vec_cell{1}=v_i;
        for jj=1:j
            u_vec_cell{jj+1}=A*u_vec_cell{jj};
        end
        
        w_vec_cell=cell(j+1,1);
        w_vec_cell{1}=D(:,1).*u_vec_cell{j+1};
        for ll=2:j
            w_vec_cell{ll}=D(:,ll).*u_vec_cell{j+2-ll} + c*(w_vec_cell{ll-1}'*A)';
        end
        w_vec_cell{j+1}=c*(w_vec_cell{j}'*A)';
        nume(i)=w_vec_cell{j+1}(i);
    end
    D(:,j+1)=one_vec_n - nume;
end


