function D=DiagEst_update_5(A,c,k,V)

%DiagEst_update_5 computes diagonal correction matrix

[n, A_col_num]=size(A);
[V_row_num, V_col_num]=size(V);
denom=(V.*V)*ones(V_col_num,1);

if A_col_num ~= V_row_num
    disp("column number of A must be equal to row number");
end

D=zeros(n,k+1);
one_vec_n=ones(n,1);
D(:,1)=one_vec_n;

for j=1:k
    nume=zeros(n,1);
    for i=1:V_col_num
        v_i=V(:,i);
        u_vec_cell=cell(j+1,1);
        u_vec_cell{1}=v_i;
        for jj=1:j
            u_vec_cell{jj+1}=A*u_vec_cell{jj};
        end
        
        w_vec_cell=cell(j+1,1);
        w_vec_cell{1}=D(:,1).*u_vec_cell{j+1};
        if j==1
            w_vec_cell{j+1}=c*(w_vec_cell{j}'*A)';
        elseif j >= 2
            for ll=2:j
                w_vec_cell{ll}=D(:,ll).*u_vec_cell{j+2-ll} + c*(w_vec_cell{ll-1}'*A)';
            end
            w_vec_cell{j+1}=c*(w_vec_cell{j}'*A)';
        end
        nume = nume + v_i.*w_vec_cell{j+1};
    end
    
    D(:,j+1)=one_vec_n - nume./denom;
end



