function Example1()

A = [0 0 1/2 0 0 0;1 0 0 1 0 0;0 1/2 0 0 0 0;0 0 1/2 0 0 1;0 1/2 0 0 0 0;0 0 0 0 1 0];
W2=[1 0;0 1;0 -1;0 1;0 -1;0 1];
W3=[1 0 0;0 1 0;0 0 1;0 0 -1;0 0 1;0 0 -1];
W4=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1;0 0 0 -1;0 0 0 1];
W5=[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1;0 0 0 0 -1];
W6=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];

% W2_column=size(W2,2);
% W3_column=size(W3,2);
% W4_column=size(W4,2);
% W5_column=size(W5,2);
% W6_column=size(W6,2);
% 
% (W2.*(A*W2))*ones(W2_column,1)
% (W3.*(A*W3))*ones(W3_column,1)
% (W4.*(A*W4))*ones(W4_column,1)
% (W5.*(A*W5))*ones(W5_column,1)
% (W6.*(A*W6))*ones(W6_column,1)

c=0.8;
k=2;
D2=ApproDiag(A,c,k,W2);
D2
D3=ApproDiag(A,c,k,W3);
D3
D4=ApproDiag(A,c,k,W4);
D4
D5=ApproDiag(A,c,k,W5);
D5
D6=ApproDiag(A,c,k,W6);
D6
D=Varied_D_all(A,c,k);
D
