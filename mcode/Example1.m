function Example1()

% test the ApproDiag algorithm

A = [0 0 1/2 0 0 0;1 0 0 1 0 0;0 1/2 0 0 0 0;0 0 1/2 0 0 1;0 1/2 0 0 0 0;0 0 0 0 1 0];
W2=[1 0;0 1;0 -1;0 1;0 -1;0 1];
W3=[1 0 0;0 1 0;0 0 1;0 0 -1;0 0 1;0 0 -1];
W4=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1;0 0 0 -1;0 0 0 1];
W5=[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1;0 0 0 0 -1];
W6=[1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];


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
