function [ft , H]=get_MPC_Matrices(A,B,Q,R,F,P,k,x_k,ref,U_K)
n=size(A,1);
p=size(B,2);

%%%%%%%%%%%%

% define M, E and C
M=[eye(n);zeros(P*n, n)];

E = zeros((P+1)*n, p);

C=zeros((P+1)*n, P*p);

tmp=eye(n); 
tmp2=tril(ones(P*p));

%¡¡Update M and C
for i=1:P 
    rows =i*n+(1:n); 
    C(rows,:)=[tmp*B, C(rows-n, 1:end-p)]; 
    E(rows,:)=tmp*B+E(rows-n,:);
    tmp= A*tmp; 
    M(rows,:)=tmp; 
end

C=C*tmp2;

Q_bar = kron(eye(P), Q);
Q_bar = blkdiag(Q_bar, F);
R_bar = kron(eye(P), R);

if k==1
    E=M*x_k-ref;
else
    E=M*x_k+E*U_K(:, k-1)-ref;
end

H=2*(C'*Q_bar*C+R_bar);
H=(H+H')/2;  % Ensure that it is a symmetric matrix
ft=2*C'*Q_bar*E;

end