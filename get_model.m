function [A, B, M, J_p, R_fi]=get_model(theta, self, w)
a=self.a;
b=self.b;
Jz=self.Jz;
Jw=self.Jw;

theta_=theta+pi/4;
Aj=self.m*self.r^2/8;
Bj=Jz*self.r^2/(16*(a+b));

F_p=0.25*[(2^0.5)*sin(theta_),   (2^0.5)*cos(theta_), (2^0.5)*cos(theta_),   (2^0.5)*sin(theta_);
    -(2^0.5)*cos(theta_), (2^0.5)*sin(theta_), (2^0.5)*sin(theta_), -(2^0.5)*cos(theta_);
    -1/(a+b),                        1/(a+b),              -1/(a+b),                    1/(a+b)];

F_=[   (2^0.5)*sin(theta_),-(2^0.5)*cos(theta_),   -(a+b);
     (2^0.5)*cos(theta_),(2^0.5)*sin(theta_), (a+b);
    (2^0.5)*cos(theta_),(2^0.5)*sin(theta_),  -(a+b);
     (2^0.5)*sin(theta_), -(2^0.5)*cos(theta_),  (a+b)];

dF_=w*[ (2^0.5)*cos(theta_),(2^0.5)*sin(theta_),  0;
     -(2^0.5)*sin(theta_),(2^0.5)*cos(theta_),  0;
      -(2^0.5)*sin(theta_), (2^0.5)*cos(theta_), 0;
     (2^0.5)*cos(theta_), -(2^0.5)*sin(theta_), 0];

M=[Aj+Bj+Jw,-Bj,Bj,Aj-Bj;
    -Bj,Aj+Bj+Jw,Aj-Bj,Bj;
    Bj,Aj-Bj,Aj+Bj+Jw,-Bj;
    Aj-Bj,Bj,-Bj,Aj+Bj+Jw];

J_p=0.25*[1,1,1,1;
    -1,1,1,-1;
    -1/(a+b), 1/(a+b), -1/(a+b), 1/(a+b)];

R_fi=[cos(theta),-sin(theta),0;
    sin(theta), cos(theta), 0;
    0,0,1];

%% Generate A and B matrices
A=zeros(6,6);
A(1:3,4:6)=eye(3);
tmp=-(F_p*dF_+self.D*F_p*inv(M)*F_);
A(4:6,4:6)=tmp;

B=zeros(6,4);
B(4:6,:)=self.r*F_p/M;

n = size(A, 1);

% Discretization of matrix A and B
A=eye(n)+self.T*A;
B=self.T*B;
end