%Compare 30 S-shaped trajectories with the other four observers and output error box plots

clear ;
close all ;
clc;
%%Define total step size N
N=1600;  % Nmax=1600
%%Define prediction step size P
P=5;
%% Initial parameter definition
self.r=0.05;  % wheel radius (m)
self.m=8;  % Platform mass (kg)
self.a=0.35;  % Half the width of the platform (m)
self.b=0.25;  % Half of the platform length (m)
self.Jz=0.1;  % Inertia moment of the platform around the center of rotation (N · m)
self.Jw=0.2;  % Moment of inertia of the wheel around the center of rotation (N · m)
self.T=0.1;  % sample time
self.D=0.2;  % Wheel viscous friction coefficient
theta=0;  % Deviation angle between world coordinate q and moving coordinate r
w=zeros(4,N);  % The speed of four wheels
%f=10*[sign(w(1));sign(w(2));sign(w(3));sign(w(4))];  % 轮子静摩擦

REF=zeros(6,1600);

sum_perror1 = zeros(3,30);
sum_perror2 = zeros(3,30);
sum_perror3 = zeros(3,30);
sum_perror4 = zeros(3,30);
sum_perror5 = zeros(3,30);

sum_ferror1 = zeros(3,30);
sum_ferror2 = zeros(3,30);
sum_ferror3 = zeros(3,30);
sum_ferror4 = zeros(3,30);
sum_ferror5 = zeros(3,30);

for t = 1:30

X_K=zeros(6, N);
U_K=zeros(4, N);
n = size(X_K, 1);
p = size(U_K, 1);

Aeq=zeros(P*p,P*p);
beq=zeros(P*p,1);

Q = diag([100,100,100,50,50,50]);

F = diag([100,100,100,50,50,50]);

R = 1*eye(p);

% Parameter Definition of Sliding Mode Observer
lambda1 = zeros(3,N);
lambda2 = zeros(3,N);
    
dK111=zeros(1,N);
dK112=zeros(1,N);
dK113=zeros(1,N);
K111=zeros(1,N);
K112=zeros(1,N);
K113=zeros(1,N);

dK121=zeros(1,N);
dK122=zeros(1,N);
dK123=zeros(1,N);
K121=zeros(1,N);
K122=zeros(1,N);
K123=zeros(1,N);

dK131=zeros(1,N);
dK132=zeros(1,N);
dK133=zeros(1,N);
K131=zeros(1,N);
K132=zeros(1,N);
K133=zeros(1,N);

dK141=zeros(1,N);
dK142=zeros(1,N);
dK143=zeros(1,N);
K141=zeros(1,N);
K142=zeros(1,N);
K143=zeros(1,N);

K211 = 1;
K212 = 1;
K213 = 1;

K221 = 0.1;
K222 = 0.1;
K223 = 0.1;

K231 = 3;
K232 = 3;
K233 = 3;

K241 = 0.05;
K242 = 0.05;
K243 = 0.05;



K31=3*[1 0 0; 0 1 0; 0 0 1];
K32=1*[1 0 0; 0 1 0; 0 0 1];
K33=3*[1 0 0; 0 1 0; 0 0 1];
K34=1*[1 0 0; 0 1 0; 0 0 1];

% Define sliding mode observer output
f = zeros(3, N-P+1);  % disturbance
f_obs = zeros(3, N-P+1);  % Disturbance observations
s = zeros(3, N-P+1);  % sliding surface
p_obs = zeros(3, N-P+1);  % Position observations

    % Comparison Observer1：
    f_obs1 = zeros(3, N);
    s1 = zeros(3, N);
    p_obs1 = zeros(3, N);
    e1 = zeros(3,N);
    
    m21=2;
    n21=1;
    alpha_g1=0.1;
    alpha_g2=0.28;
    alpha_g3=0.1;
    
    % Comparison Observer2：
    p_obs2 = zeros(3, N);
    x2g = zeros(3, N);
    f_obs2 = zeros(3, N);
    e2 = zeros(3,N);
    epsilon1=0.1;
    epsilon2=0.1;
    epsilon3=0.1;
    
    % Comparison Observer3：
    p_obs3 = zeros(3, N);
    f_obs3 = zeros(3, N);
    s3 = zeros(3,N);
    
    % Comparison Observer4:
    p_obs4 = zeros(3, N);
    f_obs4 = zeros(3, N);
    s4 = zeros(3,N);


% Initial Position
REF(1:2,1) = [5;3];
X_K(1,1) = 4-0.01*t;
X_K(2,1) = 2-0.01*t;
p_obs(:,1) = X_K(1:3,1);
p_obs1(:,1) = X_K(1:3,1);
p_obs2(:,1) = X_K(1:3,1);
p_obs3(:,1) = X_K(1:3,1);
p_obs4(:,1) = X_K(1:3,1);

X_K1=X_K; 
U_K1=U_K;
w1=w;
X_K2=X_K;
U_K2=U_K;
w2=w;
X_K3=X_K;
U_K3=U_K;
w3=w;
X_K4=X_K;
U_K4=U_K;
w4=w;
w_leader=w;

perror1 = [0;0;0];
perror2 = [0;0;0];
perror3 = [0;0;0];
perror4 = [0;0;0];
perror5 = [0;0;0];

ferror1 = [0;0;0];
ferror2 = [0;0;0];
ferror3 = [0;0;0];
ferror4 = [0;0;0];
ferror5 = [0;0;0];

k=1;
while k <= N

    
    if k<800
        w_leader(1,k)=1.2;w_leader(2,k)=1.3;w_leader(3,k)=1;w_leader(4,k)=1.5;
    else
        w_leader(1,k)=1.5;w_leader(2,k)=1;w_leader(3,k)=1.3;w_leader(4,k)=1.2;
    end
    
    [~, ~, ~, J_p_leader, R_fi_leader]=get_model(REF(3,k), self, REF(6,k));
    REF(4:6,k) = self.r*R_fi_leader*J_p_leader*w_leader(1:4,k);
    REF(1:3,k+1) = REF(1:3,k)+REF(4:6,k)*self.T;
    
t3=(k-1)*0.1;

f(:,k) = 0.01*t*[0.5+0.5*sin(0.5*t3) ;
    1*sin(0.07*t3)+0.5*sin(0.05*t3) 
    sin(0.1*t3)+0.5*cos(0.5*t3) ];
     
    [A, B, M, J_p, R_fi]=get_model(X_K(3,k), self, X_K(6,k));
    
     %% Obtain the trajectory of the leader
    % Predicting the trajectory of leader
    ref=zeros(n,P+1);
    ref(:,1)=REF(:,k);
    % Default leader speed remains unchanged
    ref(4,:)=ref(4,1);
    ref(5,:)=ref(5,1);
    ref(6,:)=ref(6,1);
    % predicted position
    for i=1:P
        %s(k)=s(k-1)+v*t
        ref(1,i+1)=ref(1,i)+ref(4,i)*self.T;
        ref(2,i+1)=ref(2,i)+ref(5,i)*self.T;
        ref(3,i+1)=ref(3,i)+ref(6,i)*self.T;
    end
    ref=reshape(ref,n*(P+1),1);

    [ft , H] = get_MPC_Matrices(A,B,Q,R,F,P,k,X_K(:, k),ref,U_K);
    
    umax = eye(P*p, 1);
    umax(:,:)=0.5;
    
    dU_k = quadprog(H,ft,[],[],Aeq,beq,-umax,umax);
    U_K(:,k+1) = U_K(:,k)+dU_k(1:4,1);

    delta_w=inv(M)*(U_K(:,k+1) - self.D*w(1:4,k));
    w(1:4,k+1)=w(1:4,k)+delta_w*self.T; 

    v_id = self.r*R_fi*J_p*w(1:4,k+1)-f_obs(:,k);
    
    v_real = v_id+f(:,k);
    X_K(4:6,k+1) = v_real;
    X_K(1:3,k+1) = X_K(1:3,k) +  X_K(4:6,k+1)*self.T;
    
        % 方法1 AVPSMO
    s(:,k) = X_K(1:3,k)-p_obs(:,k);
    
    eta0 = 1000000000000000000000;  % eta0 = inf
    eta1 = 1.5;  % eta1>1
    eta2 = 0.1;  % 0<eta2<1
    eta3 = 1.5;  % eta3>1
    
    lambda1(:,k) = [eta1*tanh(abs(s(1,k))^eta0)*abs(s(1,k))+eta2;
              eta1*tanh(abs(s(2,k))^eta0)*abs(s(2,k))+eta2;
              eta1*tanh(abs(s(3,k))^eta0)*abs(s(3,k))+eta2];
    lambda2(:,k) = [0.5*eta3+0.5+(0.5*eta3-0.5)*sign(abs(s(1,k))-1);
              0.5*eta3+0.5+(0.5*eta3-0.5)*sign(abs(s(2,k))-1);
              0.5*eta3+0.5+(0.5*eta3-0.5)*sign(abs(s(3,k))-1)];
    
    dK111(:,k)=0.02*s(1,k)*s(1,k);
    dK112(:,k)=0.02*s(2,k)*s(2,k);
    dK113(:,k)=0.02*s(3,k)*s(3,k);
    K111(:,1)=3;K112(:,1)=3;K113(:,1)=3;
    K111(:,k+1)=K111(:,k)+dK111(:,k)*self.T;
    K112(:,k+1)=K112(:,k)+dK112(:,k)*self.T;
    K113(:,k+1)=K113(:,k)+dK113(:,k)*self.T;

    dK121(:,k)=0.2*s(1,k)*sat(s(1,k));
    dK122(:,k)=0.2*s(2,k)*sat(s(2,k));
    dK123(:,k)=0.2*s(3,k)*sat(s(3,k));
    K121(:,1)=2;K122(:,1)=2;K123(:,1)=2;
    K121(:,k+1)=K121(:,k)+dK121(:,k)*self.T;
    K122(:,k+1)=K122(:,k)+dK122(:,k)*self.T;
    K123(:,k+1)=K123(:,k)+dK123(:,k)*self.T;

    dK131(:,k)=0.2*abs(s(1,k))^(lambda1(1,k));
    dK132(:,k)=0.2*abs(s(2,k))^(lambda1(2,k));
    dK133(:,k)=0.2*abs(s(3,k))^(lambda1(3,k));
    K131(:,1)=3;K132(:,1)=3;K133(:,1)=3;
    K131(:,k+1)=K131(:,k)+dK131(:,k)*self.T;
    K132(:,k+1)=K132(:,k)+dK132(:,k)*self.T;
    K133(:,k+1)=K133(:,k)+dK133(:,k)*self.T;

    dK141(:,k)=0.1*abs(s(1,k))^(lambda2(1,k));
    dK142(:,k)=0.1*abs(s(2,k))^(lambda2(2,k));
    dK143(:,k)=0.1*abs(s(3,k))^(lambda2(3,k));
    K141(:,1)=2;K142(:,1)=2;K143(:,1)=2;
    K141(:,k+1)=K141(:,k)+dK141(:,k)*self.T;
    K142(:,k+1)=K142(:,k)+dK142(:,k)*self.T;
    K143(:,k+1)=K143(:,k)+dK143(:,k)*self.T;
    
    K11 = diag([K111(:,k),K112(:,k),K113(:,k)]);
    K12 = diag([K121(:,k),K122(:,k),K123(:,k)]);
    K13 = diag([K131(:,k),K132(:,k),K133(:,k)]);
    K14 = diag([K141(:,k),K142(:,k),K143(:,k)]);
    
    s_lam1 = diag([abs(s(1,k))^lambda1(1,k),abs(s(2,k))^lambda1(2,k),abs(s(3,k))^lambda1(3,k)]);
    s_lam2 = diag([abs(s(1,k))^lambda2(1,k),abs(s(2,k))^lambda2(2,k),abs(s(3,k))^lambda2(3,k)]);
    
    dp_obs = self.r*R_fi*J_p*w(1:4,k+1);
    p_obs(:,k+1) = p_obs(:,k)+dp_obs*self.T;
    
        theta=X_K(3,k+1);
    
    R_fi=[cos(theta),-sin(theta),0;
    sin(theta), cos(theta), 0;
    0,0,1];
    
    
    s(:,k+1) = X_K(1:3,k+1)-p_obs(:,k+1);
    
    lambda1(:,k+1) = [eta1*tanh(abs(s(1,k+1))^eta0)*abs(s(1,k+1))+eta2;
              eta1*tanh(abs(s(2,k+1))^eta0)*abs(s(2,k+1))+eta2;
              eta1*tanh(abs(s(3,k+1))^eta0)*abs(s(3,k+1))+eta2];
    lambda2(:,k+1) = [0.5*eta3+0.5+(0.5*eta3-0.5)*sign(abs(s(1,k+1))-1);
              0.5*eta3+0.5+(0.5*eta3-0.5)*sign(abs(s(2,k+1))-1);
              0.5*eta3+0.5+(0.5*eta3-0.5)*sign(abs(s(3,k+1))-1)];
          
    s_lam1 = diag([abs(s(1,k+1))^lambda1(1,k+1),abs(s(2,k+1))^lambda1(2,k+1),abs(s(3,k+1))^lambda1(3,k+1)]);
    s_lam2 = diag([abs(s(1,k+1))^lambda2(1,k+1),abs(s(2,k+1))^lambda2(2,k+1),abs(s(3,k+1))^lambda2(3,k+1)]);
    K11 = diag([K111(:,k+1),K112(:,k+1),K113(:,k+1)]);
    K12 = diag([K121(:,k+1),K122(:,k+1),K123(:,k+1)]);
    K13 = diag([K131(:,k+1),K132(:,k+1),K133(:,k+1)]);
    K14 = diag([K141(:,k+1),K142(:,k+1),K143(:,k+1)]);
    
    f_obs(:,k+1) = K11*s(:,k+1)+(K12+K13*s_lam1+K14*s_lam2)*sat(s(:,k+1));
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Comparative simulation1：ATSMDO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [A1, B1, M1, J_p1, R_fi1]=get_model(X_K1(3,k), self, X_K1(6,k));

    ref=zeros(n,P+1);
    ref(:,1)=REF(:,k);

    ref(4,:)=ref(4,1);
    ref(5,:)=ref(5,1);
    ref(6,:)=ref(6,1);

    for i=1:P
        %s(k)=s(k-1)+v*t
        ref(1,i+1)=ref(1,i)+ref(4,i)*self.T;
        ref(2,i+1)=ref(2,i)+ref(5,i)*self.T;
        ref(3,i+1)=ref(3,i)+ref(6,i)*self.T;
    end
    ref=reshape(ref,n*(P+1),1);
    

    [ft1 , H1] = get_MPC_Matrices(A1,B1,Q,R,F,P,k,X_K1(:, k),ref,U_K1);
    
    umax = eye(P*p, 1);
    umax(:,:)=0.5;
    
    dU_k1 = quadprog(H1,ft1,[],[],Aeq,beq,-umax,umax);
    U_K1(:,k+1) = U_K1(:,k)+dU_k1(1:4,1);

    delta_w=inv(M1)*(U_K1(:,k+1) - self.D*w1(1:4,k));
    w1(1:4,k+1)=w1(1:4,k)+delta_w*self.T;  

    v_id1 = self.r*R_fi1*J_p1*w1(1:4,k+1)-f_obs1(:,k);
    v_real1 = self.r*R_fi1*J_p1*w1(1:4,k+1)-f_obs1(:,k)+f(:,k);
    X_K1(4:6,k+1) = v_real1;
    X_K1(1:3,k+1) = X_K1(1:3,k) +  X_K1(4:6,k+1)*self.T;
   %% 对比观测器1 ATSMDO
    dp_obs1 = self.r*R_fi1*J_p1*w1(1:4,k+1);
    p_obs1(:,k+1) = p_obs1(:,k)+dp_obs1*self.T;
   
    e1(:,k+1) = X_K1(1:3,k+1)-p_obs1(:,k+1);
    dot_e1 = (e1(:,k+1)-e1(:,k))/self.T;
    
    de1mn1 = dot_e1(1)^(m21/n21);
    de1mn2 = dot_e1(2)^(m21/n21);
    de1mn3 = dot_e1(3)^(m21/n21);
    
    s1(1,k+1) = e1(1,k+1)+K211*dot_e1(1)+K221*de1mn1;
    s1(2,k+1) = e1(2,k+1)+K212*dot_e1(2)+K222*de1mn2;
    s1(3,k+1) = e1(3,k+1)+K213*dot_e1(3)+K223*de1mn3;
    
    de1mn1 = dot_e1(1)^(m21/n21-1);
    de1mn2 = dot_e1(2)^(m21/n21-1);
    de1mn3 = dot_e1(3)^(m21/n21-1);

    dot_f_obs11 = (dot_e1(1)+K231*s1(1,k+1)+K241*sign(s1(1,k+1)))/(K211+K221*m21*de1mn1/n21)+alpha_g1;
    dot_f_obs12 = (dot_e1(2)+K232*s1(2,k+1)+K242*sign(s1(2,k+1)))/(K212+K222*m21*de1mn2/n21)+alpha_g2;
    dot_f_obs13 = (dot_e1(3)+K233*s1(3,k+1)+K243*sign(s1(3,k+1)))/(K213+K223*m21*de1mn3/n21)+alpha_g3;
    
    f_obs1(:,k+1) = f_obs1(:,k)+[dot_f_obs11;dot_f_obs12;dot_f_obs13]*self.T;
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Comparative simulation2：ESO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [A2, B2, M2, J_p2, R_fi2]=get_model(X_K2(3,k), self, X_K2(6,k));
    
    ref=zeros(n,P+1);
    ref(:,1)=REF(:,k);

    ref(4,:)=ref(4,1);
    ref(5,:)=ref(5,1);
    ref(6,:)=ref(6,1);

    for i=1:P
        %s(k)=s(k-1)+v*t
        ref(1,i+1)=ref(1,i)+ref(4,i)*self.T;
        ref(2,i+1)=ref(2,i)+ref(5,i)*self.T;
        ref(3,i+1)=ref(3,i)+ref(6,i)*self.T;
    end
    ref=reshape(ref,n*(P+1),1);
    
    [ft2 , H2] = get_MPC_Matrices(A2,B2,Q,R,F,P,k,X_K2(:, k),ref,U_K2);
    
    umax = eye(P*p, 1);
    umax(:,:)=0.5;
    
    dU_k2 = quadprog(H2,ft2,[],[],Aeq,beq,-umax,umax);
    U_K2(:,k+1) = U_K2(:,k)+dU_k2(1:4,1);

    delta_w2=inv(M2)*(U_K2(:,k+1) - self.D*w2(1:4,k));
    w2(1:4,k+1)=w2(1:4,k)+delta_w2*self.T; 

    v_id2 = self.r*R_fi2*J_p2*w2(1:4,k+1)-f_obs2(:,k); 

    v_real2 = v_id2+f(:,k);
    X_K2(4:6,k+1) = v_real2;
    X_K2(1:3,k+1) = X_K2(1:3,k) +  X_K2(4:6,k+1)*self.T;
    
    theta2=X_K2(3,k+1);
    
    R_fi2=[cos(theta2),-sin(theta2),0;
        sin(theta2), cos(theta2), 0;
        0,0,1];
    
    %% 对比观测器2 ESO
%     dzeta1 = v_id2+e2(:,k)/epsilon;
    
    dzeta1 = v_id2(1)+e2(1,k)/epsilon1;
    dzeta2 = v_id2(2)+e2(2,k)/epsilon2;
    dzeta3 = v_id2(3)+e2(3,k)/epsilon3;
    
    dzeta = [dzeta1;dzeta2;dzeta3];
    p_obs2(:,k+1) = p_obs2(:,k) + dzeta*self.T;
    e2(:,k+1) = X_K2(1:3,k+1)-p_obs2(:,k+1);
    
    x2g(1,k+1) = e2(1,k+1)/epsilon1;
    x2g(2,k+1) = e2(2,k+1)/epsilon2;
    x2g(3,k+1) = e2(3,k+1)/epsilon3;
    %x1g(:,k+1) = x1g(:,k)+x2g(:,k+1)*self.T;
    f_obs2(:,k+1) =  x2g(:,k+1);
%% %%%%%%%%%%%%%%%%%%%%%%%%Comparative simulation3：DPSMO%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A3, B3, M3, J_p3, R_fi3]=get_model(X_K3(3,k), self, X_K3(6,k));

    ref=zeros(n,P+1);
    ref(:,1)=REF(:,k);

    ref(4,:)=ref(4,1);
    ref(5,:)=ref(5,1);
    ref(6,:)=ref(6,1);

    for i=1:P
        %s(k)=s(k-1)+v*t
        ref(1,i+1)=ref(1,i)+ref(4,i)*self.T;
        ref(2,i+1)=ref(2,i)+ref(5,i)*self.T;
        ref(3,i+1)=ref(3,i)+ref(6,i)*self.T;
    end
    ref=reshape(ref,n*(P+1),1);
    
    [ft3 , H3] = get_MPC_Matrices(A3,B3,Q,R,F,P,k,X_K3(:, k),ref,U_K3);
    
    umax = eye(P*p, 1);
    umax(:,:)=0.5;
    
    dU_k3 = quadprog(H3,ft3,[],[],Aeq,beq,-umax,umax);
    U_K3(:,k+1) = U_K3(:,k)+dU_k3(1:4,1);

    delta_w3=inv(M3)*(U_K3(:,k+1) - self.D*w3(1:4,k));
    w3(1:4,k+1)=w3(1:4,k)+delta_w3*self.T; 

    v_id3 = self.r*R_fi3*J_p3*w3(1:4,k+1)-f_obs3(:,k); 
    v_real3 = v_id3+f(:,k);
    X_K3(4:6,k+1) = v_real3;
    X_K3(1:3,k+1) = X_K3(1:3,k) +  X_K3(4:6,k+1)*self.T;
    
    %% 对比观测器3 DPSMO
    s3(:,k) = X_K3(1:3,k)-p_obs3(:,k);
    
%     dp_obs = v_id+K11*s(:,k)+(K12+K13*s_lam1+K14*s_lam2)*sat(s(:,k));
    dp_obs3 = self.r*R_fi3*J_p3*w3(1:4,k+1);
    p_obs3(:,k+1) = p_obs3(:,k)+dp_obs3*self.T;
    
        theta3=X_K3(3,k+1);
    
    R_fi3=[cos(theta3),-sin(theta3),0;
    sin(theta3), cos(theta3), 0;
    0,0,1];
    
    
    s3(:,k+1) = X_K3(1:3,k+1)-p_obs3(:,k+1);
    
    aaaAA=[abs(s3(1,k+1))^1.5 0 0;
       0 abs(s3(2,k+1))^1.5 0;
       0 0 abs(s3(3,k+1))^1.5];
    aaaBB=[abs(s3(1,k+1))^0.5 0 0;
       0 abs(s3(2,k+1))^0.5 0;
       0 0 abs(s3(3,k+1))^0.5];
    
   f_obs3(:,k+1) = K31*s3(:,k+1)+(K32+K33*aaaAA+K34*aaaBB)*sat(s3(:,k+1));
   
   %% %%%%%%%%%%%%%%%%%%%%%%%%Comparative simulation4: SMDO%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A4, B4, M4, J_p4, R_fi4]=get_model(X_K4(3,k), self, X_K4(6,k));


    ref=zeros(n,P+1);
    ref(:,1)=REF(:,k);

    ref(4,:)=ref(4,1);
    ref(5,:)=ref(5,1);
    ref(6,:)=ref(6,1);

    for i=1:P
        %s(k)=s(k-1)+v*t
        ref(1,i+1)=ref(1,i)+ref(4,i)*self.T;
        ref(2,i+1)=ref(2,i)+ref(5,i)*self.T;
        ref(3,i+1)=ref(3,i)+ref(6,i)*self.T;
    end
    ref=reshape(ref,n*(P+1),1);
    
    [ft4 , H4] = get_MPC_Matrices(A4,B4,Q,R,F,P,k,X_K4(:, k),ref,U_K4);
    
    umax = eye(P*p, 1);
    umax(:,:)=0.5;
    
    dU_k4 = quadprog(H4,ft4,[],[],Aeq,beq,-umax,umax);
    U_K4(:,k+1) = U_K4(:,k)+dU_k4(1:4,1);

    delta_w4=inv(M4)*(U_K4(:,k+1) - self.D*w4(1:4,k));
    w4(1:4,k+1)=w4(1:4,k)+delta_w4*self.T;  

    v_id4 = self.r*R_fi4*J_p4*w4(1:4,k+1)-f_obs4(:,k); 
    v_real4 = v_id4+f(:,k);
    X_K4(4:6,k+1) = v_real4;
    X_K4(1:3,k+1) = X_K4(1:3,k) +  X_K4(4:6,k+1)*self.T;
    
    %% 对比观测器4 SMDO
    s4(:,k) = X_K4(1:3,k)-p_obs4(:,k);
    
%     dp_obs = v_id+K11*s(:,k)+(K12+K13*s_lam1+K14*s_lam2)*sat(s(:,k));
    dp_obs4 = self.r*R_fi4*J_p4*w4(1:4,k+1);
    p_obs4(:,k+1) = p_obs4(:,k)+dp_obs4*self.T;
    
        theta4=X_K4(3,k+1);
    
    R_fi4=[cos(theta4),-sin(theta4),0;
    sin(theta4), cos(theta4), 0;
    0,0,1];
    
    
    s4(:,k+1) = X_K4(1:3,k+1)-p_obs4(:,k+1);
    
    aaaAA=[abs(s4(1,k+1))^1.5 0 0;
       0 abs(s4(2,k+1))^1.5 0;
       0 0 abs(s4(3,k+1))^1.5];
    aaaBB=[abs(s4(1,k+1))^0.5 0 0;
       0 abs(s4(2,k+1))^0.5 0;
       0 0 abs(s4(3,k+1))^0.5];
    
    f_obs4(:,k+1) = K31*s4(:,k+1)+K32*sat(s4(:,k+1));
    %% next circle
    perror1 = perror1+abs(X_K(1:3,k)-REF(1:3,k))*self.T;
    perror2 = perror2+abs(X_K1(1:3,k)-REF(1:3,k))*self.T;
    perror3 = perror3+abs(X_K2(1:3,k)-REF(1:3,k))*self.T;
    perror4 = perror4+abs(X_K3(1:3,k)-REF(1:3,k))*self.T;
    perror5 = perror5+abs(X_K4(1:3,k)-REF(1:3,k))*self.T;
    
    ferror1 = ferror1+abs(f(1:3,k)-f_obs(1:3,k))*self.T;
    ferror2 = ferror2+abs(f(1:3,k)-f_obs1(1:3,k))*self.T;
    ferror3 = ferror3+abs(f(1:3,k)-f_obs2(1:3,k))*self.T;
    ferror4 = ferror4+abs(f(1:3,k)-f_obs3(1:3,k))*self.T;
    ferror5 = ferror5+abs(f(1:3,k)-f_obs4(1:3,k))*self.T;
    
    disp([t,k]);
    k=k+1;
end
sum_perror1(:,t) = perror1;
sum_perror2(:,t) = perror2;
sum_perror3(:,t) = perror3;
sum_perror4(:,t) = perror4;
sum_perror5(:,t) = perror5;

sum_ferror1(:,t) = ferror1;
sum_ferror2(:,t) = ferror2;
sum_ferror3(:,t) = ferror3;
sum_ferror4(:,t) = ferror4;
sum_ferror5(:,t) = ferror5;
end

%% draw pictures

set(0, 'DefaulttextInterpreter', 'latex')
colors = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880]];
figure;
hold on;
subplot(1,3,1)
boxplot([sum_perror1(1,:)', sum_perror2(1,:)', sum_perror3(1,:)', sum_perror4(1,:)', sum_perror5(1,:)'], 'Colors',colors, 'Symbol','+');
%legend("AVPSMO-MPC","ATSMDO-MPC","ESO-MPC","DPSMO-MPC","SMDO-MPC");
% ylabel('IAE','FontSize',16);
xlabel('$e_{\rm x}$(m)','FontSize',16);
subplot(1,3,2)
boxplot([sum_perror1(2,:)', sum_perror2(2,:)', sum_perror3(2,:)', sum_perror4(2,:)', sum_perror5(2,:)'], 'Colors',colors, 'Symbol','+');
xlabel('$e_{\rm y}$(m)','FontSize',16);
subplot(1,3,3)
boxplot([sum_perror1(3,:)', sum_perror2(3,:)', sum_perror3(3,:)', sum_perror4(3,:)', sum_perror5(3,:)'], 'Colors',colors, 'Symbol','+');
xlabel('$e_\varphi$(m)','FontSize',16);
hLegend = legend(findall(gca,'Tag','Box'), {'AVPSMO-MPC','ATSMDO-MPC','NESO-MPC','DPSMO-MPC','SMDO-MPC'},'Location','NorthOutside','NumColumns',3);

figure;
hold on;
subplot(1,3,1)
boxplot([sum_ferror1(1,:)', sum_ferror2(1,:)', sum_ferror3(1,:)', sum_ferror4(1,:)', sum_ferror5(1,:)'], 'Colors',colors, 'Symbol','+');
%legend("AVPSMO-MPC","ATSMDO-MPC","ESO-MPC","DPSMO-MPC","SMDO-MPC");
% ylabel('IAE','FontSize',16);
xlabel('$f_{\rm x}$','FontSize',16);
subplot(1,3,2)
boxplot([sum_ferror1(2,:)', sum_ferror2(2,:)', sum_ferror3(2,:)', sum_ferror4(2,:)', sum_ferror5(2,:)'], 'Colors',colors, 'Symbol','+');
xlabel('$f_{\rm y}$','FontSize',16);
subplot(1,3,3)
boxplot([sum_ferror1(3,:)', sum_ferror2(3,:)', sum_ferror3(3,:)', sum_ferror4(3,:)', sum_ferror5(3,:)'], 'Colors',colors, 'Symbol','+');
xlabel('$f_\omega$','FontSize',16);
hLegend = legend(findall(gca,'Tag','Box'), {'AVPSMO','ATSMDO','NESO','DPSMO','SMDO'},'NumColumns',5, 'Location', 'northoutside');