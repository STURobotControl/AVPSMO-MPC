% AVPSMO as disturbance observer output
% Record the observer convergence time for 30 S-shaped routes

clear ;
close all ;
clc;

N=1600;  % Nmax=1600

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

clock_ef1 = zeros(1,30);
clock_ef2 = zeros(1,30);
clock_ef3 = zeros(1,30);
clock_ef4 = zeros(1,30);
clock_ef5 = zeros(1,30);
clock_ef6 = zeros(1,30);

%30 times
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
    

    f = zeros(3, N);  
    f_obs = zeros(3, N);  
    s = zeros(3, N);  
    p_obs = zeros(3, N);  
    
    
    % Initial Position(changes every time)
    REF(1:2,1) = [5;3];
    X_K(1,1) = 4-0.01*t;
    X_K(2,1) = 2-0.01*t;
    p_obs(:,1) = X_K(1:3,1);
    
    w_leader=w;
    
    
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
        
        % AVPSMO
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
        
        %     dp_obs = v_id+K11*s(:,k)+(K12+K13*s_lam1+K14*s_lam2)*sat(s(:,k));
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
        disp([t,k]);disp('total 30 times')
        if k < 300
            e_f1 = abs(f_obs(1,k)-f(1,k))/f(1,k);
            e_f2 = abs(f_obs(2,k)-f(2,k))/f(2,k);
            e_f3 = abs(f_obs(3,k)-f(3,k))/f(3,k);
            if e_f1<=0.05 && clock_ef1(1,t) == 0
                clock_ef1(1,t) = k;
            end
            if e_f2<=0.05 && clock_ef2(1,t) == 0
                clock_ef2(1,t) = k;
            end
            if e_f3<=0.05 && clock_ef3(1,t) == 0
                clock_ef3(1,t) = k;
            end
            if e_f1<=0.02 && clock_ef4(1,t) == 0
                clock_ef4(1,t) = k;
            end
            if e_f2<=0.02 && clock_ef5(1,t) == 0
                clock_ef5(1,t) = k;
            end
            if e_f3<=0.02 && clock_ef6(1,t) == 0
                clock_ef6(1,t) = k;
            end
        end
        k = k+1;
    end
end

figure("name","Disturbance observation");
subplot (3,1,1);
hold on;box on;
plot(0:0.1:159.9,f(1,1:N),'LineWidth',2);
plot(0:0.1:159.9,f_obs(1,1:N),'--','LineWidth',2);
xlim([0,160]);
ylim([-0.12,0.12]);
ylabel('$f_{\rm x}$','Interpreter','latex','FontSize',16);
legend("Disturbance","Estimation",'FontSize',14,'NumColumns',2,'Interpreter', 'latex')
subplot (3,1,2);
hold on;box on;
plot(0:0.1:159.9,f(2,1:N),'LineWidth',2);
plot(0:0.1:159.9,f_obs(2,1:N),'--','LineWidth',2);
xlim([0,160]);
%ylim([-0.05,0.05]);
ylabel('$f_{\rm y}$','Interpreter','latex','FontSize',16);
subplot (3,1,3);
hold on;box on;
plot(0:0.1:159.9,f(3,1:N),'LineWidth',2);
plot(0:0.1:159.9,f_obs(3,1:N),'--','LineWidth',2);
xlim([0,160]);
%ylim([-0.05,0.05]);
ylabel('$f_{\omega}$','Interpreter','latex','FontSize',16);
xlabel('Time (s)','FontSize',16,'Fontname','Times');