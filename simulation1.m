%Rectangle trajectory
%AVPSMO as disturbance observer compensation output, compared with the control effect of MPC without observer

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

    Jv=1/self.r*[1 -1 -(self.a+self.b);
        1 1 (self.a+self.b);
        1 1 -(self.a+self.b);
        1 -1 (self.a+self.b)];
%% Set reference path matrix
REF=zeros(6,1600);

%% Define state quantity and input quantity matrix
X_K=zeros(6, N);
U_K=zeros(4, N);
n = size(X_K, 1);
p = size(U_K, 1);

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

% Define sliding mode observer output
f = zeros(3, N-P+1);  % disturbance
f_obs = zeros(3, N-P+1);  % Disturbance observations
s = zeros(3, N-P+1);  % sliding surface
p_obs = zeros(3, N-P+1);  % Position observations

% Initial Position
REF(1:2,1) = [5;3];
X_K(1,1) = 4;
X_K(2,1) = 3.2;
p_obs(:,1) = X_K(1:3,1);

X_K1=X_K;
U_K1=U_K;
w1=w;
wre=w;
wre1=w;
w_leader=w;

% constraint matrix
Aeq=zeros(P*p,P*p);
beq=zeros(P*p,1);

k=1;
while k <= N

    
    if k<400
        w_leader(1,k)=2;w_leader(2,k)=2;w_leader(3,k)=2;w_leader(4,k)=2;
    elseif 400<=k && k<800
        w_leader(1,k)=2;w_leader(2,k)=-2;w_leader(3,k)=-2;w_leader(4,k)=2;
    elseif 800<=k && k<1200
        w_leader(1,k)=-2;w_leader(2,k)=-2;w_leader(3,k)=-2;w_leader(4,k)=-2;
    else 
        w_leader(1,k)=-2;w_leader(2,k)=2;w_leader(3,k)=2;w_leader(4,k)=-2;
    end
    
    [~, ~, ~, J_p_leader, R_fi_leader]=get_model(REF(3,k), self, REF(6,k));
    REF(4:6,k+1) = self.r*R_fi_leader*J_p_leader*w_leader(1:4,k);
    REF(1:3,k+1) = REF(1:3,k)+REF(4:6,k)*self.T;
    
    t3=(k-1)*0.1;

f(:,k)=[0.08*sin(0.2*t3)+0.04*cos(0.1*t3);
         0.08*cos(0.1*t3);
         0.04*sin(0.2*t3)]; 
     
    %% Obtain the model and generate A and B matrices
    [A, B, M, J_p, R_fi]=get_model(X_K(3,k), self, X_K(6,k));
    
    %%
    % Define Q matrix, n x n
     Q = diag([300,300,300,100,100,100]);
%    Q = 0.5*eye(n);

    % Define F matrix, n x n 
    F = diag([300,300,300,100,100,100]);% = 100*eye(n);
    % Define R matrix, p x p
    R = 1*eye(p);
    
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

    
    %% get H, f matrices
    [ft , H] = get_MPC_Matrices(A,B,Q,R,F,P,k,X_K(:, k),ref,U_K);
    
    %% Solving QP problems

    umax = eye(P*p, 1);
    umax(:,:)=1;
    
    dU_k = quadprog(H,ft,[],[],Aeq,beq,-umax,umax);
    U_K(:,k+1) = U_K(:,k)+dU_k(1:4,1);

    % Update status quantity
    delta_w=inv(M)*(U_K(:,k+1) - self.D*w(1:4,k));
    w(1:4,k+1)=w(1:4,k)+delta_w*self.T;  

    v_id = self.r*R_fi*J_p*w(1:4,k+1)-R_fi*f_obs(:,k);
    v_real = self.r*R_fi*J_p*w(1:4,k+1)-R_fi*f_obs(:,k)+R_fi*f(:,k);
    X_K(4:6,k+1) = v_real;
    X_K(1:3,k+1) = X_K(1:3,k) +  X_K(4:6,k+1)*self.T;
    

    wre(:,k+1) = Jv*v_real;
    
    theta=X_K(3,k+1);
    
    R_fi=[cos(theta),-sin(theta),0;
    sin(theta), cos(theta), 0;
    0,0,1];
    
        % AVPSMO
    s(:,k) = X_K(1:3,k)-p_obs(:,k);
    
    eta0 = 1000000000000000000000;  % eta0 = inf
    eta1 = 1.3;  % eta1>1
    eta2 = 0.3;  % 0<eta2<1
    eta3 = 1.6;  % eta3>1
    
    lambda1(:,k) = [eta1*tanh(abs(s(1,k))^eta0)*abs(s(1,k))+eta2;
              eta1*tanh(abs(s(2,k))^eta0)*abs(s(2,k))+eta2;
              eta1*tanh(abs(s(3,k))^eta0)*abs(s(3,k))+eta2];
    lambda2(:,k) = [0.5*eta3+0.5+(0.5*eta3-0.5)*sign(abs(s(1,k))-1);
              0.5*eta3+0.5+(0.5*eta3-0.5)*sign(abs(s(2,k))-1);
              0.5*eta3+0.5+(0.5*eta3-0.5)*sign(abs(s(3,k))-1)];
    
    dK111(:,k)=0.2*s(1,k)*s(1,k);
    dK112(:,k)=0.2*s(2,k)*s(2,k);
    dK113(:,k)=0.2*s(3,k)*s(3,k);
    K111(:,1)=3;K112(:,1)=3;K113(:,1)=3;
    K111(:,k+1)=K111(:,k)+dK111(:,k)*self.T;
    K112(:,k+1)=K112(:,k)+dK112(:,k)*self.T;
    K113(:,k+1)=K113(:,k)+dK113(:,k)*self.T;

    dK121(:,k)=0.2*s(1,k)*sat(s(1,k));
    dK122(:,k)=0.2*s(2,k)*sat(s(2,k));
    dK123(:,k)=0.2*s(3,k)*sat(s(3,k));
    K121(:,1)=1;K122(:,1)=1;K123(:,1)=1;
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
    K141(:,1)=1;K142(:,1)=1;K143(:,1)=1;
    K141(:,k+1)=K141(:,k)+dK141(:,k)*self.T;
    K142(:,k+1)=K142(:,k)+dK142(:,k)*self.T;
    K143(:,k+1)=K143(:,k)+dK143(:,k)*self.T;
    
    K11 = diag([K111(:,k),K112(:,k),K113(:,k)]);
    K12 = diag([K121(:,k),K122(:,k),K123(:,k)]);
    K13 = diag([K131(:,k),K132(:,k),K133(:,k)]);
    K14 = diag([K141(:,k),K142(:,k),K143(:,k)]);
    
    s_lam1 = diag([abs(s(1,k))^lambda1(1,k),abs(s(2,k))^lambda1(2,k),abs(s(3,k))^lambda1(3,k)]);
    s_lam2 = diag([abs(s(1,k))^lambda2(1,k),abs(s(2,k))^lambda2(2,k),abs(s(3,k))^lambda2(3,k)]);
    
    dp_obs = v_id+K11*s(:,k)+(K12+K13*s_lam1+K14*s_lam2)*sat(s(:,k));
    p_obs(:,k+1) = p_obs(:,k)+dp_obs*self.T;
    
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
    
    f_obs(:,k+1) = inv(R_fi)*(K11*s(:,k+1)+(K12+K13*s_lam1+K14*s_lam2)*sat(s(:,k+1)));
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Comparative Simulation 1: MPC only%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [A1, B1, M1, J_p1, R_fi1]=get_model(X_K1(3,k), self, X_K1(6,k));

    Q1 = diag([300,300,300,100,100,100]);

    F1 = 100*eye(n);

    R1 = 1*eye(p);
    


    ref=zeros(n,P+1);
    ref(:,1)=REF(:,k);

    ref(4,:)=ref(4,1);
    ref(5,:)=ref(5,1);
    ref(6,:)=ref(6,1);
    % 预测位置
    for i=1:P
        %s(k)=s(k-1)+v*t
        ref(1,i+1)=ref(1,i)+ref(4,i)*self.T;
        ref(2,i+1)=ref(2,i)+ref(5,i)*self.T;
        ref(3,i+1)=ref(3,i)+ref(6,i)*self.T;
    end
    ref=reshape(ref,n*(P+1),1);

    [ft1 , H1] = get_MPC_Matrices(A1,B1,Q1,R1,F1,P,k,X_K1(:, k),ref,U_K1);

    umax = eye(P*p, 1);
    umax(:,:)=1;

    dU_k = quadprog(H1,ft1,[],[],Aeq,beq,-umax,umax);
    U_K1(:,k+1) = U_K1(:,k)+dU_k(1:4,1);

    delta_w=inv(M)*(U_K1(:,k+1) - self.D*w1(1:4,k));
    w1(1:4,k+1)=w1(1:4,k)+delta_w*self.T;  

    X_K1(4:6,k+1) = self.r*R_fi1*J_p1*w1(1:4,k+1)+R_fi1*f(:,k);
    X_K1(1:3,k+1) = X_K1(1:3,k) +  X_K1(4:6,k+1)*self.T;
    
    wre1(:,k+1) = Jv*X_K1(4:6,k+1);
    
    
    %% Enter the next cycle
    k=k+1
end

%% drew picture
% figure("name","控制量状态量",'Position',[100,560,560,420]);
% subplot (3,1,1);
% hold;
% for i = 1:3
%     plot(X_K(i,:));
% end
% legend("x1","x2","x3")
% xlabel('时间t');
% ylabel('状态变量');
% hold off;
% 
% subplot (3,1,2);
% hold;
% for i = 1:3
%     plot(X_K(i+3,:));
% %     xlim([0,1100]);
% %     ylim([-0.5,0.5]);
% end
% legend("x4","x5","x6")
% xlabel('时间t');
% ylabel('状态变量');
% hold off;
% 
% subplot(3,1,3);
% hold;
% for i = 1:size(U_K,1)
%     plot(U_K(i,:));
% end
% legend("u1","u2","u3","u4")
% xlabel('时间t');
% ylabel('输入变量');
% hold off;


figure("name","trajectory");
t=1:size(X_K,2);
plot(REF(1,:),REF(2,:),'-','LineWidth',2);
%plot3(REF(1,:),REF(2,:),t,'-','LineWidth',2);
hold on;
plot(X_K1(1,:),X_K1(2,:),'--','LineWidth',2);
%plot3(X_K(1,:),X_K(2,:),t,'--','LineWidth',2);
hold on;
plot(X_K(1,:),X_K(2,:),':','LineWidth',2);
%plot3(X_K(1,:),X_K(2,:),t,':','LineWidth',2);
plot(X_K(1,1),X_K(2,1),'rd','LineWidth',2);
xlim([3.5,10]);
ylim([-1.5,4]);
legend("Reference","MPC","AVPSMO-MPC",'NumColumns',3,'FontSize',14,'Interpreter', 'latex');
xlabel('X (m)','FontSize',16,'Fontname','Times');
ylabel('Y (m)','FontSize',16,'Fontname','Times');
%创建一个图窗内的图窗（左间距比例，下间距比例，宽，高）
subcirc = axes('Position',[0.4 0.14 0.10 0.10]); 
axis square;  subcirc.XAxis.Visible = 'off';    subcirc.YAxis.Visible = 'off';  % 去掉小图窗的坐标轴
viscircles([0.4,pi/2], 0.2,'Color','k','LineWidth',1);  % 画一个圆
set(gca,'box','off','color','none');  % 使圆的底色透明
annotation('arrow',[0.46 0.5],[0.26 0.4]);  % 画一个箭头（箭头样式、x轴指向、y轴指向）
%创建一个图窗内的图窗（左间距比例，下间距比例，宽，高）
subfig = axes('Position',[0.45 0.47 0.3 0.12]);
axes(subfig);  % 不知道有啥用
set(gca,'FontSize',12, 'FontName','Times New Roman');  % 不知道有啥用
hold on; box on;
plot(REF(1,:),REF(2,:),'-','LineWidth',2);  % 在新建的小图窗内画图
plot(X_K1(1,:),X_K1(2,:),'--','LineWidth',2);
plot(X_K(1,:),X_K(2,:),':','LineWidth',2);
set(subfig,'xlim',[5.9,6.6],'ylim',[-1.25,-0.75],'color','none');  % 定义小图窗内的坐标范围
hold off;

figure("name","Position Error")
subplot (3,1,1);
plot(0:0.1:159.9,REF(1,1:N)-X_K1(1,1:N),'LineWidth',2);
hold on;
plot(0:0.1:159.9,REF(1,1:N)-X_K(1,1:N),'--','LineWidth',2);
xlim([0,160]);
ylim([-0.2,0.61]);
ylabel('$e_{\rm x}$(m)','Interpreter','latex','FontSize',16);
legend("MPC","AVPSMO-MPC",'FontSize',14,'NumColumns',2,'Interpreter', 'latex');
subplot (3,1,2);
plot(0:0.1:159.9,REF(2,1:N)-X_K1(2,1:N),'LineWidth',2);
hold on;
plot(0:0.1:159.9,REF(2,1:N)-X_K(2,1:N),'--','LineWidth',2);
xlim([0,160]);
ylim([-0.2,0.1]);
ylabel('$e_{\rm y}$(m)','Interpreter','latex','FontSize',16);
subplot (3,1,3);
plot(0:0.1:159.9,REF(3,1:N)-X_K1(3,1:N),'LineWidth',2);
hold on;
plot(0:0.1:159.9,REF(3,1:N)-X_K(3,1:N),'--','LineWidth',2);
xlim([0,160]);
ylim([-0.05,0.05]);
ylabel('$e_{\omega}{\rm (m)}$','Interpreter','latex','FontSize',16);
xlabel('Time(s)','FontSize',16,'Fontname','Times');

figure("name","Disturbance observation error")
for i = 1:3
    subplot (3,1,i);
    hold;
    plot(f(i,1:N)-f_obs(i,1:N));
    xlim([0,N]);
    ylim([-0.05,0.05]);
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

figure4=figure('Color',[1 1 1]);
set(gca,'FontSize',12);
hold on; box on;
plot(0:0.1:159.9,wre(1,1:N)+wre(2,1:N),'LineWidth',2);
hold on;
plot(0:0.1:159.9,wre(3,1:N)+wre(4,1:N),'--','LineWidth',2);
xlabel('Time (s)','FontSize',16,'Fontname','Times');
ylabel('Velocity (rad/s)','FontSize',16,'Fontname','Times');
h=legend({'${\bar w}_{1 \rm f}+{\bar w}_{2 \rm f}$','${\bar w}_{3 \rm f}+{\bar w}_{4 \rm f}$'},'FontSize',16)
set(h,'Interpreter','latex');

figure9=figure('Color',[1 1 1]);
set(gca,'FontSize',12);
hold on; box on;
plot(0:0.1:159.9,w_leader(1,1:N),'LineWidth',2);
hold on;
plot(0:0.1:159.9,wre(1,1:N),':','LineWidth',2);
hold on;
plot(0:0.1:159.9,wre1(1,1:N),'--','LineWidth',2);
hold on;
ylim([-5,15]);
xlabel('Time (s)','FontSize',16,'Fontname','Times');
ylabel('Velocity (rad/s)','FontSize',16,'Fontname','Times');
h=legend({'$w_{1\rm r}$','MPC ${\bar w}_{1}$','AVPSMO-MPC ${\bar w}_{1}$'},'FontSize',14);
set(h,'Interpreter','latex')

figure10=figure('Color',[1 1 1]);
set(gca,'FontSize',12);
hold on; box on;
plot(0:0.1:159.9,w_leader(2,1:N),'LineWidth',2);
hold on;
plot(0:0.1:159.9,wre(2,1:N),':','LineWidth',2);
hold on;
plot(0:0.1:159.9,wre1(2,1:N),'--','LineWidth',2);
hold on;
ylim([-5,15]);
xlabel('Time (s)','FontSize',16,'Fontname','Times');
ylabel('Velocity (rad/s)','FontSize',16,'Fontname','Times');
h=legend({'$w_{2\rm r}$','MPC ${\bar w}_{2}$','AVPSMO-MPC ${\bar w}_{2}$'},'FontSize',14);
set(h,'Interpreter','latex')

figure11=figure('Color',[1 1 1]);
set(gca,'FontSize',12);
hold on; box on;
plot(0:0.1:159.9,w_leader(3,1:N),'','LineWidth',2);
hold on;
plot(0:0.1:159.9,wre(3,1:N),':','LineWidth',2);
hold on;
plot(0:0.1:159.9,wre1(3,1:N),'--','LineWidth',2);
hold on;
ylim([-5,15]);
xlabel('Time (s)','FontSize',16,'Fontname','Times');
ylabel('Velocity (rad/s)','FontSize',16,'Fontname','Times');
h=legend({'$w_{3\rm r}$','MPC ${\bar w}_{3}$','AVPSMO-MPC ${\bar w}_{3}$'},'FontSize',14);
set(h,'Interpreter','latex')

figure12=figure('Color',[1 1 1]);
set(gca,'FontSize',12);
hold on; box on;
plot(0:0.1:159.9,w_leader(4,1:N),'LineWidth',2);
hold on;
plot(0:0.1:159.9,wre(4,1:N),':','LineWidth',2);
hold on;
plot(0:0.1:159.9,wre1(4,1:N),'--','LineWidth',2);
hold on;
ylim([-5,15]);
xlabel('Time (s)','FontSize',16,'Fontname','Times');
ylabel('Velocity (rad/s)','FontSize',16,'Fontname','Times');
h=legend({'$w_{4 \rm r}$','MPC ${\bar w}_{4}$','AVPSMO-MPC ${\bar w}_{4}$'});
set(h,'Interpreter','latex');