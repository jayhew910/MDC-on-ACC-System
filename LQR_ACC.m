%% ACC LQR MATLAB
clear all;
clc;
%% 系统特征
T_eng    = 0.460;
K_eng    = 0.732;
A_f      = -1/T_eng;
B_f      = -K_eng/T_eng;
T_hw     = 1.6;
Ts       = 0.05;
T_total  = 10; %仿真时间
T        = T_total/Ts;
v0       = 15;

vh = host_velocity(v0, T_total);

%% 状态空间&离散化系统

At    = [0 1 -T_hw; 0 0 -1; 0 0 A_f];
Bt    = [0; 0; B_f];
C_f   = eye(3);
D     = zeros(3,1);

sys1  = ss(At,Bt,C_f,D);
sys2  = c2d(sys1,Ts,'zoh');
A     = sys2.A;
B     = sys2.B;
C     = sys2.C;
D     = sys2.D;

%% LQR

x0 = [5;5;15];

% 参考序列
r = [ zeros(1,T);
      zeros(1,T);
      zeros(1,T);
      zeros(1,T)];
  
Q = transpose(C)*C;
R = 1;
N = 0;
[K,S,e] = dlqr(A,B,Q,R,N); 

%根据LQR控制律定义闭环系统
sysd_cl_unnormalized = ss(A-B*K, B*K, C, [], Ts);


x = zeros(length(A(:,1)),T);
u = zeros(length(B(1,:)),T);
y = zeros(length(C(:,1)),T);
t = zeros(1,T);

x(:,1) = x0(:,1);
sat = @(s) min(max(s, -3), 5);

for k = 1:1:T
    t(k) = (k-1)*Ts;
    % 计算
    u(:,k) = -K*x(:,k);  
    u(:,k) = sat(u(:,k));   
    % 应用
    x(:,k+1) = A*x(:,k) + B*u(:,k);
    y(:,k)   = C*x(:,k);
end

states_trajectory = y';

%% 结果
plot_mpc(u,x,t);
hold on;
ACC_MPC_Final;
for i = 1:4
    subplot(4,1,i);
    legend({'LQR Output','MPC Output'});
end

subplot(411);
axis([0 10 -1 6]);

subplot(412);
axis([0 10 -1 6]);

subplot(413);
axis([0 10 -3 6]);
