clc;
clear all;

T_eng     = 0.460;
K_eng     = 0.732;
A_f       = -1/T_eng;
B_f       = -K_eng/T_eng;
C_f       = eye(3);
T_hw      = 1.6;
Ts        = 0.05;
T_total   = 10;
T         = T_total/Ts;
v0        = 15;           % 初始目标速度
init_dist = 5;            % 初始距离


%% 离散化系统
At    = [0 1 -T_hw; 0 0 -1; 0 0 A_f];
Bt    = [0; 0; B_f];
sys1  = ss(At,Bt,C_f,0);
sys2  = c2d(sys1,Ts,'zoh');
A     = sys2.A;
B     = sys2.B;
C     = sys2.C;
% 自适应巡航车辆参数
vh = host_velocity(v0,T);               % 速度
sec_row = [init_dist*ones(1,length(vh))];
x0      = [sec_row;                     % 距离误差
           sec_row;                     % 速度误差
           vh];                         % 初始速度           
u  = zeros(1,length(vh));               % 输入加速度

xr = [zeros(1,length(vh));       
      zeros(1,length(vh));       
      zeros(1,length(vh));];

x     = [zeros(1,T_total);zeros(1,T_total);zeros(1,T_total)];
s_lb  = [0;0;15];
s_ub  = [2;2.5;40]; 

LTI.A = A;
LTI.B = B;
LTI.C = eye(3);

%系统维度的定义
dim.nx = length(A);     % 状态维度
dim.ny = length(B);     % 输出维度
dim.nu = 1;             % 输入维度
dim.N  = 20;            
N      = dim.N; 

Q = C'*C;   % 输出权重
R = 1;      % 输入权重

%% 预测模型和成本函数
[P,S] = predmodgen(LTI,dim);
[H,h] = costgen(P,S,Q,R,dim);
%% 模型预测控制

umin    = -3*ones(1,N);
umax    = 5*ones(1,N);
xr(:,1) = x0(:,1);
y(:,1)  = C*x0(:,1);
for i = 1:T
    t(i) = (i-1)*Ts;
    f = h*xr(:,i);
    Aueq = C*x0(:,1:N);
    bueq = y(:,1);
    options = optimoptions('quadprog','Display','off');
    warning off;
    [ures,~,exitflag] = quadprog(H,f,[],[],Aueq,bueq,umin,umax,[],options);
    u(i) = ures(1);
    xr(:,i+1)  = A*xr(:,i) + B*u(i);
        y(:,i) = C*xr(:,i+1);
        x(:,1) = xr(:,i);
        for j = 1:N-1
            x(:,j+1)   = A*x(:,j)+B*ures(j);
            y_res(:,j) = C*x(:,j);
        end
end

%% 结果
plot_mpc(u,xr,t);

for i = 1:4
    subplot(4,1,i);
    legend({'init\_dist=5m'});
end

subplot(412);
axis([0 10 -2 6]);

subplot(413);
axis([0 10 -4 6]);

subplot(414);
axis([0 10 -2 15]);