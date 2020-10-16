clear
close all
clc

%% Example 1 info
teta     = pi/3;
T        = 3 ;
num_vec  = 4; % dimension of state vector
num_meas = 2; % number of measurment

A =[1 0 T 0  % process Equation
    0 1 0 T
    0 0 1 0
    0 0 0 1];

B =[1 0 0 0  % meaurment equation
    0 1 0 0];

initial_x = [1;1;0;0]; % initial of state vector X
initial_P = diag([4,4,3,3]);% initial of cov matrix P
PEst      = initial_P;  % (MCC_KF)

%% Mix_Noise
mu_n1_x = [-3;-3;-3;-3]; % mean of noise1 state vector
mu_n2_x = [2;2;2;2]; % mean of noise2 state vector
mu_n1_z = [-2;-2]; % mean of noise1 measurment
mu_n2_z = [2;2]; % mean of noise2 measurment

R_n1 = diag([0.1,0.1]);% cov of Noise1 measurment
R_n2 = diag([0.1,0.1]); % cov of Noise2 measurment
Q_n1 = diag([0.1,0.1,0.1,0.1]); % cov of Noise1 process
Q_n2 = diag([0.1,0.1,0.1,0.1]); % cov of Noise2 process

F = @(x)[x(1)+T*x(3);x(2)+T*x(4);x(3);x(4)];  % nonlinear state equations for UKF
H = @(x)[x(1);x(2)];                               % measurement equation for UKF

%% initiallization
tf                 = 300;
iter               = floor(tf/T); % length of signal
num_of_exprement   = 100; % number of itterations

%% Shot_Noise
num_shot_noise     = 15;
start_of_shotnoise = 60;
index_rand_shot    = [randi([start_of_shotnoise/T iter],1,num_shot_noise-1) 21];

%% Choice the type of noise
while true
    flag_noise = input('Please choose the type of noise (1 = shot noise, 2 = Gaussian mixture noise): ');
    if flag_noise == 1 || flag_noise == 2
        break;
    end
end

MSE = zeros(num_of_exprement,num_vec,iter);

%% 100 Monte Carlo simulation
for Numexper = 1:num_of_exprement
    
    disp(['Simulation # ',num2str(Numexper.'),'/',num2str(num_of_exprement)]);
    
   %% Simulation Init Paramter
    % noise
    Q = Q_n1;  
    R = R_n1;
    MeasErrX = sqrt(Q)*randn(num_vec,iter);
    MeasErrZ = sqrt(R)*randn(num_meas,iter);

    % initial conditions for system and estimator(x0 , xhat0)
    xTrue = initial_x; % real states
    xEst  = initial_x; % (CF): C-Filter By cinar 2012

    % define some arrays to store the main signal and the esitimation signals
    xEst_MCC = zeros(num_vec,iter); %(C-Filter)
    X_TRUE   = zeros(num_vec,iter);
    
    %% type of noise
    if flag_noise == 1
        [MeasErrX, Q, MeasErrZ, R] = Shot_Noise(MeasErrX, MeasErrZ, index_rand_shot, num_shot_noise, num_vec, num_meas);
    else
        [MeasErrX, Q, MeasErrZ, R] = Mix_Noise(MeasErrX, MeasErrZ, mu_n1_x, mu_n2_x, mu_n1_z, mu_n2_z, R_n1, R_n2, Q_n1, Q_n2, num_vec, num_meas, iter);
    end
    
    %% one time simulation
    for t = 1 : 1 : iter
       %%  ======================= Observation =======================
        z = B * xTrue;
        z = z + MeasErrZ(:,t);

       %%  ======================= UKF ========================
        [xEst, PEst] = UKF(F, xEst, PEst, H, z, Q, R);
        xEst_MCC(:, t) = xEst;

       %% ======================= The System ========================
        xTrue = A * xTrue + MeasErrX(:, t);
        X_TRUE(:, t + 1) = xTrue;
    end
    
    X_TRUE(:, iter)   = [];
    MSE(Numexper,:,:) = (X_TRUE - xEst_MCC).^2;
end

%% RMSE
RMSE = zeros(num_vec,iter);

for i = 1 : iter
    RMSE(:,i) = sqrt(mean(MSE(:,:,i)))';
end

MMSE  = mean(RMSE,2);
MMMSE = mean(MMSE);
max   = max(RMSE.');

%% Plot data
DrowGraph(T, tf, flag_noise, num_shot_noise, index_rand_shot, mu_n1_x, mu_n2_x, mu_n1_z, mu_n2_z, R_n1, R_n2, Q_n1, Q_n2, MeasErrZ, RMSE);

%% Table result
disp('Mean square error : ');
disp('           x1          x2          x3        x4');
disp(['CF     : ',num2str(MMSE.'),'']);

%% UKF
function [x, P] = UKF(fstate, x, P, hmeas, z, Q, R)
% UKF   Unscented Kalman Filter for nonlinear dynamic systems
% [x, P] = ukf(f,x,P,h,z,Q,R) returns state estimate, x and state covariance, P 
% for nonlinear dynamic system (for simplicity, noises are assumed as additive):
%           x_k+1 = f(x_k) + w_k
%           z_k   = h(x_k) + v_k
% where w ~ N(0,Q) meaning w is gaussian noise with covariance Q
%       v ~ N(0,R) meaning v is gaussian noise with covariance R
% Inputs:   f: function handle for f(x)
%           x: "a priori" state estimate
%           P: "a priori" estimated state covariance
%           h: fanction handle for h(x)
%           z: current measurement
%           Q: process noise covariance 
%           R: measurement noise covariance
% Output:   x: "a posteriori" state estimate
%           P: "a posteriori" state covariance
%
% Example:
%{
n=3;      %number of state
q=0.1;    %std of process 
r=0.1;    %std of measurement
Q=q^2*eye(n); % covariance of process
R=r^2;        % covariance of measurement  
f=@(x)[x(2);x(3);0.05*x(1)*(x(2)+x(3))];  % nonlinear state equations
h=@(x)x(1);                               % measurement equation
s=[0;0;1];                                % initial state
x=s+q*randn(3,1); %initial state          % initial state with noise
P = eye(n);                               % initial state covraiance
N=20;                                     % total dynamic steps
xV = zeros(n,N);          %estmate        % allocate memory
sV = zeros(n,N);          %actual
zV = zeros(1,N);
for k=1:N
  z = h(s) + r*randn;                     % measurments
  sV(:,k)= s;                             % save actual state
  zV(k)  = z;                             % save measurment
  [x, P] = ukf(f,x,P,h,z,Q,R);            % ekf 
  xV(:,k) = x;                            % save estimate
  s = f(s) + q*randn(3,1);                % update process 
end
for k=1:3                                 % plot results
  subplot(3,1,k)
  plot(1:N, sV(k,:), '-', 1:N, xV(k,:), '--')
end
%}
% Reference: Julier, SJ. and Uhlmann, J.K., Unscented Filtering and
% Nonlinear Estimation, Proceedings of the IEEE, Vol. 92, No. 3,
% pp.401-422, 2004. 

L=numel(x);                                 %numer of states
m=numel(z);                                 %numer of measurements
alpha=1e-3;                                 %default, tunable
ki=0;                                       %default, tunable
beta=2;                                     %default, tunable
lambda = alpha^2*(L+ki)-L;                   %scaling factor
c = L+lambda;                                %scaling factor
Wm=[lambda/c 0.5/c+zeros(1,2*L)];           %weights for means
Wc=Wm;
Wc(1)=Wc(1)+(1-alpha^2+beta);               %weights for covariance
c= sqrt(c);

X=sigmas(x,P,c);                            %sigma points around x
[x1,X1,P1,X2]=ut(fstate,X,Wm,Wc,L,Q);          %unscented transformation of process
% X1=sigmas(x1,P1,c);                         %sigma points around x1
% X2=X1-x1(:,ones(1,size(X1,2)));             %deviation of X1
[z1,Z1,P2,Z2]=ut(hmeas,X1,Wm,Wc,m,R);       %#ok<ASGLU> %unscented transformation of measurments
P12=X2*diag(Wc)*Z2';                        %transformed cross-covariance
K=P12/(P2);
x=x1+K*(z-z1);                              %state update
P=P1-K*P2*K';                                %covariance update
end

function [y,Y,P,Y1]=ut(f,X,Wm,Wc,n,R)
%Unscented Transformation
%Input:
%        f: nonlinear map
%        X: sigma points
%       Wm: weights for mean
%       Wc: weights for covraiance
%        n: numer of outputs of f
%        R: additive covariance
%Output:
%        y: transformed mean
%        Y: transformed smapling points
%        P: transformed covariance
%       Y1: transformed deviations

L=size(X,2);
y=zeros(n,1);
Y=zeros(n,L);
for k=1:L                   
    Y(:,k)=f(X(:,k));       
    y=y+Wm(k)*Y(:,k);       
end
Y1=Y-y(:,ones(1,L));
P=Y1*diag(Wc)*Y1'+R;          
end

function X=sigmas(x,P,c)
%Sigma points around reference point
%Inputs:
%       x: reference point
%       P: covariance
%       c: coefficient
%Output:
%       X: Sigma points

A = c*chol(P)';      % A'A=P  (sqrt of P)
Y = x(:,ones(1,numel(x)));
X = [x Y+A Y-A]; 
end

%% Shot Noise
function [MeasErrX, Q, MeasErrZ, R] = Shot_Noise(MeasErrX, MeasErrZ, index_rand_shot, num_shot_noise, num_vec, num_meas)
    %%%% Generate Non-guassian Noise(Shot Noise)
    % In process equation
    for i=1:num_shot_noise
        temp = randi([1 5],num_vec,1);
        MeasErrX(:,index_rand_shot(i)) =  MeasErrX(:,index_rand_shot(i)) + temp + 0.6*randn(num_vec,1);
    end
    Q = cov(MeasErrX');

    % In Measurment equation
    for i=1:num_shot_noise
        temp = randi([1 5],num_meas,1);
        MeasErrZ(:,index_rand_shot(i)) =  MeasErrZ(:,index_rand_shot(i)) + temp + 0.6*randn(num_meas,1); 
    end
    R = cov(MeasErrZ');
end

%% Mix_Noise
function [MeasErrX, Q, MeasErrZ, R] = Mix_Noise(MeasErrX, MeasErrZ, mu_n1_x, mu_n2_x, mu_n1_z, mu_n2_z, R_n1, R_n2, Q_n1, Q_n2, num_vec, num_meas, iter)
    %%%%% Generate Non-guassian Noise (guassian mixture)
    % In  process equation
    for i = 21 : iter
        Mixt = [mu_n1_x+ sqrt(Q_n1)*randn(num_vec,1),mu_n1_x + sqrt(Q_n1)*randn(num_vec,1),mu_n1_x+ sqrt(Q_n1)*randn(num_vec,1),mu_n1_x+ sqrt(Q_n1)*randn(num_vec,1),mu_n1_x+ sqrt(Q_n1)*randn(num_vec,1),mu_n2_x+ sqrt(Q_n2)*randn(num_vec,1),mu_n2_x+ sqrt(Q_n2)*randn(num_vec,1),mu_n2_x+ sqrt(Q_n2)*randn(num_vec,1),mu_n2_x+ sqrt(Q_n2)*randn(num_vec,1),mu_n2_x+ sqrt(Q_n2)*randn(num_vec,1)];
        rand_index = randi([1 10],1,1);
        MeasErrX(:,i) = Mixt(:,rand_index); 
    end
    Q = cov(MeasErrX');

    % In Measurment equation
    for i = 21 : iter
        Mixt = [ mu_n1_z+ sqrt(R_n1)*randn(num_meas,1),mu_n1_z+ sqrt(R_n1)*randn(num_meas,1),mu_n1_z+ sqrt(R_n1)*randn(num_meas,1),mu_n1_z+ sqrt(R_n1)*randn(num_meas,1),mu_n1_z+ sqrt(R_n1)*randn(num_meas,1),mu_n2_z+ sqrt(R_n2)*randn(num_meas,1),mu_n2_z+ sqrt(R_n2)*randn(num_meas,1),mu_n2_z+ sqrt(R_n2)*randn(num_meas,1),mu_n2_z+ sqrt(R_n2)*randn(num_meas,1),mu_n2_z+ sqrt(R_n2)*randn(num_meas,1)];
        rand_index = randi([1 10],1,1);
        MeasErrZ(:,i) = Mixt(:,rand_index); 
    end
    R = cov(MeasErrZ');
end

%% Plot data
function [] = DrowGraph(T, tf, flag_noise, num_shot_noise, index_rand_shot, mu_n1_x, mu_n2_x, mu_n1_z, mu_n2_z, R_n1, R_n2, Q_n1, Q_n2, MeasErrZ, MSE_CF)
    SetPlotOptions;
    t = 0 : T : tf-T;
    c = 3;

    if flag_noise == 1 % ploting in Shot noise case
        figure(1)
        xlabel('time');
        ylabel('Value of shot noise');
        hold on
        X = index_rand_shot;
        Y = MeasErrZ(1,index_rand_shot(1:end));
        for i=1:num_shot_noise
            bar(X(i)*T,Y(i),1,'b','EdgeColor','b');
        end

        figure(2)
        hold on
        plot(t,MSE_CF(c,:) ,'r');
        legend('C-Filter');
        xlabel('time'); 
        ylabel('Root mean square error')

    else % ploting in mixture guassian noise case

        f1 =  @(s) 1/sqrt(2*pi*sqrt(Q_n1(1,1)))*(exp((-(s-mu_n1_x(1,1))*(s-mu_n1_x(1,1)))/(2*sqrt(Q_n1(1,1))))+exp((-(s-mu_n2_x(1,1))*(s-mu_n2_x(1,1)))/(2*sqrt(Q_n2(1,1)))));
        f2 =  @(s) 1/sqrt(2*pi*sqrt(R_n1(1,1)))*(exp((-(s-mu_n1_z(1,1))*(s-mu_n1_z(1,1)))/(2*sqrt(R_n1(1,1))))+exp((-(s-mu_n2_z(1,1))*(s-mu_n2_z(1,1)))/(2*sqrt(R_n2(1,1)))));

        figure(1)
        subplot(2,1,1);
        fplot(f1,[-6,6]);title('pdf of process noise');
        subplot(2,1,2);
        fplot(f2,[-6,6]);title('pdf of measurment noise');

        figure(2)
        hold on
        plot(t,MSE_CF(c,:) ,'r'); 
        legend('C-Filter');
        xlabel('time');
        ylabel('Root mean square error');
    end
end
