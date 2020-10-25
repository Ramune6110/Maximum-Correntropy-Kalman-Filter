clear;
close all;
clc

%% Init
n  = 2;        %number of states
m  = 1;        %number of measurements
N  = 100;     %total number of time steps
F  = [cos(pi/18),-sin(pi/18);sin(pi/18),cos(pi/18)];  %state transition matrix
H  = [1,1];    %observation matrix
xTrue = [0;0];    %initial true state
xEst  = [1;1];    %initial estimated state
PEst  = eye(2);   %initial state covariance

Q  = diag([0.01,0.01]);     %covariance matrix of process noise
R  = 0.01;         %covariance matrix of measurement noi

sV = zeros(n,N);  %true state value
xV = zeros(n,N);  %estimated state value
zV = zeros(m,N);  %true measurement value
bV = zeros(1,N);  %iteration number

mu_n1_x = [-2;-2]; % mean of noise1 state vector
mu_n2_x = [2;2]; % mean of noise2 state vector
mu_n1_z = -2; % mean of noise1 measurment
mu_n2_z = 2; % mean of noise2 measurment

Q_n1 = [0.1 0 ;0 0.1]; % cov of Noise1 measurment
Q_n2 = [0.1 0 ;0 0.1 ]; % cov of Noise2 measurment
R_n1 = 0.1; % cov of Noise1 process
R_n2 = 0.1; % cov of Noise2 process

MeasErrX = sqrt(Q)*randn(n,N);
MeasErrZ = sqrt(R)*randn(m,N);

%% Generate Non-guassian Noise (guassian mixture)
% In  process equation
for i = 1 : N
    Mixt = [mu_n1_x+ sqrt(Q_n1)*randn(n,1),mu_n1_x + sqrt(Q_n1)*randn(n,1),mu_n1_x+ sqrt(Q_n1)*randn(n,1),mu_n1_x+ sqrt(Q_n1)*randn(n,1),mu_n1_x+ sqrt(Q_n1)*randn(n,1),mu_n2_x+ sqrt(Q_n2)*randn(n,1),mu_n2_x+ sqrt(Q_n2)*randn(n,1),mu_n2_x+ sqrt(Q_n2)*randn(n,1),mu_n2_x+ sqrt(Q_n2)*randn(n,1),mu_n2_x+ sqrt(Q_n2)*randn(n,1)];
    rand_index = randi([1 10],1,1);
    MeasErrX(:,i) = Mixt(:,rand_index); 
end
Q = cov(MeasErrX');

% In Measurment equation
for i = 1 : N
    Mixt = [ mu_n1_z+ sqrt(R_n1)*randn(m,1),mu_n1_z+ sqrt(R_n1)*randn(m,1),mu_n1_z+ sqrt(R_n1)*randn(m,1),mu_n1_z+ sqrt(R_n1)*randn(m,1),mu_n1_z+ sqrt(R_n1)*randn(m,1),mu_n2_z+ sqrt(R_n2)*randn(m,1),mu_n2_z+ sqrt(R_n2)*randn(m,1),mu_n2_z+ sqrt(R_n2)*randn(m,1),mu_n2_z+ sqrt(R_n2)*randn(m,1),mu_n2_z+ sqrt(R_n2)*randn(m,1)];
    rand_index = randi([1 10],1,1);
    MeasErrZ(:,i) = Mixt(:,rand_index);
end
R = cov(MeasErrZ');

%% simulate Kalman FIlter
for k=1:N
    
    % Observation
    xTrue = F * xTrue + MeasErrX(:, k); %true state update process
    z = H * xTrue + MeasErrZ(:, k);     %true measurement update process
    
    % MCKF
    [xEst,PEst,b] = MCKF(F,xEst,PEst,H,z,Q,R); % KF

    % save data
    sV(:,k) = xTrue;           %save true state value
    zV(:,k) = z;               %save true measurement value
    xV(:,k) = xEst;            %save estimated state value
    bV(:,k) = b;               %save iteration number
end

%% plot results
figure(1);
plot(1:N,sV(1,:),'b-',1:N,xV(1,:),'r--');
legend('true','MCKF');
xlabel('{\ittime step k}','Fontsize',15,'Fontname','Times new roman');
title('state estimate x(1)','Fontsize',15,'Fontname','Times new roman');

figure(2);
plot(1:N,sV(2,:),'b-',1:N,xV(2,:),'r--')
legend('true','MCKF');
xlabel('{\ittime step k}','Fontsize',15,'Fontname','Times new roman');
title('state estimate x(2)','Fontsize',15,'Fontname','Times new roman');

f1 =  @(s) 1/sqrt(2*pi*sqrt(Q_n1(1,1)))*(exp((-(s-mu_n1_x(1,1))*(s-mu_n1_x(1,1)))/(2*sqrt(Q_n1(1,1))))+exp((-(s-mu_n2_x(1,1))*(s-mu_n2_x(1,1)))/(2*sqrt(Q_n2(1,1)))));
f2 =  @(s) 1/sqrt(2*pi*sqrt(R_n1(1,1)))*(exp((-(s-mu_n1_z(1,1))*(s-mu_n1_z(1,1)))/(2*sqrt(R_n1(1,1))))+exp((-(s-mu_n2_z(1,1))*(s-mu_n2_z(1,1)))/(2*sqrt(R_n2(1,1)))));
    
figure(3);
subplot(2,1,1);
fplot(f1,[-6,6]);title('pdf of process noise');
subplot(2,1,2);
fplot(f2,[-6,6]);title('pdf of measurment noise');

figure(4);
plot(1:N,bV(1,:),'b--');