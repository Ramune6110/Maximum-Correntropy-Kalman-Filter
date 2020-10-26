clear;
close all;
clc

%% Init
n = 2;        %number of states
m = 1;        %number of measurements
N = 100;      %total number of time steps
F = [cos(pi/18),-sin(pi/18);sin(pi/18),cos(pi/18)];  %state transition matrix
H = [1,1];    %observation matrix
xTrue     = [0;0];    %initial true state
xEst_KF   = [1;1];    %initial estimated state
xEst_MCKF = [1;1];    %initial estimated state
xEst_LRKF = [1;1];    %initial estimated state
PEst_KF   = eye(2);   %initial state covariance
PEst_MCKF = eye(2);   %initial state covariance
PEst_LRKF = eye(2);   %initial state covariance

Q = diag([0.01,0.01]);     %covariance matrix of process noise
R = 0.01;         %covariance matrix of measurement noi

% LRKF parameter
W = 1;

%% save array
XTRUE     = zeros(n,N);  %true state value
XEst_KF   = zeros(n,N);  %estimated state value
XEst_MCKF = zeros(n,N);  %estimated state value
XEst_LRKF = zeros(n,N);  %estimated state value
zV = zeros(m,N);  %true measurement value
bV = zeros(1,N);  %iteration number
LRKF_ita = zeros(1,N);  %iteration number
LRKF_W   = zeros(1, N);

LRKF_Gain  = zeros(n, N);
MCCKF_Gain = zeros(n, N);

% square error
SE_KF   = zeros(n,N);
SE_CF   = zeros(n,N);
SE_LRKF = zeros(n,N);

% noise array
MeasErrX = sqrt(Q)*randn(n,N);
MeasErrZ = sqrt(R)*randn(m,N);

%% type of noise
while true
    flag_noise = input('Please choose the type of noise (1 = Laplace noise, 2 = Gaussian mixture noise, 3 = Gaussian noise, 4 = Shot noise): ');
    if flag_noise == 1 || flag_noise == 2 || flag_noise == 3 || flag_noise == 4
        break;
    end
end

if flag_noise == 1
    Laplace_Noise;
elseif flag_noise == 2
    Mix_Noise; 
elseif flag_noise == 3
    Gaussian_Noise;
elseif flag_noise == 4
    Shot_Noise;
end

%% simulate Kalman FIlter
for k = 1 : N
    
    % Observation
    xTrue = F * xTrue + MeasErrX(:, k); %true state update process
    z = H * xTrue + MeasErrZ(:, k);     %true measurement update process
    
    % KF
    [xEst_KF, PEst_KF] = KF(xEst_KF, PEst_KF, z, F, H, Q, R, n);
    % MCKF
    %[xEst_MCKF, PEst_MCKF, b] = MCKF(F, xEst_MCKF ,PEst_MCKF ,H, z, Q, R); 
    [xEst_MCKF, PEst_MCKF, Gain_MCCKF, b] = MCCCF(xEst_MCKF, PEst_MCKF, z, F, H, Q, R, n);
    % LRKF
    [xEst_LRKF, PEst_LRKF, Gain_LRKF, in, W] = LRKF(xEst_LRKF, PEst_LRKF, z, F, H, Q, R, W);
    
    % save data
    XTRUE(:,k) = xTrue;        %save true state value
    zV(:,k) = z;               %save true measurement value
    XEst_KF(:,k)   = xEst_KF;            %save estimated state value
    XEst_MCKF(:,k) = xEst_MCKF;            %save estimated state value
    XEst_LRKF(:,k) = xEst_LRKF;            %save estimated state value
    bV(:,k) = b;               %save iteration number
    LRKF_ita(:,k) = in;        %save iteration number
    LRKF_W(:, k)  = W;
    LRKF_Gain(:, k)  = Gain_LRKF;
    MCCKF_Gain(:, k) = Gain_MCCKF;
end

SE_KF(:,:)   = (XTRUE - XEst_KF).^2;
SE_CF(:,:)   = (XTRUE - XEst_MCKF).^2;
SE_LRKF(:,:) = (XTRUE - XEst_LRKF).^2;

MSE_KF   = zeros(n,N);
MSE_CF   = zeros(n,N);
MSE_LRKF = zeros(n,N);

for i = 1 : N
    MSE_KF(:,i)   = sqrt(mean(SE_KF(:,i)))';
    MSE_CF(:,i)   = sqrt(mean(SE_CF(:,i)))';
    MSE_LRKF(:,i) = sqrt(mean(SE_LRKF(:,i)))';
end

MMSE_KF   = mean(MSE_KF,2);
MMSE_CF   = mean(MSE_CF,2);
MMSE_LRKF = mean(MSE_LRKF,2);

%% plot results
figure(1);
plot(1:N, XTRUE(1,:),'b-', 1:N,XEst_MCKF(1,:),'r--', 1:N,XEst_LRKF(1,:),'m--');
legend('true','MCCKF', 'LRKF');
xlabel('{\ittime step k}','Fontsize',15,'Fontname','Times new roman');
title('state estimate x(1)','Fontsize',15,'Fontname','Times new roman');

figure(2);
plot(1:N, XTRUE(2,:),'b-', 1:N,XEst_MCKF(2,:),'r--', 1:N,XEst_LRKF(2,:),'m--')
legend('true','MCCKF', 'LRKF');
xlabel('{\ittime step k}','Fontsize',15,'Fontname','Times new roman');
title('state estimate x(2)','Fontsize',15,'Fontname','Times new roman');

if flag_noise == 4
    figure(3);
    xlabel('time');
    ylabel('Value of shot noise');
    hold on
    X = index_rand_shot;
    Y = MeasErrZ(1,index_rand_shot(1:end));
    for i=1:num_shot_noise
        bar(X(i)*T,Y(i),1,'b','EdgeColor','b');
    end
else 
    f1 =  @(s) 1/sqrt(2*pi*sqrt(Q_n1(1,1)))*(exp((-(s-mu_n1_x(1,1))*(s-mu_n1_x(1,1)))/(2*sqrt(Q_n1(1,1))))+exp((-(s-mu_n2_x(1,1))*(s-mu_n2_x(1,1)))/(2*sqrt(Q_n2(1,1)))));
    f2 =  @(s) 1/sqrt(2*pi*sqrt(R_n1(1,1)))*(exp((-(s-mu_n1_z(1,1))*(s-mu_n1_z(1,1)))/(2*sqrt(R_n1(1,1))))+exp((-(s-mu_n2_z(1,1))*(s-mu_n2_z(1,1)))/(2*sqrt(R_n2(1,1)))));

    figure(3);
    subplot(2,1,1);
    fplot(f1,[-6,6]);title('pdf of process noise');
    subplot(2,1,2);
    fplot(f2,[-6,6]);title('pdf of measurment noise');
end

figure(4);
% plot(1:N,bV(1,:),'b--');
% plot(1:N,LRKF_ita(1,:),'m--');
plot(1:N,LRKF_W(1,:),'b--');

figure(5);
plot(1:N,SE_CF(1,:),'r--',1:N,SE_LRKF(1,:),'m--');
legend('KF','MCCKF', 'LRKF');

figure(6);
subplot(2, 1, 1);
plot(1:N, LRKF_Gain(1,:),'m--', 1:N, MCCKF_Gain(1,:),'r--');
ylim([0 0.8])
legend('LRKF', 'MCCKF');
xlabel('{\ittime step k}','Fontsize',15,'Fontname','Times new roman');
title('state estimate x(1)','Fontsize',15,'Fontname','Times new roman');

subplot(2, 1, 2);
plot(1:N, LRKF_Gain(2,:),'m--', 1:N, MCCKF_Gain(2,:),'r--');
ylim([0 0.8])
legend('LRKF', 'MCCKF');
xlabel('{\ittime step k}','Fontsize',15,'Fontname','Times new roman');
title('state estimate x(2)','Fontsize',15,'Fontname','Times new roman');

%% Table result

disp('Mean square error : ');
disp('           x1          x2');
disp(['KF     : ',num2str(MMSE_KF.'),'']);
disp(['CF     : ',num2str(MMSE_CF.'),'']);
disp(['LRKF   : ',num2str(MMSE_LRKF.'),'']);
