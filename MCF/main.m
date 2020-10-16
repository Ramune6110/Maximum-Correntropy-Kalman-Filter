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
    xEst_MCF = zeros(num_vec,iter); %(C-Filter)
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

       %%  ======================= MC_Filter ========================
        xEst = MCF(xEst, z, A, B, num_vec);
        xEst_MCF(:, t) = xEst;

       %% ======================= The System ========================
        xTrue = A * xTrue + MeasErrX(:, t);
        X_TRUE(:, t + 1) = xTrue;
    end
    
    X_TRUE(:, iter)   = [];
    MSE(Numexper,:,:) = (X_TRUE - xEst_MCF).^2;
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

%% MCF
function xEst = MCF(xEst, z, A, B, num_vec)
    xEst       = A * xEst;
    innov      = z - B * xEst;
    norm_innov = sqrt((innov)'*(innov));
    sigma      = 1 * norm_innov;
    K          = exp(-(norm_innov)^2/(2 * sigma^2));
    Gain       = pinv(eye(num_vec) + K * B'*B ) * K * B';
    xEst       = xEst + Gain *(innov);
    innov      = z - B * xEst;
    norm_innov = sqrt((innov)'*(innov));
    sigma      = 1 * norm_innov;
    K          = exp(-(norm_innov)^2/(2 * sigma^2));
    Gain       = pinv(eye(num_vec) + K * B' * B ) * K * B';
    xEst       = xEst + Gain * (innov);
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
