%% Example 3 info
teta = pi/3;
T = 3;
num_vec = 2; % dimension of state vector
num_meas = 1; % number of measurment
A=[cos(teta) sin(teta);   % process Equation   [0.5000    0.8660
    -sin(teta) cos(teta)];    %                 -0.8660    0.5000]
B =[1 1];  % measurment equation
initial_x = [1;1]; % initial of state vector X
initial_P = diag([4,4]);% initial of cov matrix P
mu_n1_x = [0;0]; % mean of noise1 state vector
mu_n2_x = [0;0]; % mean of noise2 state vector
mu_n1_z = [0]; % mean of noise1 measurment
mu_n2_z = [0]; % mean of noise2 measurment
Q_n1 = [0.01 0 ;0 0.01]; % cov of Noise1 measurment
Q_n2 = [0.01 0 ;0 0.01 ]; % cov of Noise2 measurment
R_n1 = [0.01 ]; % cov of Noise1 process
R_n2 = [0.01 ]; % cov of Noise2 process
F = @(x)[cos(teta)*x(1)+sin(teta)*x(2);-sin(teta)*x(1)+cos(teta)*x(2)];  % nonlinear state equations for UKF
H = @(x)[x(1)+x(2)];                               % measurement equation for UKF
%% initiallization
tf                 = 300;
iter               = floor(tf/T); % length of signal
num_of_exprement   = 100; % number of itterations
num_shot_noise     = 15;
start_of_shotnoise = 60;
index_rand_shot    = [randi([start_of_shotnoise/T iter],1,num_shot_noise-1) 21];