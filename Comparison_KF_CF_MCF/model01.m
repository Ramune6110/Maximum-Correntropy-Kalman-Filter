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
num_shot_noise     = 15;
start_of_shotnoise = 60;
index_rand_shot    = [randi([start_of_shotnoise/T iter],1,num_shot_noise-1) 21];