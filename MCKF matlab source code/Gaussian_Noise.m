%%%% Generate guassian Noise
mu_n1_x = [0;0]; % mean of noise1 state vector
mu_n2_x = [0;0]; % mean of noise2 state vector
mu_n1_z = 0; % mean of noise1 measurment
mu_n2_z = 0; % mean of noise2 measurment

Q_n1 = [0.1 0 ;0 0.1]; % cov of Noise1 measurment
Q_n2 = [0.1 0 ;0 0.1 ]; % cov of Noise2 measurment
R_n1 = 0.01; % cov of Noise1 process
R_n2 = 0.01; % cov of Noise2 process
%% In process equation
MeasErrX = sqrt(Q)*randn(n,N);
        
%% In Measurment equation
MeasErrZ = sqrt(R)*randn(m,N);