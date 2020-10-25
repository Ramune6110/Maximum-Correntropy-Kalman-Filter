% %% Generate Non-guassian Noise (guassian mixture)
mu_n1_x = [0;0]; % mean of noise1 state vector
mu_n2_x = [0;0]; % mean of noise2 state vector
mu_n1_z = 0; % mean of noise1 measurment
mu_n2_z = 0; % mean of noise2 measurment

Q_n1 = [0.1 0 ;0 0.1]; % cov of Noise1 measurment
Q_n2 = [0.1 0 ;0 0.1 ]; % cov of Noise2 measurment
R_n1 = 0.01; % cov of Noise1 process
R_n2 = 0.5; % cov of Noise2 process

% %In  process equation
% for i = 1 : N
%     Mixt = [mu_n1_x+ sqrt(Q_n1)*randn(n,1),mu_n1_x + sqrt(Q_n1)*randn(n,1),mu_n1_x+ sqrt(Q_n1)*randn(n,1),mu_n1_x+ sqrt(Q_n1)*randn(n,1),mu_n1_x+ sqrt(Q_n1)*randn(n,1),mu_n2_x+ sqrt(Q_n2)*randn(n,1),mu_n2_x+ sqrt(Q_n2)*randn(n,1),mu_n2_x+ sqrt(Q_n2)*randn(n,1),mu_n2_x+ sqrt(Q_n2)*randn(n,1),mu_n2_x+ sqrt(Q_n2)*randn(n,1)];
%     rand_index = randi([1 10],1,1);
%     MeasErrX(:,i) = Mixt(:,rand_index); 
% end
% Q = cov(MeasErrX');

% In Measurment equation
for i = 1 : N
    Mixt = [ mu_n1_z+ sqrt(R_n1)*randn(m,1),mu_n1_z+ sqrt(R_n1)*randn(m,1),mu_n1_z+ sqrt(R_n1)*randn(m,1),mu_n1_z+ sqrt(R_n1)*randn(m,1),mu_n1_z+ sqrt(R_n1)*randn(m,1),mu_n2_z+ sqrt(R_n2)*randn(m,1),mu_n2_z+ sqrt(R_n2)*randn(m,1),mu_n2_z+ sqrt(R_n2)*randn(m,1),mu_n2_z+ sqrt(R_n2)*randn(m,1),mu_n2_z+ sqrt(R_n2)*randn(m,1)];
    rand_index = randi([1 10],1,1);
    MeasErrZ(:,i) = Mixt(:,rand_index);
end
R = cov(MeasErrZ');