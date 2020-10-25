%%%%% Generate Non-guassian Noise (guassian mixture)
%% In  process equation
for i = 1 : iter
    Mixt = [mu_n1_x+ sqrt(Q_n1)*randn(num_vec,1),mu_n1_x + sqrt(Q_n1)*randn(num_vec,1),mu_n1_x+ sqrt(Q_n1)*randn(num_vec,1),mu_n1_x+ sqrt(Q_n1)*randn(num_vec,1),mu_n1_x+ sqrt(Q_n1)*randn(num_vec,1),mu_n2_x+ sqrt(Q_n2)*randn(num_vec,1),mu_n2_x+ sqrt(Q_n2)*randn(num_vec,1),mu_n2_x+ sqrt(Q_n2)*randn(num_vec,1),mu_n2_x+ sqrt(Q_n2)*randn(num_vec,1),mu_n2_x+ sqrt(Q_n2)*randn(num_vec,1)];
    rand_index = randi([1 10],1,1);
    MeasErrX(:,i) = Mixt(:,rand_index); %#ok<SAGROW>
end
Q = cov(MeasErrX');

%% In Measurment equation
for i = 1 : iter
    Mixt = [ mu_n1_z+ sqrt(R_n1)*randn(num_meas,1),mu_n1_z+ sqrt(R_n1)*randn(num_meas,1),mu_n1_z+ sqrt(R_n1)*randn(num_meas,1),mu_n1_z+ sqrt(R_n1)*randn(num_meas,1),mu_n1_z+ sqrt(R_n1)*randn(num_meas,1),mu_n2_z+ sqrt(R_n2)*randn(num_meas,1),mu_n2_z+ sqrt(R_n2)*randn(num_meas,1),mu_n2_z+ sqrt(R_n2)*randn(num_meas,1),mu_n2_z+ sqrt(R_n2)*randn(num_meas,1),mu_n2_z+ sqrt(R_n2)*randn(num_meas,1)];
    rand_index = randi([1 10],1,1);
    MeasErrZ(:,i) = Mixt(:,rand_index); %#ok<SAGROW>
end
R = cov(MeasErrZ');
