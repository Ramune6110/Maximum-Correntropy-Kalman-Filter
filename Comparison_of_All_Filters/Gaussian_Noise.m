%%%% Generate guassian Noise
%% In process equation
Q = Q_n1;  
MeasErrX = sqrt(Q)*randn(num_vec,iter);
        
%% In Measurment equation
R = R_n1;
MeasErrZ = sqrt(R)*randn(num_meas,iter);