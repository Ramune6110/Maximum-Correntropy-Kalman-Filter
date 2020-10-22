% Simulation Init Paramter
% noise
Q = Q_n1;  
R = R_n1;
MeasErrX = sqrt(Q)*randn(num_vec,iter);
MeasErrZ = sqrt(R)*randn(num_meas,iter);

% initial conditions for system and estimator(x0 , xhat0)
x     = initial_x; % real states
xhat2 = initial_x; % (CF): C-Filter By cinar 2012
xhat3 = initial_x; % (MCC_KF) : second Algorithm (Weighted maximum correntropy)
xhat4 = initial_x; % (UKF)
xhat6 = initial_x; % (EnKF)
xhat7 = initial_x; % (GSF)
xhat8 = initial_x; % (MCF) : first Algorithm (Modified C-filter)

% elements of initial estimation error covariance (P0hat)
P_MCC_KF = initial_P;  % (MCC_KF)
P_UKF    = initial_P;  % (UKF)
P_EnKF   = initial_P;  % (EnKF)
P_GSF    = initial_P;  % (GSF)

%%initialization of the estimators
% number of samples in Enkf
N_Enkf = 100;
no = sqrt(P_EnKF) * repmat(randn(size(xhat6)),1,N_Enkf);
X_enkf = repmat(xhat6,1,N_Enkf) + no;