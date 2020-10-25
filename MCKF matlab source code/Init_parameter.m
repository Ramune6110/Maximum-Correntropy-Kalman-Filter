%% Init
xTrue = [0;0];    %initial true state
xEst_KF    = [1;1];    %initial estimated state
xEst_MCKF  = [1;1];    %initial estimated state
xEst_LRKF  = [1;1];    %initial estimated state
PEst_KF    = eye(2);   %initial state covariance
PEst_MCKF  = eye(2);   %initial state covariance
PEst_LRKF  = eye(2);   %initial state covariance

Q  = diag([0.01,0.01]);     %covariance matrix of process noise
R  = 0.01;         %covariance matrix of measurement noi

W = 1;

MeasErrX = sqrt(Q)*randn(n,N);
MeasErrZ = sqrt(R)*randn(m,N);