clear
close all
clc

%% call system dynamic
sysinfo;

%% initiallization
tf = 300;
iter = floor(tf/T); % length of signal
num_of_exprement = 100 ; % number of itterations
num_shot_noise = 15;
start_of_shotnoise = 60;
index_rand_shot = [randi([start_of_shotnoise/T iter],1,num_shot_noise-1) 21];

%% Choice the type of noise
while true
    flag_noise = input('Please choose the type of noise (1 = shot noise, 2 = Gaussian mixture noise): ');
    if flag_noise == 1 || flag_noise == 2
        break;
    end
end

SE_CF     = zeros(num_of_exprement,num_vec,iter);
SE_MCC_KF = zeros(num_of_exprement,num_vec,iter);
SE_UKF    = zeros(num_of_exprement,num_vec,iter);
SE_EnKF   = zeros(num_of_exprement,num_vec,iter);
SE_GSF    = zeros(num_of_exprement,num_vec,iter);
SE_MCF    = zeros(num_of_exprement,num_vec,iter);

% 100 Monte Carlo simulation
for Numexper = 1:num_of_exprement
    
    disp(['Simulation # ',num2str(Numexper.'),'/',num2str(num_of_exprement)]);
   
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
    
    %% define some arrays to store the main signal and the esitimation signals
    xhat_CF     = zeros(num_vec,iter); %(C-Filter)
    xhat_MCC_KF = zeros(num_vec,iter); % (MCC_KF)
    xhat_UKF    = zeros(num_vec,iter); % (UKF)
    xhat_EnKF   = zeros(num_vec,iter); % (EnKF)
    xhat_GSF    = zeros(num_vec,iter); % (GSF)
    xhat_MCF    = zeros(num_vec,iter); % (Modified C-Filter)
    x_main      = zeros(num_vec,iter);
    %x_main(:,1) = x ; %  (main signal)
    
    %% type of noise
    if flag_noise == 1
        Shot_Noise; 
    else
        Mix_Noise; 
    end
    
    %% one time simulation
    for t = 1 : 1 : iter
        % make the measurement signal
        z = B*x;
        z = z + MeasErrZ(:,t);

        %%  ======================= C_Filter ========================
        xhat2         = CF(xhat2, z, A, B);
        xhat_CF(:, t) = xhat2;
        
        %% ======================= MCC_kf ========================
        [xhat3, P_MCC_KF] = MCCCF(xhat3, P_MCC_KF, z, A, B, Q, R, num_vec);
        xhat_MCC_KF(:, t) = xhat3;
        
        %% ======================= Unscented Kalman Filter (UKF) ========================
        [xhat4, P_UKF] = UKF(F, xhat4, P_UKF, H, z, Q, R);
        xhat_UKF(:, t) = xhat4;
        
        %% ======================= Ensemble Kalman Filter(EnKF) =====================
        [xhat6, X_enkf] = EnKF(X_enkf, z, A, B, Q, R, N_Enkf);
        xhat_EnKF(:, t) = xhat6;

        %% ========================= Adaptive Gaussian Sum Filter (GSF) ==================
        [xhat7, P_GSF] = GSF(xhat7, P_GSF, z, A, B, Q, R, num_vec, num_meas);
        xhat_GSF(:, t) = xhat7;

        %% ===========================  MC_Filter ==============================
%         xhat8      = A * xhat8;
%         innov      = z - B * xhat8;
%         norm_innov = sqrt((innov)'*(innov));
%         sigma      = 1 * norm_innov;
%         K1         = exp(-(norm_innov)^2/(2*sigma^2));
%         Gain       = pinv(eye(num_vec)+ K1 * B'*B )*K1 * B';
%         xhat8      = xhat8 + Gain *(innov);
%         innov      = z - B * xhat8;
%         norm_innov = sqrt((innov)'*(innov));
%         sigma      = 1 * norm_innov;
%         K1         = exp(-(norm_innov)^2/(2*sigma^2));
%         Gain       = pinv(eye(num_vec)+ K1* B'*B )*K1 * B';
%         xhat8      = xhat8 + Gain *(innov);
        xhat8 = MCF(xhat8, z, A, B, Q, R, num_vec);
        
        xhat_MCF(:, t) = xhat8;
        %% Simulate the system.
        x = A*x + MeasErrX(:,t);
        x_main(:, t + 1) = x;
        
    end
    x_main(:,iter) = [];
    
    SE_CF(Numexper,:,:)     = (x_main - xhat_CF).^2;
    SE_MCC_KF(Numexper,:,:) = (x_main - xhat_MCC_KF).^2;
    SE_UKF(Numexper,:,:)    = (x_main - xhat_UKF).^2;
    SE_EnKF(Numexper,:,:)   = (x_main - xhat_EnKF).^2;
    SE_GSF(Numexper,:,:)    = (x_main - xhat_GSF).^2;
    SE_MCF(Numexper,:,:)    = (x_main - xhat_MCF).^2;
end

MSE_CF     = zeros(num_vec,iter);
MSE_MCC_KF = zeros(num_vec,iter);
MSE_UKF    = zeros(num_vec,iter);
MSE_EnKF   = zeros(num_vec,iter);
MSE_GSF    = zeros(num_vec,iter);
MSE_MCF    = zeros(num_vec,iter);

for i = 1 : iter
    MSE_CF(:,i)     = sqrt(mean(SE_CF(:,:,i)))';
    MSE_MCC_KF(:,i) = sqrt(mean(SE_MCC_KF(:,:,i)))';
    MSE_UKF(:,i)    = sqrt(mean(SE_UKF(:,:,i)))';
    MSE_EnKF(:,i)   = sqrt(mean(SE_EnKF(:,:,i)))';
    MSE_GSF(:,i)    = sqrt(mean(SE_GSF(:,:,i)))';
    MSE_MCF(:,i)    = sqrt(mean(SE_MCF(:,:,i)))';
end

MMSE_CF     = mean(MSE_CF,2);
MMSE_MCC_KF = mean(MSE_MCC_KF,2);
MMSE_UKF    = mean(MSE_UKF,2);
MMSE_EnKF   = mean(MSE_EnKF,2);
MMSE_GSF    = mean(MSE_GSF,2);
MMSE_MCF    = mean(MSE_MCF,2);

MMMSE_UKF    = mean(MMSE_UKF);
MMMSE_GSF    = mean(MMSE_GSF);
MMMSE_EnKF   = mean(MMSE_EnKF);
MMMSE_CF     = mean(MMSE_CF);
MMMSE_MCF    = mean(MMSE_MCF);
MMMSE_MCC_KF = mean(MMSE_MCC_KF);
%% Normlized RMSE
max_UKF    = max(MSE_UKF.');
max_CF     = max(MSE_CF.');
max_MCF    = max(MSE_MCF.');
max_MCC_KF = max(MSE_MCC_KF.');
max_EnKF   = max(MSE_EnKF.');
max_GSF    = max(MSE_GSF.');
NRMSE      = zeros(1,6);

j = 1;
for i = 1:4
    NRMSE(1,j) = NRMSE(1,j) + MMSE_UKF(i,1)/ max_UKF(1,i);
end
j = 2;
for i = 1:4
    NRMSE(1,j) =NRMSE(1,j) + MMSE_GSF(i,1)/ max_GSF(1,i);
end
j = 3;
for i = 1:4
    NRMSE(1,j) = NRMSE(1,j) + MMSE_EnKF(i,1)/ max_EnKF(1,i);
end
j = 4;
for i = 1:4
    NRMSE(1,j) = NRMSE(1,j) + MMSE_CF(i,1)/ max_CF(1,i);
end
j = 5;
for i = 1:4
    NRMSE(1,j) = NRMSE(1,j) + MMSE_MCF(i,1)/ max_MCF(1,i);
end
j = 6;
for i = 1:4
    NRMSE(1,j) = NRMSE(1,j) + MMSE_MCC_KF(i,1)/ max_MCC_KF(1,i);
end

%% Plot data.
close all
SetPlotOptions;
t = 0 : T : tf-T;
c = 3;

if flag_noise == 1 % ploting in Shot noise case
    figure,
      xlabel('time'), ylabel('Value of shot noise')

    hold on
    X = index_rand_shot;
    Y = MeasErrZ(1,index_rand_shot(1:end));
    for i=1:num_shot_noise
        bar(X(i)*T,Y(i),1,'b','EdgeColor','b');
    end
    
    figure
    hold on
    plot(t,MSE_CF(c,:) ,'r',t,MSE_MCC_KF(c,:) ,'g');legend('C-Filter','MCC-KF'),xlabel('time'), ylabel('Root mean square error')

else % ploting in mixture guassian noise case
      
    f1 =  @(s) 1/sqrt(2*pi*sqrt(Q_n1(1,1)))*(exp((-(s-mu_n1_x(1,1))*(s-mu_n1_x(1,1)))/(2*sqrt(Q_n1(1,1))))+exp((-(s-mu_n2_x(1,1))*(s-mu_n2_x(1,1)))/(2*sqrt(Q_n2(1,1)))));
    f2 =  @(s) 1/sqrt(2*pi*sqrt(R_n1(1,1)))*(exp((-(s-mu_n1_z(1,1))*(s-mu_n1_z(1,1)))/(2*sqrt(R_n1(1,1))))+exp((-(s-mu_n2_z(1,1))*(s-mu_n2_z(1,1)))/(2*sqrt(R_n2(1,1)))));
    
    figure
    subplot(2,1,1),fplot(f1,[-6,6]);title('pdf of process noise')
    subplot(2,1,2),fplot(f2,[-6,6]);title('pdf of measurment noise')
    
    figure
    hold on
    plot(t,MSE_CF(c,:) ,'r',t,MSE_MCC_KF(c,:) ,'g'); legend('C-Filter','MCC-KF'),xlabel('time'), ylabel('Root mean square error')
end

%% Table result

disp('Mean square error : ');
disp('           x1          x2          x3        x4');
disp(['UKF    : ',num2str(MMSE_UKF.'),'']);
disp(['GSF    : ',num2str(MMSE_GSF.'),'']);
disp(['EnKF   : ',num2str(MMSE_EnKF.'),'']);
disp(['CF     : ',num2str(MMSE_CF.'),'']);
disp(['MCF    : ',num2str(MMSE_MCF.'),'']);
disp(['MCC_KF : ',num2str(MMSE_MCC_KF.'),'']);
