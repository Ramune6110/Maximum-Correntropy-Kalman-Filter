clear
close all
clc

%% call system dynamic
Init_System;

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

%% 100 Monte Carlo simulation
for Numexper = 1:num_of_exprement
    
    disp(['Simulation # ',num2str(Numexper.'),'/',num2str(num_of_exprement)]);
    
    %% Simulation Init Paramter
    Init_Parameter;
    
    %% type of noise
    if flag_noise == 1
        Shot_Noise; 
    else
        Mix_Noise; 
    end
    
    %% one time simulation
    for t = 1 : 1 : iter
       %%  ======================= Observation =======================
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
        xhat8          = MCF(xhat8, z, A, B, Q, R, num_vec);
        xhat_MCF(:, t) = xhat8;

       %% Simulate the system.
        x = A * x + MeasErrX(:, t);
        x_main(:, t + 1) = x;
    end
    
    x_main(:, iter) = [];
    
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

%% Plot data.
close all
SetPlotOptions;
t = 0 : T : tf-T;
c = 3;

if flag_noise == 1 % ploting in Shot noise case
    figure(1);
    xlabel('time');
    ylabel('Value of shot noise');
    hold on
    X = index_rand_shot;
    Y = MeasErrZ(1,index_rand_shot(1:end));
    for i=1:num_shot_noise
        bar(X(i)*T,Y(i),1,'b','EdgeColor','b');
    end
    
    figure(2);
    hold on
    plot(t,MSE_CF(c,:) ,'r','LineWidth',1.0);
    plot(t,MSE_MCF(c,:) ,'b','LineWidth',1.0);
    plot(t,MSE_MCC_KF(c,:) ,'g','LineWidth',1.0);
    plot(t,MSE_GSF(c,:) ,'m','LineWidth',1.0);
    plot(t,MSE_UKF(c,:) ,'k','LineWidth',1.0);
    plot(t,MSE_EnKF(c,:) ,'c','LineWidth',1.0);
    grid on;
    legend('CF','MCF','MCCKF', 'GSF', 'UKF', 'EnKF','Location','best');
    xlabel('time');
    ylabel('Root mean square error');

else % ploting in mixture guassian noise case
      
    f1 =  @(s) 1/sqrt(2*pi*sqrt(Q_n1(1,1)))*(exp((-(s-mu_n1_x(1,1))*(s-mu_n1_x(1,1)))/(2*sqrt(Q_n1(1,1))))+exp((-(s-mu_n2_x(1,1))*(s-mu_n2_x(1,1)))/(2*sqrt(Q_n2(1,1)))));
    f2 =  @(s) 1/sqrt(2*pi*sqrt(R_n1(1,1)))*(exp((-(s-mu_n1_z(1,1))*(s-mu_n1_z(1,1)))/(2*sqrt(R_n1(1,1))))+exp((-(s-mu_n2_z(1,1))*(s-mu_n2_z(1,1)))/(2*sqrt(R_n2(1,1)))));
    
    figure(1);
    subplot(2,1,1);
    fplot(f1,[-6,6]);title('pdf of process noise');
    subplot(2,1,2);
    fplot(f2,[-6,6]);title('pdf of measurment noise');
    
    figure(2);
    hold on
    plot(t,MSE_CF(c,:) ,'r','LineWidth',1.0);
    plot(t,MSE_MCF(c,:) ,'b','LineWidth',1.0);
    plot(t,MSE_MCC_KF(c,:) ,'g','LineWidth',1.0);
    plot(t,MSE_GSF(c,:) ,'m','LineWidth',1.0);
    plot(t,MSE_UKF(c,:) ,'k','LineWidth',1.0);
    plot(t,MSE_EnKF(c,:) ,'c','LineWidth',1.0);
    grid on;
    legend('CF','MCF','MCCKF', 'GSF', 'UKF', 'EnKF','Location','best');
    xlabel('time');
    ylabel('Root mean square error');
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
