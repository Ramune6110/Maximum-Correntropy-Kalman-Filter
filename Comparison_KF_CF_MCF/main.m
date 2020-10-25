clear
close all
clc

%% Choice the type of System
while true
    model = input('Please choose the type of System (1 = Example1, 2 = Example2, 3 = Example3): ');
    switch model
        case 1
            model01;
            break;
        case 2
            model02;
            break;
        case 3
            model03;
            break;
        otherwise
            disp('No System');
    end
end

%% Choice the type of noise
while true
    flag_noise = input('Please choose the type of noise (1 = shot noise, 2 = Gaussian mixture noise, 3 = Gaussian noise): ');
    if flag_noise == 1 || flag_noise == 2 || flag_noise == 3
        break;
    end
end

%% Square Error
SE_KF     = zeros(num_of_exprement,num_vec,iter);
SE_CF     = zeros(num_of_exprement,num_vec,iter);
SE_MCF    = zeros(num_of_exprement,num_vec,iter);
SE_MCC_KF = zeros(num_of_exprement,num_vec,iter);

%% 100 Monte Carlo simulation
for Numexper = 1:num_of_exprement
    
    disp(['Simulation # ',num2str(Numexper.'),'/',num2str(num_of_exprement)]);
    
    %% Simulation Init Paramter
    Init_Parameter;
    
    %% define some arrays to store the main signal and the esitimation signals
    xhat_KF     = zeros(num_vec,iter); %(C-Filter)
    xhat_CF     = zeros(num_vec,iter); %(C-Filter)
    xhat_MCF    = zeros(num_vec,iter); % (Modified C-Filter)
    xhat_MCC_KF = zeros(num_vec,iter); % (Modified C-Filter)
    x_main      = zeros(num_vec,iter);

    %% type of noise
    if flag_noise == 1
        Shot_Noise; 
    elseif flag_noise == 2
        Mix_Noise; 
    elseif flag_noise == 3
        Gaussian_Noise;
    end
    
    %% one time simulation
    for t = 1 : 1 : iter
       %%  ======================= Observation =======================
        z = B*x;
        z = z + MeasErrZ(:,t);

       %%  ======================= C_Filter ========================
        xhat2 = CF(xhat2, z, A, B);
        xhat_CF(:, t) = xhat2;
        
        %% ===========================  KF ==============================
        [xhat3, P_KF] = KF(xhat3, P_KF, z, A, B, Q, R, num_vec);
        xhat_KF(:, t) = xhat3;
        
        %% ======================= MCC_kf ========================
        [xhat4, P_MCC_KF] = MCCCF(xhat4, P_MCC_KF, z, A, B, Q, R, num_vec);
        xhat_MCC_KF(:, t) = xhat4;
        
       %% ===========================  MC_Filter ==============================
        xhat8 = MCF(xhat8, z, A, B, Q, R, num_vec);
        xhat_MCF(:, t) = xhat8;
        
       %% Simulate the system.
        x = A * x + MeasErrX(:, t);
        x_main(:, t + 1) = x;
    end
    
    x_main(:, iter) = [];
    
    SE_KF(Numexper,:,:)     = (x_main - xhat_KF).^2;
    SE_CF(Numexper,:,:)     = (x_main - xhat_CF).^2;
    SE_MCF(Numexper,:,:)    = (x_main - xhat_MCF).^2;
    SE_MCC_KF(Numexper,:,:) = (x_main - xhat_MCC_KF).^2;
    
end

%% Mean Square Error
MSE_KF     = zeros(num_vec,iter);
MSE_CF     = zeros(num_vec,iter);
MSE_MCF    = zeros(num_vec,iter);
MSE_MCC_KF = zeros(num_vec,iter);

for i = 1 : iter
    MSE_KF(:,i)     = sqrt(mean(SE_KF(:,:,i)))';
    MSE_CF(:,i)     = sqrt(mean(SE_CF(:,:,i)))';
    MSE_MCF(:,i)    = sqrt(mean(SE_MCF(:,:,i)))';
    MSE_MCC_KF(:,i) = sqrt(mean(SE_MCC_KF(:,:,i)))';
end

MMSE_KF     = mean(MSE_KF,2);
MMSE_CF     = mean(MSE_CF,2);
MMSE_MCF    = mean(MSE_MCF,2);
MMSE_MCC_KF = mean(MSE_MCC_KF,2);

MMMSE_KF     = mean(MMSE_KF);
MMMSE_CF     = mean(MMSE_CF);
MMMSE_MCF    = mean(MMSE_MCF);
MMMSE_MCC_KF = mean(MMSE_MCC_KF);

%% Normlized RMSE
max_KF     = max(MSE_KF.');
max_CF     = max(MSE_CF.');
max_MCF    = max(MSE_MCF.');
max_MCC_KF = max(MSE_MCC_KF.');

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
    
    for i = 1:num_vec
        figure(i + 1);
        hold on;
        plot(t,MSE_KF(i, :) ,'g','LineWidth',1.0);
        plot(t,MSE_CF(i, :) ,'r','LineWidth',1.0);
        plot(t,MSE_MCF(i, :) ,'b','LineWidth',1.0);
        plot(t,MSE_MCC_KF(i, :) ,'m','LineWidth',1.0);
        grid on;
        legend('KF','CF','MCF','MCCKF','Location','best');
        xlabel('time');
        ylabel('Root mean square error');
    end

else % ploting in mixture guassian noise case
      
    f1 =  @(s) 1/sqrt(2*pi*sqrt(Q_n1(1,1)))*(exp((-(s-mu_n1_x(1,1))*(s-mu_n1_x(1,1)))/(2*sqrt(Q_n1(1,1))))+exp((-(s-mu_n2_x(1,1))*(s-mu_n2_x(1,1)))/(2*sqrt(Q_n2(1,1)))));
    f2 =  @(s) 1/sqrt(2*pi*sqrt(R_n1(1,1)))*(exp((-(s-mu_n1_z(1,1))*(s-mu_n1_z(1,1)))/(2*sqrt(R_n1(1,1))))+exp((-(s-mu_n2_z(1,1))*(s-mu_n2_z(1,1)))/(2*sqrt(R_n2(1,1)))));
    
    figure(1);
    subplot(2,1,1);
    fplot(f1,[-6,6]);title('pdf of process noise');
    subplot(2,1,2);
    fplot(f2,[-6,6]);title('pdf of measurment noise');
    
    for i = 1:num_vec
        figure(i + 1);
        hold on;
        plot(t,MSE_KF(i, :) ,'g','LineWidth',1.0);
        plot(t,MSE_CF(i, :) ,'r','LineWidth',1.0);
        plot(t,MSE_MCF(i, :) ,'b','LineWidth',1.0);
        plot(t,MSE_MCC_KF(i, :) ,'m','LineWidth',1.0);
        grid on;
        legend('KF','CF','MCF','MCCKF','Location','best');
        xlabel('time');
        ylabel('Root mean square error');
    end
end

%% Table result

disp('Mean square error : ');
disp('           x1          x2          x3        x4');
disp(['KF     : ',num2str(MMSE_KF.'),'']);
disp(['CF     : ',num2str(MMSE_CF.'),'']);
disp(['MCF    : ',num2str(MMSE_MCF.'),'']);
disp(['MCC_KF : ',num2str(MMSE_MCC_KF.'),'']);