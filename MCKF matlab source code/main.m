clear;
close all;
clc

%% type of noise
while true
    flag_noise = input('Please choose the type of noise (1 = Laplace noise, 2 = Gaussian mixture noise, 3 = Gaussian noise): ');
    if flag_noise == 1 || flag_noise == 2 || flag_noise == 3
        break;
    end
end

%% Init
Init_system;

num_of_exprement   = 100; % number of itterations

SE_KF     = zeros(num_of_exprement,n,N);
SE_CF     = zeros(num_of_exprement,n,N);
SE_LRKF   = zeros(num_of_exprement,n,N);

%% 100 Monte Carlo simulation
for Numexper = 1:num_of_exprement
    
    disp(['Simulation # ',num2str(Numexper.'),'/',num2str(num_of_exprement)]);
    
    % Simulation Init Paramter
    Init_parameter;
    
    sV = zeros(n,N);  %true state value
    XEST_KF = zeros(n,N);  %estimated state value
    XEST_CF = zeros(n,N);  %estimated state value
    XEST_LRKF = zeros(n,N);  %estimated state value
    zV = zeros(m,N);  %true measurement value
    bV = zeros(1,N);  %iteration number
    LRKF_ita = zeros(1,N);  %iteration number

    if flag_noise == 1
        Laplace_Noise;
    elseif flag_noise == 2
        Mix_Noise; 
    elseif flag_noise == 3
        Gaussian_Noise;
    end

    %% simulate Kalman FIlter
    for k=1:N

        % Observation
        xTrue = F * xTrue + MeasErrX(:, k); %true state update process
        z = H * xTrue + MeasErrZ(:, k);     %true measurement update process

        % KF
        [xEst_KF, PEst_KF] = KF(xEst_KF, PEst_KF, z, F, H, Q, R, n);
        % MCKF
        [xEst_MCKF, PEst_MCKF, b] = MCKF(F, xEst_MCKF ,PEst_MCKF ,H, z, Q, R); % KF
        % LRKF
        [xEst_LRKF, PEst_LRKF, in, W] = LRKF(xEst_LRKF, PEst_LRKF, z, F, H, Q, R, W);

        % save data
        sV(:,k + 1) = xTrue;           %save true state value
        zV(:,k) = z;               %save true measurement value
        XEST_KF(:,k)   = xEst_KF;            %save estimated state value
        XEST_CF(:,k) = xEst_MCKF;            %save estimated state value
        XEST_LRKF(:,k) = xEst_LRKF;            %save estimated state value
        bV(:,k) = b;               %save iteration number
        LRKF_ita(:,k) = in;               %save iteration number
    end

    sV(:, N) = [];
    
    SE_KF(Numexper,:,:)     = (sV - XEST_KF).^2;
    SE_CF(Numexper,:,:)     = (sV - XEST_CF).^2;
    SE_LRKF(Numexper,:,:)    = (sV - XEST_LRKF).^2;
    
end

MSE_KF     = zeros(n,N);
MSE_CF     = zeros(n,N);
MSE_LRKF   = zeros(n,N);

for i = 1 : N
    MSE_KF(:,i)     = sqrt(mean(SE_KF(:, :, i)))';
    MSE_CF(:,i)     = sqrt(mean(SE_CF(:, :, i)))';
    MSE_LRKF(:,i)   = sqrt(mean(SE_LRKF(:, :, i)))';
end

MMSE_KF     = mean(MSE_KF,2);
MMSE_CF     = mean(MSE_CF,2);
MMSE_LRKF   = mean(MSE_LRKF,2);

%% plot results
figure(1);
plot(1:N, sV(1,:),'b-', 1:N,XEST_CF(1,:),'r--', 1:N,XEST_KF(1,:),'g--', 1:N,XEST_LRKF(1,:),'m--');
legend('true','MCKF', 'KF', 'LRKF');
xlabel('{\ittime step k}','Fontsize',15,'Fontname','Times new roman');
title('state estimate x(1)','Fontsize',15,'Fontname','Times new roman');

figure(2);
plot(1:N, sV(2,:),'b-', 1:N,XEST_CF(2,:),'r--', 1:N,XEST_KF(2,:),'g--', 1:N,XEST_LRKF(2,:),'m--')
legend('true','MCKF', 'KF', 'LRKF');
xlabel('{\ittime step k}','Fontsize',15,'Fontname','Times new roman');
title('state estimate x(2)','Fontsize',15,'Fontname','Times new roman');

f1 =  @(s) 1/sqrt(2*pi*sqrt(Q_n1(1,1)))*(exp((-(s-mu_n1_x(1,1))*(s-mu_n1_x(1,1)))/(2*sqrt(Q_n1(1,1))))+exp((-(s-mu_n2_x(1,1))*(s-mu_n2_x(1,1)))/(2*sqrt(Q_n2(1,1)))));
f2 =  @(s) 1/sqrt(2*pi*sqrt(R_n1(1,1)))*(exp((-(s-mu_n1_z(1,1))*(s-mu_n1_z(1,1)))/(2*sqrt(R_n1(1,1))))+exp((-(s-mu_n2_z(1,1))*(s-mu_n2_z(1,1)))/(2*sqrt(R_n2(1,1)))));

figure(3);
subplot(2,1,1);
fplot(f1,[-6,6]);title('pdf of process noise');
subplot(2,1,2);
fplot(f2,[-6,6]);title('pdf of measurment noise');

figure(4);
plot(1:N,bV(1,:),'b--');
plot(1:N,LRKF_ita(1,:),'m--');

figure(5);
for i = 1:n
        figure(i + 4);
        hold on;
        plot(1:N,MSE_KF(i, :) ,'g','LineWidth',1.0);
        plot(1:N,MSE_CF(i, :) ,'r','LineWidth',1.0);
        plot(1:N,MSE_LRKF(i, :) ,'m','LineWidth',1.0);
        grid on;
        legend('KF','CF','LRKF','Location','best');
        xlabel('time');
        ylabel('Root mean square error');
 end

%% Table result

disp('Mean square error : ');
disp('           x1          x2');
disp(['KF     : ',num2str(MMSE_KF.'),'']);
disp(['CF     : ',num2str(MMSE_CF.'),'']);
disp(['LRKF   : ',num2str(MMSE_LRKF.'),'']);
