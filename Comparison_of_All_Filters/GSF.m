function [xEst, PEst] = GSF(xEst, PEst, z, A, B, Q, R, num_vec, num_meas)
    % initialization for GSF
    Ng      = 2; 
    miu     = ones(num_meas,Ng); 
    lambda  = 100;
    epsilon = 0.5; 
    a_i_1   = epsilon;
    a_i_2   = 1-epsilon;
    m       = num_meas;
    P_zi    = zeros(num_meas,num_meas,Ng); % weights for GSF
    
    xPred       = A * xEst;
    PPred       = A * PEst * A' + Q;
    zhat_i      = repmat(B * xPred,1,Ng) + miu;
    P_zi(:,:,1) = B * PPred * B' + R;
    P_zi(:,:,2) = B * PPred * B' + lambda * R;
    c_k         = ((2*pi)^-m )/(det(P_zi(:,:,1))) * exp(-0.5 * (z-zhat_i(:,1))' /(P_zi(:,:,1))* (z-zhat_i(:,1))) * a_i_1 + ((2*pi)^-m )/(det(P_zi(:,:,2))) * exp(-0.5 * (z-zhat_i(:,2))'/(P_zi(:,:,2))* (z-zhat_i(:,2))) * a_i_2 ;
    w_i(:,1)    = ((2*pi)^-m /(det(P_zi(:,:,1))) * exp(-0.5 * (z-zhat_i(:,1))' /(P_zi(:,:,1))* (z-zhat_i(:,1))) * a_i_1) / c_k ;
    w_i(:,2)    = ((2*pi)^-m /(det(P_zi(:,:,2))) * exp(-0.5 * (z-zhat_i(:,2))' /(P_zi(:,:,2))* (z-zhat_i(:,2))) * a_i_2) / c_k ;
    zhat        = w_i(:,1) * zhat_i(:,1) + w_i(:,2) * zhat_i(:,2);
    Pz          = (P_zi(:,:,1) + (zhat-zhat_i(:,1))*(zhat-zhat_i(:,1))') * w_i(:,1) + (P_zi(:,:,2) + (zhat-zhat_i(:,2))*(zhat-zhat_i(:,2))') * w_i(:,2);
    K           = (PPred * B' )/(Pz);
    xEst        = xPred + K * (z - zhat);
    PEst        = (eye(num_vec)-K*B)*PPred;
end