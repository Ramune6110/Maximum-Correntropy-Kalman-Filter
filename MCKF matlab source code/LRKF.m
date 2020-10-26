function [xEst, PEst, K, b, W] = LRKF(xEst, PEst, z, A, B, Q, R, W)
    % Predict
    xPred = A * xEst;
    PPred = A * PEst * A' + Q;
    
    b = 0; %initial iteration number
    bias = 100;
    epsilon = 1e-6;
    
    % Now iterate
    while (norm((xEst - xPred)/norm(xPred))<=bias && (b<=7))
        % Update
        R_overline = (sqrt(2) / 2) * sqrt(R) * W * sqrt(R);
        K    = (PPred * B') / (B * PPred * B' + R_overline);
        xEst = xPred + K * (z - B * xPred);
        W = abs(sqrt(R) * (z - B * xPred)) + epsilon;
        % loop counter
        b   = b + 1;
    end
    
    PEst = (eye(size(xEst,1)) - K * B) * PPred;
    
end