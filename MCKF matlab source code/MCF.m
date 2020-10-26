function [xEst, b] = MCF(xEst, z, A, B, Q, R, num_vec)
    xEst       = A * xEst;
    innov      = z - B * xEst;
    norm_innov = sqrt((innov)'*(innov));
    sigma      = 1 * norm_innov;
    K          = exp(-(norm_innov)^2/(2 * sigma^2));
    Gain       = pinv(eye(num_vec) + K * B'*B ) * K * B';
    xEst       = xEst + Gain *(innov);
    innov      = z - B * xEst;
    norm_innov = sqrt((innov)'*(innov));
    sigma      = 1 * norm_innov;
    K          = exp(-(norm_innov)^2/(2 * sigma^2));
    Gain       = pinv(eye(num_vec) + K * B' * B ) * K * B';
    xEst       = xEst + Gain * (innov);
    b = 0;
end