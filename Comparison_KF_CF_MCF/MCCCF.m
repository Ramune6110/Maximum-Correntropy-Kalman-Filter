function [xEst, PEst] = MCCCF(xEst, PEst, z, A, B, Q, R, num_vec)
    xEst       = A * xEst;
    PPred      = A * PEst * A' + Q;
    invers_R   = pinv(R);
    innov      = z - B * xEst;
    norm_innov = sqrt((innov)' * invers_R * (innov));
    sigma      = 1 * norm_innov;
    K          = exp(-(norm_innov^2) /(2 * sigma^2));
    Gain       = pinv(pinv(PPred) + K * B' * invers_R * B) * K * B'*invers_R;
    xEst       = xEst + Gain *(innov);
    PEst       = (eye(num_vec) - Gain*B) * PPred *(eye(num_vec) - Gain * B)' + Gain * R * Gain';
end