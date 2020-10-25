function [xEst, PEst] = KF(xEst, PEst, z, A, B, Q, R, num_vec)
    xPred      = A * xEst;
    PPred      = A * PEst * A' + Q;
    invers_R   = pinv(R);
    innov      = z - B * xPred;
    Gain       = pinv(pinv(PPred) + B' * invers_R * B) * B'*invers_R;
    xEst       = xPred + Gain *(innov);
    PEst       = (eye(num_vec) - Gain*B) * PPred *(eye(num_vec) - Gain * B)' + Gain * R * Gain';
end