function xEst = CF(xEst, z, A, B)
    xEst       = A * xEst;
    innov      = z - B * xEst;
    norm_innov = sqrt(innov' * innov);
    sigma      = 1 * norm_innov;
    K          = exp(-(norm_innov^2) /(2 * sigma^2));
    xEst       = xEst +  K * B' * (innov);
end