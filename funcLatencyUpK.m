function t_upk = funcLatencyUpK(d_uk, h_uk, p_uk, H_other, P_other, ...
    d_upu, h_upu, p_upu, Q, sigma, BW, n_request)
v_light = 3e8; % The velocity of light

SINR_uk = ( abs(h_uk)^2*p_uk )/( sigma + sum(abs(H_other).^2.*P_other) );
R_uk = BW/log(2)*log(1+ SINR_uk);

R_upu = BW/log(2)*log(1+ abs(h_upu)^2*p_upu/sigma);
t_upk = 2*(d_uk + d_upu)/v_light + Q/R_uk + Q*n_request/R_upu;
end