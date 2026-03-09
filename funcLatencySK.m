function t_sk = funcLatencySK(d_uk, h_uk, p_uk, H_other, P_other, ...
    d_su, h_su, precoding_su, p_su, Q, sigma, BW, n_request)
v_light = 3e8; % The velocity of light

SINR_uk = ( abs(h_uk)^2*p_uk )/( sigma + sum(abs(H_other).^2.*P_other) );
R_uk = BW/log(2)*log(1+ SINR_uk);

R_su = BW/log(2)*log(1+ abs(squeeze(h_su)'*squeeze(precoding_su))^2*p_su/sigma);
t_sk = 2*(d_uk + d_su)/v_light + Q/R_uk + Q*n_request/R_su;
end