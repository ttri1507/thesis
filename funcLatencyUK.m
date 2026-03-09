function t_uk = funcLatencyUK(d_uk, h_uk, p_uk, H_other, P_other, Q, sigma, BW)
v_light = 3e8; % The velocity of light

SINR_uk = ( abs(h_uk)^2*p_uk )/( sigma + sum(abs(H_other).^2.*P_other) );
R_uk = BW/log(2)*log(1+ SINR_uk);

t_uk = 2*d_uk/v_light + Q/R_uk;
end