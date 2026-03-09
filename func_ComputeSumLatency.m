function [t_total, T] = func_ComputeSumLatency(MODEL, A, B,...
              Area, Q, Noise_var, BW, numUAVeachSAT, numUEeachUAV, numSAT, numUAV, numUE, PmaxSAT, PmaxUAV)
T = [];
for k = 1:numUE
    f = MODEL.RP(k);
    [P_sat_temp, P_uav_uav_temp, P_uav_ue_temp, ASU, AUU] = FUNC_PA(MODEL, A, B, ...
            numSAT, numUAV, numUE, numUAVeachSAT, numUEeachUAV, PmaxSAT, PmaxUAV);
    [tk, case_tk] = FUNCLatency(MODEL, k, f, A, B, P_sat_temp, P_uav_uav_temp, P_uav_ue_temp,...
            Area, Q, Noise_var, BW, ASU, AUU);
    T = [T, [tk;case_tk]];
end
t_total = sum(T(1,:));
end