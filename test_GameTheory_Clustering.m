clc
clear
load('model.mat');

% [P_sat, P_uav_uav, P_uav_ue] = FUNC_PA(MODEL, A, B, ...
% numSAT, numUAV, numUE, numUAVeachSAT, numUEeachUAV, PmaxSAT, PmaxUAV);

for k = 1:numUE
    f = RP(k);
    s = ceil(k/(numUAVeachSAT*numUEeachUAV));
    Us = numUAVeachSAT*(s-1) + 1:numUAVeachSAT*s;
    UT_min = 1e6; A_min = A;
    for id_u = 1:numUAVeachSAT
        u = Us(id_u);
        A_temp = A; A_temp(:,k) = 0; A_temp(u,k) = 1;
        if(all(sum(A_temp, 2) <= N_U)) % Check whether there is any UAV being overloading
            % Update temporary Power
            [P_sat_temp, P_uav_uav_temp, P_uav_ue_temp] = FUNC_PA(MODEL, A_temp, B, ...
    numSAT, numUAV, numUE, numUAVeachSAT, numUEeachUAV, PmaxSAT, PmaxUAV);
            UT = FUNCLatency(MODEL, k, f, A_temp, B, P_sat_temp, P_uav_uav_temp, P_uav_ue_temp,...
                  Area, Q, Noise_var, BW);
        else
            UT = 1e6; % a very large number
        end
        % Find best response (minimum)
        if(UT < UT_min)
            UT_min = UT;
            A_min = A_temp;
        end
    end
    A = A_min;
end

