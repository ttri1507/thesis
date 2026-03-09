function A_out = FUNC_GameTheory_Clustering(MODEL, A, B, numSAT, numUAV,...
    numUE, numUAVeachSAT, numUEeachUAV, N_U, PmaxSAT, PmaxUAV, Area, Q, Noise_var, BW)
numLoop = 0;
A_old = A;
while(true)
    numLoop = numLoop + 1; %display(numLoop);
    A_min = A_old;
    for k = 1:numUE
        f = MODEL.RP(k);
        s = ceil(k/(numUAVeachSAT*numUEeachUAV));
        Us = numUAVeachSAT*(s-1) + 1:numUAVeachSAT*s;
        UT_min = 1e6;
        for id_u = 1:numUAVeachSAT
            u = Us(id_u);
            A_temp = A_min; A_temp(:,k) = 0; A_temp(u,k) = 1;
            if(all(sum(A_temp, 2) <= N_U)) % Check whether there is any UAV being overloading
                % Update temporary Power
                [P_sat_temp, P_uav_uav_temp, P_uav_ue_temp, ASU, AUU] = FUNC_PA(MODEL, A_temp, B, ...
                    numSAT, numUAV, numUE, numUAVeachSAT, numUEeachUAV, PmaxSAT, PmaxUAV);
                [UT, case_tk] = FUNCLatency(MODEL, k, f, A_temp, B, P_sat_temp, P_uav_uav_temp, P_uav_ue_temp,...
                    Area, Q, Noise_var, BW, ASU, AUU);
            else
                UT = 1e6; % a very large number
            end
            % Find best response (minimum)
            if(UT < UT_min)
                UT_min = UT;
                A_min = A_temp;
            end
        end
    end
    if((sum(sum(A_min==A_old))>0.99*numUAV*numUE)||(numLoop==100))
        A_out = A_min;
        break
    else
        A_old = A_min;
    end
end
end