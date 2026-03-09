function [Fitness, Fit_min, bestEle] = func_GA_Evaluate(Parent, P, MODEL, A, B_full,...
              Area, Q, Noise_var, BW, s, numUAVeachSAT, numUEeachUAV, numM, numSAT, numUAV, numUE, PmaxSAT, PmaxUAV)
numUE_s = numUAVeachSAT*numUEeachUAV;
UEs = (numUE_s*(s-1) + 1):(numUE_s*s);
Us = numUAVeachSAT*(s-1) + 1:numUAVeachSAT*s;
Fitness = zeros(1,P);
for par = 1:P
    T = [];
    for idk = 1:size(UEs,2)
        k = UEs(idk);
        f = MODEL.RP(k);
        B = Parent(:,:,par);
        if(any(sum(B,2) > numM)) % any of element breaks the constraint
            tk =1e6; case_tk = 0;
        else
            % Update temporary Power
            B_full_temp = B_full; B_full_temp(Us,:) = B;
            [P_sat_temp, P_uav_uav_temp, P_uav_ue_temp, ASU, AUU] = FUNC_PA(MODEL, A, B_full_temp, ...
        numSAT, numUAV, numUE, numUAVeachSAT, numUEeachUAV, PmaxSAT, PmaxUAV);
            [tk, case_tk] = FUNCLatency(MODEL, k, f, A, B_full_temp, P_sat_temp, P_uav_uav_temp, P_uav_ue_temp,...
                      Area, Q, Noise_var, BW, ASU, AUU);
        end
        T = [T, [tk;case_tk]];
        Fitness(par) = Fitness(par) + tk;
    end
end
[Fit_min, bestEle] = min(Fitness);
end