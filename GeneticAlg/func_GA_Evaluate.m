function [Fitness, Fit_min, bestEle] = func_GA_Evaluate(Parent, P, MODEL, A, P_sat, P_uav_uav, P_uav_ue,...
              Area, Q, sigma, BW, s, numUAVeachSAT, numUEeachUAV, numM)
numUE_s = numUAVeachSAT*numUEeachUAV;
UEs = (numUE_s*(s-1) + 1):(numUE_s*s);
Fitness = zeros(1,P);
for par = 1:P
    for idk = 1:size(UEs,2)
        k = UEs(idk);
        f = MODEL.RP(k);
        B = Parent(:,:,par);
        if(any(sum(B,2) > numM)) % any of element breaks the constraint
            tk =1e6;
        else
            tk = FUNCLatency(MODEL, k, f, A, B, P_sat, P_uav_uav, P_uav_ue,...
                      Area, Q, sigma, BW);
        end
        Fitness(par) = Fitness(par) + tk;
    end
end
[Fit_min, bestEle] = min(Fitness);
end